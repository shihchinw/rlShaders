// Disney Principled BRDF based on
// [1]  Physically Based Shading at Disney. Brent Burley. SIGGRAPH 2012 course notes.
//      http://blog.selfshadow.com/publications/s2012-shading-course/
//      https://github.com/wdas/brdf/blob/master/src/brdfs/disney.brdf
// [2]  Microfacet Models for Refraction through Rough Surfaces. Bruce Walter et al. EGSR'07
//      http://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html
// [3] 	Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals.
//      Eric Heitz and Eugene d'Eon. EGSR'14
//      https://hal.inria.fr/hal-00996995/en
//      https://hal.inria.fr/hal-00996995v2/file/supplemental1.pdf
#include <vector>
#include <set>
#include <cassert>

#include <ai.h>
//#include "rlUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(DisneyMethod);

inline float colorToLuminance(const AtRGB &col)
{
    return col.r * 0.212671f + col.g * 0.715160f + col.b * 0.072169f;
}

enum    DisneyParams
{
    p_base_color,
    p_subsurface,
    p_metallic,
    p_Ks,
    p_specular_tint,
    p_roughness,
    p_anisotropic,
    p_sheen,
    p_sheen_tint,
    p_clearcoat,
    p_clearcoat_gloss,
    p_opacity,

    p_indirect_diffuse,
    p_indirect_specular,

    p_aov_direct_diffuse,
    p_aov_direct_specular,
    p_aov_indirect_diffuse,
    p_aov_indirect_specular,
};

namespace
{

struct  ShaderData
{
    ShaderData() : mDiffuseSampler(nullptr), mGlossySampler(nullptr)
    {
    }

    virtual ~ShaderData()
    {
        AiSamplerDestroy(mDiffuseSampler);
    }

    void    update()
    {
        AtNode *options = AiUniverseGetOptions();
        mMaxRayDepth = AiNodeGetInt(options, "GI_total_depth");
        mDiffuseDepth = AiNodeGetInt(options, "GI_diffuse_depth");
        mGlossyDepth = AiNodeGetInt(options, "GI_glossy_depth");

        auto diffuseSampleNum = AiNodeGetInt(options, "GI_diffuse_samples");
        AiSamplerDestroy(mDiffuseSampler);
        mDiffuseSampler = AiSampler(diffuseSampleNum, 2);

        auto glossySampleNum = AiNodeGetInt(options, "GI_glossy_samples");
        AiSamplerDestroy(mGlossySampler);
        mGlossySampler = AiSampler(glossySampleNum, 2);
    }

    bool    shouldTraceDiffuse(AtShaderGlobals *sg) const
    {
        return (sg->Rr_diff < mDiffuseDepth && sg->Rr < mMaxRayDepth);
    }

    bool    shouldTraceGlossy(AtShaderGlobals *sg) const
    {
        return (sg->Rr_gloss < mGlossyDepth && sg->Rr < mMaxRayDepth);
    }

    AtSamplerIterator *getDiffuseSamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mDiffuseSampler, sg);
    }

    AtSamplerIterator *getGlossySamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mGlossySampler, sg);
    }

private:
    AtSampler   *mDiffuseSampler;
    AtSampler   *mGlossySampler;
    int         mMaxRayDepth;
    int         mDiffuseDepth;
    int         mGlossyDepth;
};

}

class   DisneySampler
{
public:
    //! Return sample direction
    static AtVector evalSample(const void *brdf, float rx, float ry)
    {
        const auto *disneyBrdf = static_cast<const DisneySampler*>(brdf);
        if (disneyBrdf->mSampleType == AI_RAY_DIFFUSE) {
            return disneyBrdf->sampleDiffuseDirection(rx, ry);
        }

        return disneyBrdf->sampleSpecularDirection(rx, ry);
    }

    //! Return the reflectance
    static AtColor  evalBrdf(const void *brdf, const AtVector *indir)
    {
        const auto L = *indir;

        if (AiV3IsZero(L)) {
            // TODO: Check how did this happen...
            return AI_RGB_BLACK;
        }

        const auto *disneyBrdf = static_cast<const DisneySampler*>(brdf);
        auto NdotL = AiV3Dot(disneyBrdf->mAxisN, L);

        if (disneyBrdf->mSampleType == AI_RAY_DIFFUSE) {
            return disneyBrdf->evalDiffuse(L) * NdotL;
        }

        return disneyBrdf->evalSpecular(L) * NdotL;
    }

    static float evalPdf(const void *brdf, const AtVector *indir)
    {
        if (AiV3IsZero(*indir)) {
            // TODO: Check how did this happen...
            return 0.0f;
        }

        const auto *disneyBrdf = static_cast<const DisneySampler*>(brdf);
        if (disneyBrdf->mSampleType == AI_RAY_DIFFUSE) {
            return disneyBrdf->evalDiffusePdf(*indir);
        }

        return disneyBrdf->evalSpecularPdf(*indir);
    }

public:
    DisneySampler(AtNode *node, AtShaderGlobals *sg)
    {
        mBaseColor = AiShaderEvalParamRGB(p_base_color);
        mRoughness = AiShaderEvalParamFlt(p_roughness);
        mSubsurface = AiShaderEvalParamFlt(p_subsurface);

        // Specular represents the fresnel relectance at normal direction.
        // Remap the range to [0, 0.08] for IOR [1.0, 1.8]
        mSpecular = AiShaderEvalParamFlt(p_Ks) * 0.08f;
        mSpecularTint = AiShaderEvalParamFlt(p_specular_tint);
        mMetallic = AiShaderEvalParamFlt(p_metallic);
        mSheen = AiShaderEvalParamFlt(p_sheen);
        mSheenTint = AiShaderEvalParamFlt(p_sheen_tint);
        mAnisotropic = AiShaderEvalParamFlt(p_anisotropic);
        mClearcoat = AiShaderEvalParamFlt(p_clearcoat) * 0.25f;
        mClearcoatGloss = AiShaderEvalParamFlt(p_clearcoat_gloss);

        mViewDir = -sg->Rd;
        mAxisN = sg->Nf;
        AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);

        //float aniso = fabs(mAnisotropic - 0.5f) * 2.0f;
        float aspect = sqrt(1.0f - mAnisotropic * 0.9f);
        mAlphaX = MAX(1e-4f, SQR(mRoughness) / aspect);
        mAlphaY = MAX(1e-4f, SQR(mRoughness) * aspect);

        mSpecularRoughness = SQR(mRoughness);

        // Compute specularF0 (i.e. reflectance at normal incidence) whiich is used for
        // the materials whose IOR is in [1.0, 1.8].
        float luminance = colorToLuminance(mBaseColor);
        AtRGB tintColor = luminance > 0.0f ? mBaseColor / luminance : AI_RGB_WHITE;
        AtRGB metallicColor = mSpecular * LERP(mSpecularTint, AI_RGB_WHITE, tintColor);
        mSpecularF0 = LERP(mMetallic, metallicColor, mBaseColor);

        mSheenColor = LERP(mSheenTint, AI_RGB_WHITE, tintColor) * mSheen;
        mSampleFromVisibleNormal = true;
    }

    inline void setSampleType(AtUInt16 type)
    {
        mSampleType = type;
    }

    AtColor     evalDiffuse(const AtVector &L) const
    {
        float LdotN = AiV3Dot(L, mAxisN);
        float VdotN = AiV3Dot(mViewDir, mAxisN);

        if (LdotN < AI_EPSILON || VdotN < AI_EPSILON) {
            return AI_RGB_BLACK;
        }

        const auto H = AiV3Normalize(L + mViewDir);
        float LdotH = AiV3Dot(L, H);
        float NdotH = AiV3Dot(mViewDir, H);

        if (NdotH < AI_EPSILON || LdotH < AI_EPSILON) {
            return AI_RGB_BLACK;
        }

        float LdotH2 = SQR(LdotH);

        // Diffuse
        float FL = powf(CLAMP(1.0f - LdotN, 0.0f, 1.0f), 5.0f);
        float FV = powf(CLAMP(1.0f - VdotN, 0.0f, 1.0f), 5.0f);
        assert(FL >= 0.0f);
        assert(FV >= 0.0f);

        float F90 = 0.5f + 2.0f * mRoughness * LdotH2;
        float diffuseFactor = LERP(FL, 1.0f, F90) * LERP(FV, 1.0f, F90);

        // Subsurface
        // Hanrahan-Krueger brdf approximation of isotropic bssrdf
        // 1.25 scale is used to (roughly) preserve albedo
        float Fss90 = mRoughness * LdotH2;
        float Fss = LERP(FL, 1.0f, Fss90) * LERP(FV, 1.0f, Fss90);
        float ssFactor = 1.25f * (Fss * (1.0f / (LdotN + VdotN) - 0.5f) + 0.5f);

        auto diffuse = mBaseColor * AI_ONEOVERPI * LERP(mSubsurface, diffuseFactor, ssFactor);
        return diffuse * (1.0f - mMetallic);
    }

    //! Integrate indirect diffuse reflectance.
    //! Use importance smapling according to cosTheta.
    AtColor     integrateDiffuse(AtShaderGlobals *sg, const ShaderData *data)
    {
        setSampleType(AI_RAY_DIFFUSE);
        return AiBRDFIntegrate(sg, this, evalSample, evalBrdf, evalPdf, AI_RAY_DIFFUSE);

        AtRay       ray;
        AtScrSample scrs;
        AiMakeRay(&ray, AI_RAY_DIFFUSE, &sg->P, nullptr, AI_BIG, sg);

        AtSamplerIterator *sampleIter = data->getDiffuseSamplerIter(sg);
        float   samples[2];
        AtColor result = AI_RGB_BLACK;

        while (AiSamplerGetSample(sampleIter, samples)) {
            ray.dir = sampleDiffuseDirection(samples[0], samples[1]);

            if (AiTrace(&ray, &scrs)) {
                result += scrs.color * evalDiffuse(ray.dir);
            }
        }

        return result * AiSamplerGetSampleInvCount(sampleIter) * AI_PI;
    }

    AtColor     evalDiffuseLightSample(AtShaderGlobals *sg)
    {
        setSampleType(AI_RAY_DIFFUSE);
        return AiEvaluateLightSample(sg, this, evalSample, evalBrdf, evalPdf);
    }

    //! Evaluate MIS with light and glossy lobes.
    AtColor     evalSpecularLightSample(AtShaderGlobals *sg)
    {
        setSampleType(AI_RAY_GLOSSY);
        return AiEvaluateLightSample(sg, this, evalSample, evalBrdf, evalPdf);
    }

    AtColor     integrateGlossy(AtShaderGlobals *sg)
    {
        setSampleType(AI_RAY_GLOSSY);
        return AiBRDFIntegrate(sg, this, evalSample, evalBrdf, evalPdf, AI_RAY_GLOSSY);
    }

    //! Another way to integrate the indirect glossy.
    //! @note: This is used to study what is the behavior of AiBRDFIntegrate.
    AtColor     integrateGlossy(AtShaderGlobals *sg, const ShaderData *data)
    {
        setSampleType(AI_RAY_GLOSSY);

        AtRay       ray;
        AtScrSample scrs;
        AiMakeRay(&ray, AI_RAY_GLOSSY, &sg->P, nullptr, AI_BIG, sg);

        AtSamplerIterator *sampleIter = data->getGlossySamplerIter(sg);
        float   samples[2];
        AtColor result = AI_RGB_BLACK;

        while (AiSamplerGetSample(sampleIter, samples)) {
            ray.dir = evalSample(this, samples[0], samples[1]);

            if (!AiTrace(&ray, &scrs)) {
                //AiTraceBackground(&ray, &scrs);
                continue;
            }

            auto pdf = evalPdf(this, &ray.dir);
            if (pdf > AI_EPSILON) {
                result += scrs.color * evalBrdf(this, &ray.dir) / pdf;
            }
        }

        return result * AiSamplerGetSampleInvCount(sampleIter);
    }

private:
    AtColor     evalSpecular(const AtVector &L) const
    {
        float LdotN = AiV3Dot(L, mAxisN);
        float VdotN = AiV3Dot(mViewDir, mAxisN);

        if (LdotN < AI_EPSILON || VdotN < AI_EPSILON) {
            return AI_RGB_BLACK;
        }

        // The half-vector is the microfacet surface normal
        AtVector M = AiV3Normalize(L + mViewDir);

        float LdotM = AiV3Dot(L, M);
        float NdotM = AiV3Dot(mAxisN, M);

        if (NdotM < AI_EPSILON || LdotM < AI_EPSILON) {
            return AI_RGB_BLACK;
        }

        float NdotM2 = SQR(AiV3Dot(mAxisN, M));

        float Ds = D_GTR2Aniso(M, NdotM2);
        float FH = powf(CLAMP(1.0f - LdotM, 0.0f, 1.0f), 5.0f);
        AtRGB Fs = LERP(FH, mSpecularF0, AI_RGB_WHITE);
        float Gs = smithG_GGX(LdotN, mSpecularRoughness) * smithG_GGX(VdotN, mSpecularRoughness);

        // The IOR for clearcoat is fixed to 1.5 (F0 = 0.04).
        static const float clearcoatF0 = 0.04f;
        static const float clearcoatRoughness = 0.25f;
        float Dr = D_GTR1(M, NdotM2);
        float Fr = LERP(FH, clearcoatF0, 1.0f);
        float Gr = smithG_GGX(LdotN, clearcoatRoughness) * smithG_GGX(VdotN, clearcoatRoughness);

        auto Fsheen = FH * mSheenColor * (1.0f - mMetallic);

        // The denominator for specular reflection (4 * LdotN * VdotN) has been already
        // separated into smithG_GGX. We don't have to divide it again.
        return (Ds * Fs * Gs + mClearcoat * Dr * Fr * Gr) + Fsheen;
    }

    static AtVector sphericalDirection(float cosTheta, float phi)
    {
        AtVector omega;
        omega.z = cosTheta;
        float r = sqrtf(1.0f - SQR(omega.z));
        omega.x = r * cosf(phi);
        omega.y = r * sinf(phi);
        return omega;
    }

    static AtVector getReflectDirection(const AtVector &i, const AtVector &m)
    {
        return 2.0f * ABS(AiV3Dot(i, m)) * m - i;
    }

    AtVector   concentricDiskSample(float rx, float ry) const
    {
        rx = rx * 2.0f - 1.0f;
        ry = ry * 2.0f - 1.0f;

        AtVector result;
        if (rx == 0.0f && ry == 0.0f) {
            // Handle degenerate case
            result.x = result.y = 0.0f;
            return result;
        }

        float r, phi;
        if (ABS(rx) > ABS(ry)) {
            r = rx;
            phi = AI_PIOVER2 * 0.5f * ry / rx;
        } else {
            r = ry;
            phi = AI_PIOVER2 * (1.0f - 0.5f * rx / ry);
        }

        result.x = r * cosf(phi);
        result.y = r * sinf(phi);
        return result;
    }

    //! Sample hemisphere according to cosine-weighted solidangle.
    AtVector    sampleDiffuseDirection(float rx, float ry) const
    {
        AtVector omega = concentricDiskSample(rx, ry);
        omega.z = sqrtf(MAX(0.0f, 1.0f - SQR(omega.x) - SQR(omega.y)));
        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return omega;
    }

    AtVector    sampleSpecularDirection(float rx, float ry) const
    {
        // Select lobe by the relative weights.
        // Sample the microfacet normal first, then compute the reflect direction.
        AtVector M;

        float gtr2Weight = 1.0f / (mClearcoat + 1.0f);

        if (rx < gtr2Weight) {
            rx /= gtr2Weight;
            M = mSampleFromVisibleNormal
                ? sampleGTR2AnisoDirectionFromSlope(rx, ry)
                : sampleGTR2AnisoDirection(rx, ry);
        } else {
            rx = (rx - gtr2Weight) / (1.0f - gtr2Weight);
            M = sampleGTR1Direction(rx, ry);
        }

        if (AiV3Dot(mAxisN, M) < 0.0f) {
            return AI_V3_ZERO;
        }

        return getReflectDirection(mViewDir, M);
    }

    //! Sample the direction according to GTR1
    AtVector    sampleGTR1Direction(float rx, float ry) const
    {
        auto phiH = AI_PITIMES2 * rx;
        auto a2 = SQR(mRoughness);
        auto cosThetaH = a2 == 1.0f
            ? sqrtf(1.0f - ry)
            : sqrtf((1.0f - pow(a2, 1.0f - ry)) / (1.0f - a2));

        AtVector omega = sphericalDirection(cosThetaH, phiH);
        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return AiV3Normalize(omega);
    }

    AtVector    sampleGTR2AnisoDirection(float rx, float ry) const
    {
        auto g = sqrtf(ry / (1.0f - ry));
        float phi = AI_PITIMES2 * rx;
        AtVector omega;
        AiV3Create(omega, g * mAlphaX * cosf(phi), g * mAlphaY * sinf(phi), 1.0f);
        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return AiV3Normalize(omega);
    }

    AtVector2    sampleGTR2AnisoSlope(float theta, float rx, float ry) const
    {
        AtVector2 slope;

        auto uniformSample = [] (float rx, float ry) {
            AtVector2 slope;
            float r = sqrtf(rx / (1.0f - rx));
            float phi = AI_PITIMES2 * ry;
            slope.x = r * cosf(phi);
            slope.y = r * sinf(phi);
            return slope;
        };

        if (theta < AI_EPSILON) {
            return uniformSample(rx, ry);
        }

        float B = tanf(theta);
        float B2 = SQR(B);
        float G1 = 2.0f / (1.0f + sqrtf(1.0f + B2));

        // sample slope_x
        float A = 2.0f * rx / G1 - 1.0f;
        float A2 = SQR(A);
        if (ABS(A2 - 1.0f) < AI_EPSILON) {
            return uniformSample(rx, ry);
        }

        float tmp = 1.0f / (A2 - 1.0f);
        float D = sqrtf(MAX(0.0f, B2 * SQR(tmp) - (A2 - B2) * tmp));
        float slopeX1 = B * tmp - D;
        float slopeX2 = B * tmp + D;
        slope.x = (A < 0.0f || slopeX2 > 1.0f / B) ? slopeX1 : slopeX2;

        // sample slope_y
        float sign = 1.0f;
        if (ry > 0.5f) {
            ry = 2.0f * (ry - 0.5f);
        } else {
            sign = -1.0f;
            ry = 2.0f * (0.5f - ry);
        }

        float z = (ry *(ry *(ry * 0.27385f - 0.73369f) + 0.46341f))
                / (ry *(ry *(ry * 0.093073f + 0.309420f) - 1.0f) + 0.597999f);
        slope.y = sign * z * sqrtf(1.0f + SQR(slope.x));
        return slope;
    }

    //! The sampling routine from Heitz and d'Eon.
    //! https://hal.inria.fr/hal-00996995v2/file/supplemental1.pdf
    AtVector    sampleGTR2AnisoDirectionFromSlope(float rx, float ry) const
    {
        auto V = mViewDir;

        // Transform to local frame
        float cosThetaV = CLAMP(AiV3Dot(mAxisN, V), -1.0f, 1.0f);
        assert(cosThetaV >= 0.0f);
        float phiV = atan2f(AiV3Dot(mAxisV, V), AiV3Dot(mAxisU, V));
        V = sphericalDirection(cosThetaV, phiV);

        // Stretch view direction
        V.x *= mAlphaX;
        V.y *= mAlphaY;
        V = AiV3Normalize(V);

        float theta = 0.0f;
        float phi = 0.0f;

        if (V.z < (1.0f - AI_EPSILON)) {
            theta = acosf(V.z);
            phi = atan2f(V.y, V.x);
        }

        auto slope = sampleGTR2AnisoSlope(theta, rx, ry);

        // Rotate & unstretch
        float cosPhi = cosf(phi);
        float sinPhi = sinf(phi);
        AtVector omega;
        omega.x = -(cosPhi * slope.x - sinPhi * slope.y) * mAlphaX;
        omega.y = -(sinPhi * slope.x + cosPhi * slope.y) * mAlphaY;
        omega.z = 1.0f;

        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return AiV3Normalize(omega);
    }

    AtVector    sampleGTR2Direction(float rx, float ry) const
    {
        auto thetaM = atanf(mRoughness * sqrtf(rx / (1.0f - rx)));
        auto phiM = AI_PITIMES2 * ry;

        AtVector omega = sphericalDirection(cosf(thetaM), phiM);
        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return omega;
    }

    //! Sample over cosine
    float   evalDiffusePdf(const AtVector &i) const
    {
        return MAX(1e-4f, AiV3Dot(i, mAxisN) * AI_ONEOVERPI);
    }

    float   evalSpecularPdf(const AtVector &i) const
    {
        auto m = AiV3Normalize(i + mViewDir);
        auto IdotM = ABS(AiV3Dot(i, m));
        auto MdotN = AiV3Dot(m, mAxisN);

        if (MdotN < 0.0f) {
            return 0.0f;
        }

        float MdotN2 = SQR(MdotN);

        auto clearcoatWeight = mClearcoat / (mClearcoat + 1.0f);

        if (mSampleFromVisibleNormal) {
            auto VdotN = MAX(1e-4f, AiV3Dot(mViewDir, mAxisN));
            auto Dw = smithG_GGX(IdotM, mSpecularRoughness) * D_GTR2Aniso(m, MdotN2) * 2.0f * IdotM / VdotN;
            auto D = LERP(clearcoatWeight, Dw, D_GTR1(m, MdotN2) * ABS(MdotN) / IdotM);
            return D * 0.25f;
        }

        auto D = LERP(clearcoatWeight, D_GTR2Aniso(m, MdotN2), D_GTR1(m, MdotN2));
        return D * ABS(MdotN) * 0.25f / IdotM;
    }

    float   D_GTR1(const AtVector &m, float MdotN2) const
    {
        auto alpha = LERP(mClearcoatGloss, 0.1f, 0.001f);
        auto a2 = SQR(alpha);
        float denominator = log(a2) * (1.0f + (a2 - 1.0f) * MdotN2);
        return (a2 - 1.0f) * AI_ONEOVERPI / denominator;
    }

    float   D_GTR2(const AtVector &m) const
    {
        float MdotN = AiV3Dot(m, mAxisN);
        float a2 = SQR(mRoughness);
        float denominator = AI_PI * SQR(1.0f + (a2 - 1.0f) * SQR(MdotN));
        return a2 / denominator;
    }

    float   D_GTR2Aniso(const AtVector &m, float MdotN2) const
    {
        float HdotU = AiV3Dot(m, mAxisU);
        float HdotV = AiV3Dot(m, mAxisV);

        float denominator = mAlphaX * mAlphaY * SQR(SQR(HdotU / mAlphaX) + SQR(HdotV / mAlphaY) + MdotN2);
        return AI_ONEOVERPI / denominator;
    }

    float smithG_GGX(float NdotV, float alphaG) const
    {
        // This is the original G term in Walter et al. [3] divided by 2 * NdotV.
        float a = alphaG * alphaG;
        float b = NdotV * NdotV;
        assert(a + b - a * b > 0.0f);
        return 1.0f / (NdotV + sqrt(a + b - a * b));
    }

private:
    AtVector    mViewDir;
    AtVector    mAxisU, mAxisV, mAxisN;
    AtColor     mSpecularF0;
    AtColor     mSheenColor;

    AtColor     mBaseColor;
    float       mRoughness;
    float       mSubsurface;
    float       mSpecular;
    float       mSpecularTint;
    float       mMetallic;
    float       mSheen;
    float       mSheenTint;
    float       mAnisotropic;
    float       mClearcoat;
    float       mClearcoatGloss;

    float       mSpecularRoughness;
    float       mAlphaX, mAlphaY;

    AtUInt16    mSampleType;
    bool        mSampleFromVisibleNormal;
};

node_parameters
{
    AiParameterRGB("base_color", 1.0f, 1.0f, 1.0f);

    std::vector<std::string> scalarAttrList = {
        "subsurface", "metallic", "specular", "specular_tint", "roughness",
        "anisotropic", "sheen", "sheen_tint", "clearcoat", "clearcoat_gloss" };

    std::set<std::string> softMaxAttrSet = {
        "specular", "roughness", "sheen" };

    for (auto attr : scalarAttrList) {
        AiParameterFlt(attr.c_str(), 0.0f);
        AiMetaDataSetFlt(mds, attr.c_str(), "min", 0.0f);
        auto maxName = softMaxAttrSet.find(attr) != softMaxAttrSet.end() ? "max" : "softmax";
        AiMetaDataSetFlt(mds, attr.c_str(), maxName, 1.0f);
    }

    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("indirectDiffuseScale", 1.0f);
    AiMetaDataSetFlt(mds, "indirectDiffuseScale", "min", 0.0f);
    AiMetaDataSetFlt(mds, "indirectDiffuseScale", "max", 1.0f);
    AiParameterFlt("indirectSpecularScale", 1.0f);
    AiMetaDataSetFlt(mds, "indirectSpecularScale", "min", 0.0f);
    AiMetaDataSetFlt(mds, "indirectSpecularScale", "max", 1.0f);

    AiParameterSTR("aov_direct_diffuse", "direct_diffuse");
    AiMetaDataSetInt(mds, "aov_direct_diffuse", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_direct_specular", "direct_specular");
    AiMetaDataSetInt(mds, "aov_direct_specular", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_indirect_diffuse", "indirect_diffuse");
    AiMetaDataSetInt(mds, "aov_indirect_diffuse", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_indirect_specular", "indirect_specular");
    AiMetaDataSetInt(mds, "aov_indirect_specular", "aov.type", AI_TYPE_RGB);
}

node_initialize
{
    // Dump samples
    //auto sg = AiShaderGlobals();
    //sg->Rd = -rls::sphericalDirection(cosf(1.5f), AI_PIOVER2 * 0.0f);
    //AiV3Create(sg->Nf, 0.0f, 0.0f, 1.0f);
    //DisneySampler brdf(node, sg);
    //rls::SampleWriter samplerWiter(512, 256);
    ////brdf.setSampleType(AI_RAY_DIFFUSE);
    //brdf.setSampleType(AI_RAY_GLOSSY);

    //samplerWiter.writeRadiance(sg, brdf);
    //samplerWiter.writeSample(sg, brdf, 50);

    ShaderData *data = new ShaderData;
    AiNodeSetLocalData(node, data);
}

node_update
{
    ShaderData *data = static_cast<ShaderData*>(AiNodeGetLocalData(node));
    if (data) {
        data->update();
    }
}

node_finish
{
    ShaderData *data = static_cast<ShaderData*>(AiNodeGetLocalData(node));

    if (data) {
        AiNodeSetLocalData(node, nullptr);
        delete data;
    }
}

shader_evaluate
{
    AtRGB opacity = AiShaderEvalParamFlt(p_opacity) * AI_RGB_WHITE;

    if (AiShaderGlobalsApplyOpacity(sg, opacity)) {
        return;
    }

    if (sg->Rt & AI_RAY_SHADOW) {
        sg->out_opacity = opacity;
        return;
    }

    DisneySampler sampler(node, sg);

    AtColor diffuse = AI_RGB_BLACK;
    AtColor specular = AI_RGB_BLACK;

    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg)) {
        if (AiLightGetAffectDiffuse(sg->Lp)) {
            diffuse += sampler.evalDiffuseLightSample(sg);
        }

        if (AiLightGetAffectSpecular(sg->Lp)) {
            specular += sampler.evalSpecularLightSample(sg);
        }
    }

    if (sg->Rt & (AI_RAY_DIFFUSE | AI_RAY_GLOSSY)) {
        diffuse *= AiShaderEvalParamFlt(p_indirect_diffuse);
        specular *= AiShaderEvalParamFlt(p_indirect_specular);
    }

    AtColor result = diffuse + specular;
    AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_direct_diffuse), diffuse);
    AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_direct_specular), specular);

    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));
    auto indirectDiffuse = data->shouldTraceDiffuse(sg) ? sampler.integrateDiffuse(sg, data) : AI_RGB_BLACK;
    auto indirectGlossy = data->shouldTraceGlossy(sg) ? sampler.integrateGlossy(sg) : AI_RGB_BLACK;
    //auto indirectGlossy = data->shouldTraceGlossy(sg) ? sampler.integrateGlossy(sg, data) : AI_RGB_BLACK;

    result += indirectDiffuse + indirectGlossy;
    AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_indirect_diffuse), indirectDiffuse);
    AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_indirect_specular), indirectGlossy);

    sg->out.RGB = result;
}