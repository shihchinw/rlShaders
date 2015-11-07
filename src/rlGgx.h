// Microfacet BSDF based on
// Microfacet Models for Refraction through Rough Surfaces. Bruce Walter et al. EGSR'07
// http://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html
#ifndef RLS_GGX_H
#define RLS_GGX_H

#include <algorithm>
#include <ai.h>

namespace rls
{

class   GgxSampler
{
public:
    //! Return sample direction
    static AtVector evalSample(const void *brdfData, float rx, float ry)
    {
        auto *brdf = static_cast<const GgxSampler*>(brdfData);
        auto M = brdf->sampleMicrofacetNormal(rx, ry);
        return brdf->getReflectDirection(brdf->mViewDir, M);
    }

    //! Return the reflectance
    static AtColor  evalBrdf(const void *brdfData, const AtVector *indir)
    {
        if (*indir == AI_V3_ZERO) {
            // TODO: Check how did this happen...
            return AI_RGB_BLACK;
        }

        auto *brdf = static_cast<const GgxSampler*>(brdfData);
        return brdf->evalReflectance(*indir);
    }

    static float evalPdf(const void *brdfData, const AtVector *indir)
    {
        auto *brdf = static_cast<const GgxSampler*>(brdfData);
        auto V = brdf->mViewDir;
        auto H = AiV3Normalize(V + *indir);
        return brdf->getReflectPdf(V, H);
    }

public:
    GgxSampler(AtShaderGlobals *sg, const AtColor &specColor, float ior, float roughness)
        : mSpecularColor(specColor)
    {
        bool isEntering = AiV3Dot(sg->N, sg->Rd) < AI_EPSILON;
        mIorIn = 1.0f;
        mIorOut = ior;
        if (!isEntering) {
            std::swap(mIorIn, mIorOut);
        }

        mViewDir = -sg->Rd;
        mAxisN = sg->Nf;
        AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);

        // Remap the roughness value to approximate the highlight range in aiStandard.
        mRoughness = MAX(1e-5f, SQR(roughness));
    }

    AtColor     evalReflectance(const AtVector &L) const
    {
        if (AiColorIsSmall(mSpecularColor)) {
            return AI_RGB_BLACK;
        }

        return mSpecularColor * reflection(mViewDir, L, mAxisN) * AiV3Dot(L, mAxisN);
    }

    AtColor     evalLightSample(AtShaderGlobals *sg)
    {
        return AiEvaluateLightSample(sg, this, evalSample, evalBrdf, evalPdf);
    }

    AtColor     integrateGlossy(AtShaderGlobals *sg)
    {
        if (AiColorIsSmall(mSpecularColor)) {
            return AI_RGB_BLACK;
        }

        return AiBRDFIntegrate(sg, this, evalSample, evalBrdf, evalPdf, AI_RAY_GLOSSY);
    }

    template <typename ShaderData>
    AtColor     traceShadow(AtShaderGlobals *sg, const ShaderData *data)
    {
        AtScrSample scrs;

        AtRay ray;
        AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Nf, AI_BIG, sg);
        AiRefractRay(&ray, &sg->Nf, mIorIn, mIorOut, sg);

        if (!data->shouldTraceRefract(sg) || !AiTrace(&ray, &scrs)) {
            // Sample environment? Or use exit color
            AiTraceBackground(&ray, &scrs);
            return scrs.color;
        }

        return scrs.color;
    }

    //! Sample indirect radiance.
    template <typename ShaderData>
    AtColor     integrateRefract(AtShaderGlobals *sg, const ShaderData *data)
    {
        AtRay       ray;
        AtScrSample scrs;
        AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, nullptr, AI_BIG, sg);

        if (!data->shouldTraceRefract(sg)) {
            bool isRefracted = AiRefractRay(&ray, &sg->Nf, mIorIn, mIorOut, sg);

            if (isRefracted) {
                auto weight = SQR(mIorOut / mIorIn) * ABS(AiV3Dot(sg->Nf, ray.dir));
                AiTraceBackground(&ray, &scrs);
                return scrs.color * weight;
            }

            return AI_RGB_BLACK;    // Total internal reflection!
        }

        AtSamplerIterator *sampleIter = data->getSamplerIter(sg);
        float   samples[2];
        AtColor result = AI_RGB_BLACK;

        while (AiSamplerGetSample(sampleIter, samples)) {
            AtVector m = sampleMicrofacetNormal(samples[0], samples[1]);
            bool isRefracted = AiRefractRay(&ray, &m, mIorIn, mIorOut, sg);

            if (!isRefracted) {
                // Total internal reflection
                AiReflectRay(&ray, &m, sg);
                //ray.dir = getReflectDirection(mViewDir, m);
            }

            if (!AiTrace(&ray, &scrs)) {
                AiTraceBackground(&ray, &scrs);
            }

            result += scrs.color * getSampleWeight(mViewDir, ray.dir, m);
        }

        return result * AiSamplerGetSampleInvCount(sampleIter);
    }

    //! Fresnel for dielectrics with unpolarized light. (Eq.22)
    float   fresnel(const AtVector &i, const AtVector &m) const
    {
        auto eta = mIorOut / mIorIn;
        auto f0 = SQR((eta - 1.0f) / (eta + 1.0f));
        auto m1 = -1.0f * SGN(AiV3Dot(i, m)) * m;
        return AiFresnelWeight(m1, i, f0);

        // TODO: Check the weird case when the IOR is less than one.
        float c = ABS(AiV3Dot(i, m));
        float gSqr = SQR(mIorOut / mIorIn) - 1.0f + c * c;

        if (gSqr < 0.0f) {
            // Total internal reflection
            return 1.0f;
        }

        float g = sqrtf(gSqr);
        float gmc = g - c;
        float gpc = g + c;

        return 0.5f * SQR(gmc / gpc) * (1.0f + SQR((c * gpc - 1.0f) / (c * gmc + 1.0f)));
    }

    float   G(const AtVector &i, const AtVector &o, const AtVector &m, const AtVector &n) const
    {
        return G1(i, m, n) * G1(o, m, n);
    }

private:
    //! Generate the normal of microfacet according to GGX distribution.
    AtVector    sampleMicrofacetNormal(float rx, float ry) const
    {
        auto thetaM = atanf(mRoughness * sqrtf(rx / (1.0f - rx)));
        auto phiM = AI_PITIMES2 * ry;

        AtVector omega;

        omega.z = cosf(thetaM);
        float r = sqrtf(1.0f - SQR(omega.z));
        omega.x = r * cosf(phiM);
        omega.y = r * sinf(phiM);

        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return omega;
    }

    AtVector    getReflectDirection(const AtVector &i, const AtVector &m) const
    {
        return 2.0f * ABS(AiV3Dot(i, m)) * m - i;
    }

    bool    getRefractDirection(const AtVector &m, const AtVector &i, AtVector &dir) const
    {
        int  sign = SGN(AiV3Dot(i, mAxisN));
        auto IdotM = (AiV3Dot(i, m));
        auto eta = mIorIn / mIorOut;

        float cosThetaTSqr = 1.0f + eta * (SQR(IdotM) - 1.0f);
        if (cosThetaTSqr < 0.0f) {
            // Total internal reflection
            return false;
        }

        dir = (eta * IdotM - sign * sqrtf(cosThetaTSqr)) * m - eta * i;
        return true;
    }

    //! The PDF of reflection. eq(38)
    float   getReflectPdf(const AtVector &i, const AtVector &m) const
    {
        auto IdotM = ABS(AiV3Dot(i, m));
        auto MdotN = ABS(AiV3Dot(m, mAxisN));
        return D(m, mAxisN) * MdotN * 0.25f / IdotM;
    }

    //! Return the sample weight for importance sampling of BSDF. (Eq.41)
    float   getSampleWeight(const AtVector &i, const AtVector &o, const AtVector &m) const
    {
        auto IdotH = (AiV3Dot(i, m));
        auto MdotN = ABS(AiV3Dot(m, mAxisN));
        auto IdotN = ABS(AiV3Dot(i, mAxisN));

        return G(i, o, m, mAxisN) * ABS(IdotH / (IdotN * MdotN));
    }

    //! The reflection term. (Eq.20)
    float   reflection(const AtVector &i, const AtVector &o, const AtVector &n) const
    {
        auto hr = static_cast<float>(SGN(AiV3Dot(i, n))) * AiV3Normalize(o + i);
        float reflectWeight = fresnel(i, hr);

        float LdotN = ABS(AiV3Dot(o, n));
        float VdotN = ABS(AiV3Dot(i, n));

        return reflectWeight * G(i, o, hr, n) * D(hr, n) * 0.25f / (LdotN * VdotN);
    }

    //! The refraction term. (Eq.21)
    float   refraction(const AtVector &i, const AtVector &o, const AtVector &n) const
    {
        auto ht = -AiV3Normalize(mIorIn * i + mIorOut * o);
        float refractWeight = 1.0f - fresnel(i, ht);

        float OdotN = ABS(AiV3Dot(o, n));
        float IdotN = ABS(AiV3Dot(i, n));
        float OdotH = (AiV3Dot(o, ht));
        float IdotH = (AiV3Dot(i, ht));

        float denominator = OdotN * IdotN * SQR(mIorIn * IdotH + mIorOut * OdotH);
        return ABS(OdotH * IdotH) * SQR(mIorOut) * refractWeight * G(i, o, ht, n) * D(ht, n) / denominator;
    }

    //! The normal distribution with GGX formula. (Eq.33)
    float   D(const AtVector &m, const AtVector &n) const
    {
        float MdotN = AiV3Dot(m, n);

        if (MdotN < 0.0f) {
            return 0.0f;
        }

        float roughnessSqr = SQR(mRoughness);
        float cosSqr = SQR(MdotN);
        float tanSqr = 1.0f / cosSqr - 1.0f;
        float denominator = AI_PI * SQR(cosSqr) * SQR(roughnessSqr + tanSqr);
        return roughnessSqr / denominator;
    }

    //! The shadowing/masking term. (Eq.34)
    float   G1(const AtVector &v, const AtVector &m, const AtVector &n) const
    {
        float VdotM = AiV3Dot(v, m);
        float VdotN = AiV3Dot(v, n);

        if (VdotM * VdotN < 0.0f) {
            // Replace division with multiplication, since we only need the sign.
            return 0.0f;
        }

        float cosSqr = SQR(VdotN);
        float tanSqr = 1.0f / cosSqr - 1.0f;
        float denominator = 1.0f + sqrtf(1.0f + SQR(mRoughness) * tanSqr);
        return 2.0f / denominator;
    }

private:
    AtColor     mSpecularColor;
    AtVector    mViewDir;
    AtVector    mAxisU, mAxisV, mAxisN;

    float       mRoughness;
    float       mIorIn, mIorOut;
};

} // Namespace of rls
#endif