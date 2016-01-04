// GGX BRDF/BSDF based on
// [1]  Microfacet Models for Refraction through Rough Surfaces. Bruce Walter et al. EGSR'07
//      http://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html
// [2]  Physically Based Shading at Disney. Brent Burley. SIGGRAPH 2012 course notes.
//      http://blog.selfshadow.com/publications/s2012-shading-course/
// [3] 	Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals.
//      Eric Heitz and Eugene d'Eon. EGSR'14
//      https://hal.inria.fr/hal-00996995/en
//      https://hal.inria.fr/hal-00996995v2/file/supplemental1.pdf
#ifndef RLS_GGX_H
#define RLS_GGX_H

#include <algorithm>
#include <memory>
#include <cassert>

#include <ai.h>

#include "rlUtil.h"

namespace rls
{

class   NDFKernel
{
public:
    NDFKernel(const AtVector &v, float ax, float ay, const CoordBasis &frame)
        : mBasis(frame), mViewDir(v), mAlphaX(ax), mAlphaY(ay)
    {
    }

    //! Sample the normal vector of microfacet. [2] Eq.14 p.26
    AtVector    evalSample(float rx, float ry) const
    {
        auto g = sqrtf(rx / (1.0f - rx));
        float phi = AI_PITIMES2 * ry;
        AtVector omega;
        AiV3Create(omega, g * mAlphaX * cosf(phi), g * mAlphaY * sinf(phi), 1.0f);
        AiV3RotateToFrame(omega, mBasis.U, mBasis.V, mBasis.N);
        return AiV3Normalize(omega);
    }

    //! The PDF for GGX reflection.  [1] Eq.38
    template <typename Sampler>
    float   evalPdf(const Sampler *ggxSampler, const AtVector &i, const AtVector &m)
    {
        auto IdotM = ABS(AiV3Dot(i, m));
        auto MdotN = ABS(AiV3Dot(m, mBasis.N));
        return ggxSampler->D(m, mBasis.N) * MdotN * 0.25f / IdotM;
    }

protected:
    const CoordBasis    &mBasis;
    AtVector            mViewDir;
    float               mAlphaX, mAlphaY;
};

//! The sampling mechanism with the visible normal distribution fuction. [3]
class   VNDFKernel
{
public:
    VNDFKernel(const AtVector &v, float ax, float ay, const CoordBasis &frame)
        : mBasis(frame), mViewDir(v), mAlphaX(ax), mAlphaY(ay)
    {
    }

    //! Sample the normal vector of microfacet.
    AtVector    evalSample(float rx, float ry) const;

    //! Return PDF with given incoming direction and microfacet normal.
    template <typename Sampler>
    float   evalPdf(const Sampler *ggxSampler, const AtVector &i, const AtVector &m)
    {
        const auto &N = mBasis.N;
        auto IdotM = ABS(AiV3Dot(i, m));
        auto MdotN = ABS(AiV3Dot(m, N));
        auto IdotN = ABS(AiV3Dot(i, N));
        auto pdf = ggxSampler->D(m, N) * ggxSampler->G1(i, m, N) / IdotN * 0.25f;
        return MAX(pdf, AI_EPSILON);
    }

private:
    AtVector2   sampleSlope(float theta, float rx, float ry) const;

private:
    const CoordBasis    &mBasis;
    AtVector            mViewDir;
    float               mAlphaX, mAlphaY;
};


template <typename NDSampler>
class   GgxSamplerT
{
public:
    //! Return sample direction
    static AtVector evalSample(const void *brdfData, float rx, float ry)
    {
        auto brdf = static_cast<const GgxSamplerT *>(brdfData);
        auto M = brdf->mNormalSampler->evalSample(rx, ry);
        auto L = rls::reflectDirection(brdf->mViewDir, M);

        brdf->mReflectWeight += brdf->fresnel(L, M);
        brdf->mMisSampleCount++;

        return L;
    }

    //! Return the reflectance
    static AtColor  evalBrdf(const void *brdfData, const AtVector *indir)
    {
        if (*indir == AI_V3_ZERO) {
            // TODO: Check how did this happen...
            return AI_RGB_BLACK;
        }

        auto brdf = static_cast<const GgxSamplerT *>(brdfData);
        return brdf->evalReflectance(*indir);
    }

    static float evalPdf(const void *brdfData, const AtVector *indir)
    {
        auto brdf = static_cast<const GgxSamplerT *>(brdfData);
        auto V = brdf->mViewDir;
        auto H = AiV3Normalize(V + *indir);
        return brdf->mNormalSampler->evalPdf(brdf, V, H);
    }

public:
    GgxSamplerT(AtShaderGlobals *sg, const AtColor &specColor,
               float ior, float roughness, float anisotropic = 0.0f)
        : mSpecularColor(specColor)
        , mNormalSampler(nullptr)
        , mReflectWeight(0.0f)
        , mMisSampleCount(0.0f)
    {
        bool isEntering = AiV3Dot(sg->N, sg->Rd) < AI_EPSILON;
        mIorIn = 1.0f;
        mIorOut = ior;
        if (!isEntering) {
            std::swap(mIorIn, mIorOut);
        }

        mViewDir = -sg->Rd;
        mBasis.N = mAxisN = sg->Nf;
        AiBuildLocalFramePolar(&mBasis.U, &mBasis.V, &mBasis.N);

        float aspect = sqrt(1.0f - anisotropic * 0.9f);
        mAlphaX = MAX(1e-4f, SQR(roughness) / aspect);
        mAlphaY = MAX(1e-4f, SQR(roughness) * aspect);

        mNormalSampler = std::make_shared<NDSampler>(mViewDir, mAlphaX, mAlphaY, mBasis);

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

    float       getAvgReflectWeight() const
    {
        return mMisSampleCount > 0.0f ? mReflectWeight / mMisSampleCount : 1.0f;
    }

    //template <typename ShaderData>
    //AtColor     traceShadow(AtShaderGlobals *sg, const ShaderData *data)
    //{
    //    AtScrSample scrs;

    //    AtRay ray;
    //    AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Nf, AI_BIG, sg);
    //    AiRefractRay(&ray, &sg->Nf, mIorIn, mIorOut, sg);

    //    if (!data->shouldTraceRefract(sg) || !AiTrace(&ray, &scrs)) {
    //        // Sample environment? Or use exit color
    //        AiTraceBackground(&ray, &scrs);
    //        return scrs.color;
    //    }

    //    return scrs.color;
    //}

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
            AtVector m = mNormalSampler->evalSample(samples[0], samples[1]);
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

    //! Fresnel for dielectrics with unpolarized light. [1] Eq.22
    float   fresnel(const AtVector &i, const AtVector &m) const
    {
        auto eta = mIorOut / mIorIn;
        auto f0 = SQR((eta - 1.0f) / (eta + 1.0f));
        auto m1 = -1.0f * SGN(AiV3Dot(i, m)) * m;
        //return AiFresnelWeight(m1, i, f0);    // This would return negative value sometimes?!

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

    //! Return the sample weight for importance sampling of BSDF. [1] Eq.41
    float   getSampleWeight(const AtVector &i, const AtVector &o, const AtVector &m) const
    {
        auto IdotH = (AiV3Dot(i, m));
        auto MdotN = ABS(AiV3Dot(m, mAxisN));
        auto IdotN = ABS(AiV3Dot(i, mAxisN));

        return G(i, o, m, mAxisN) * ABS(IdotH / (IdotN * MdotN));
    }

    //! The reflection term. [1] Eq.20
    float   reflection(const AtVector &i, const AtVector &o, const AtVector &n) const
    {
        auto hr = static_cast<float>(SGN(AiV3Dot(i, n))) * AiV3Normalize(o + i);
        float reflectWeight = fresnel(i, hr);

        float LdotN = ABS(AiV3Dot(o, n));
        float VdotN = ABS(AiV3Dot(i, n));

        return reflectWeight * G(i, o, hr, n) * D(hr, n) * 0.25f / (LdotN * VdotN);
    }

    //! The refraction term. [1] Eq.21
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

    //! The anisotropic normal distribution with GGX formula. [2] Eq.13, p.25
    //! @note The isotropic version can be referred to [1] Eq.33
    float   D(const AtVector &m, const AtVector &n) const
    {
        float MdotU = AiV3Dot(m, mBasis.U);
        float MdotV = AiV3Dot(m, mBasis.V);
        float MdotN2 = SQR(AiV3Dot(mAxisN, m));

        float denominator = mAlphaX * mAlphaY * SQR(SQR(MdotU / mAlphaX) + SQR(MdotV / mAlphaY) + MdotN2);
        return AI_ONEOVERPI / denominator;
    }

    //! The shadowing/masking term. [1] Eq.34
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
    std::shared_ptr<NDSampler>    mNormalSampler;

    CoordBasis      mBasis;
    AtColor         mSpecularColor;
    AtVector        mViewDir;
    AtVector        mAxisN;

    float           mRoughness;
    float           mAlphaX, mAlphaY;
    float           mIorIn, mIorOut;

    mutable float   mReflectWeight;
    mutable float   mMisSampleCount;
};

using GgxSampler = GgxSamplerT<VNDFKernel>;

} // Namespace of rls
#endif