// BSSRDF importance sampling for normalized diffusion.
// [1]  Approximate Reflectance Profiles for Efficient Subsurface Scattering.
//      Per H. Christensen, Brent Burley.
//      http://graphics.pixar.com/library/ApproxBSSRDF/index.html
// [2]  Extending the Disney BRDF to a BSDF with Integrated Subsurface Scattering.
//      Brent Burley. SIG'15 course notes.
//      http://blog.selfshadow.com/publications/s2015-shading-course/burley/s2015_pbs_disney_bsdf_notes.pdf
// [3]  BSSRDF Importance Sampling. Alan King et al. SIG'13 Talks.
//      https://www.solidangle.com/research/s2013_bssrdf_slides.pdf
#include <algorithm>
#include <vector>
#include <array>
#include <assert.h>

#include <ai.h>

#include "rlUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(SssNDMethod);

namespace
{

enum   SssParams
{
    p_base_color,
    p_distance_multiplier,
    p_scatter_distance,
    p_ior,
    p_opacity,
};

union FloatInt
{
    float   fValue;
    unsigned int     iValue;
};

const char *gRlsRayType = "rls_raytype";
const char *gRlsSssMsgData = "rls_sss_msg_data";
const char *gRlsSssOrigObj = "rls_sss_orig_obj";

struct  ShaderData
{
    ShaderData()
        : mDiffuseSampler(nullptr)
        , mSssSampler(nullptr)
    {
    }

    virtual ~ShaderData()
    {
        AiSamplerDestroy(mDiffuseSampler);
        AiSamplerDestroy(mSssSampler);
    }

    void    update()
    {
        AtNode *options = AiUniverseGetOptions();

        mMaxRayDepth = AiNodeGetInt(options, "GI_total_depth");
        mDiffuseDepth = AiNodeGetInt(options, "GI_diffuse_depth");

        auto diffuseSampleNum = AiNodeGetInt(options, "GI_diffuse_samples");
        AiSamplerDestroy(mDiffuseSampler);
        mDiffuseSampler = AiSampler(diffuseSampleNum, 2);

        auto sampleNum = AiNodeGetInt(options, "sss_bssrdf_samples");
        AiSamplerDestroy(mSssSampler);
        mSssSampler = AiSampler(sampleNum, 2);
    }

    AtSamplerIterator *getDiffuseSamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mDiffuseSampler, sg);
    }

    AtSamplerIterator *getSssSamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mSssSampler, sg);
    }

    bool    shouldTraceDiffuse(AtShaderGlobals *sg) const
    {
        return (sg->Rr_diff < mDiffuseDepth && sg->Rr < mMaxRayDepth);
    }

    bool    shouldTraceSss(AtShaderGlobals *sg) const
    {
        return sg->Rr < mMaxRayDepth;
    }

private:
    AtSampler   *mDiffuseSampler;
    AtSampler   *mSssSampler;
    int         mMaxRayDepth;
    int         mDiffuseDepth;
    int         mSssDepth;
};

}

//-----------------------------------------------------------------------------

//! Normalized diffuse profile
class   NDProfile
{
public:
    static int     selectDistLobe(float &x)
    {
        if (x < 0.3333f) {
            x = LINEARSTEP(0.0f, 0.3333f, x);
            return 0;
        } else if (x > 0.6666f) {
            x = LINEARSTEP(0.6666f, 1.0f, x);
            return 2;
        }

        x = LINEARSTEP(0.3333f, 0.6666f, x);
        return 1;
    }

    float maxRadius() const
    {
        return mMaxRadius;
    }

    void    setDistance(const AtVector &dist)
    {
        /*auto ad = dist.x + dist.y + dist.z;
        auto len = AiV3Length(dist);

        for (auto i = 0; i < 3; i++) {
            mDistance[i] = (1.0f - dist[i] / ad) * len;
        }*/

        mDistance = dist;
        // Scale the max tracing distance empirically.
        mMaxRadius = MAX(mDistance.x, MAX(mDistance.y, mDistance.z)) * 3.0f;

        for (auto i = 0; i < 3; i++) {
            auto d = mDistance[i];
            mC1[i] = 1.0f - exp(-mMaxRadius / d);
            mC2[i] = 1.0f - exp(-mMaxRadius / d / 3.0f);
        }
    }

    float   getRadius(float rx) const
    {
        if (mMaxRadius < AI_EPSILON) {
            return 0.0f;
        }

        // Since the profile is integrated to one for any positive d over the disk,
        // thus we could uniformly select the distance value.
        int distIdx = selectDistLobe(rx);
        auto d = mDistance[distIdx];
        if (d < AI_EPSILON) {
            return 0.0f;
        }

        /*auto w1 = 1.0f - exp(-mMaxDist / d);
        auto w2 = 1.0f - exp(-mMaxDist / d / 3.0f);*/
        auto w1 = mC1[distIdx];
        auto w2 = mC2[distIdx];
        auto w = w1 / (w1 + w2 * 3.0f);

        auto r = 1.0f;
        if (rx > w) {
            rx = LINEARSTEP(w, 1.0f, rx);
            r = log(1.0f - rx * w2) * (-d * 3.0f);
        } else {
            rx = LINEARSTEP(0.0f, w, rx);
            r = log(1.0f - rx * w1) * (-d);
        }

        return r;
    }

    float   getPdf(float r) const
    {
        if (mMaxRadius < AI_EPSILON) {
            return 1.0f;
        }

        float pdf = 0.0f;

        for (auto i = 0u; i < 3; i++) {
            auto d = MAX(mDistance[i], AI_EPSILON);
            auto p1 = exp(-r / d);
            auto p2 = exp(-r / d / 3.0f);
            pdf += (p1 + p2) / d / (mC1[i] + mC2[i] * 3.0f);
        }

        return pdf / (AI_PITIMES2 * r * 3.0f);
    }

    AtRGB   evalProfile(float r) const
    {
        if (mMaxRadius < AI_EPSILON) {
            // TODO: Replace with simple BRDF?
            return AI_RGB_BLACK;
        } else if (r < AI_EPSILON) {
            return AI_RGB_WHITE;
        }

        auto denom = 8.0f * AI_PI * r;

        AtRGB result = AI_RGB_BLACK;
        for (auto i = 0u; i < 3; i++) {
            // Remap the distance for artistic control
            auto d = mDistance[i];
            result[i] = d < AI_EPSILON ? 1.0f
                : (exp(-r / d) + exp(-r / (3.0f * d))) / (denom * d);
        }

        return result;
    }

private:
    AtVector    mDistance;
    AtVector    mC1, mC2;
    float       mMaxRadius;
};


class   GaussianProfile
{
public:
    float maxRadius() const
    {
        return mMaxRadius;
    }

    void    setDistance(const AtVector &dist)
    {
        mMaxRadius = dist.x;
        mVariance = SQR(mMaxRadius) / 12.46f;
        mNorm = 1.0f - fast_exp(-SQR(mMaxRadius) * 0.5f / mVariance);
    }

    float   getRadius(float rx) const
    {
        return sqrt(-2.0f * mVariance * log(1.0f - rx * mNorm));
    }

    float   getPdf(float r) const
    {
        return evalProfile(r) / mNorm;
    }

    float   evalProfile(float r) const
    {
        return AI_ONEOVER2PI / mVariance * fast_exp(-r * r * 0.5f / mVariance);
    }

private:
    float   mVariance;
    float   mMaxRadius;
    float   mNorm;
};



//! Sub-surface scattering based on normalized diffusion.
class   SssSampler
{
private:
    static const int kMaxProbeDepth = 16;

    enum SssRayType
    {
        kSssRayUndefined = 0x00,
        kSssProbeRay,
    };

    struct ProbeSample
    {
        AtRGB       irradiance;
        AtVector    N;
        AtVector    disp;
        float       r;
    };

    struct MsgData
    {
        AtVector        No;
        AtVector        Po;
        double          maxDist;
        unsigned int    probeDepth;
        std::array<ProbeSample, kMaxProbeDepth> probeSampleArray;
    };

public:
    SssSampler(AtNode *node, AtShaderGlobals *sg)
    {
        mBaseColor = AiShaderEvalParamRGB(p_base_color);
        auto distanceScale = AiShaderEvalParamFlt(p_distance_multiplier);
        auto scatterDist = AiShaderEvalParamVec(p_scatter_distance);
        mProfile.setDistance(scatterDist * distanceScale);

        // We use geometry normal as the up-axis during probe tracing.
        mAxisN = sg->Ns;

        //AiBuildLocalFrameShirley(&mAxisU, &mAxisV, &mAxisN);
        AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);
        AiM4Frame(mWorldToLocalMat, &sg->P, &mAxisU, &mAxisV, &mAxisN);

        auto eta = AiShaderEvalParamFlt(p_ior);
        mEnableFresnel = eta > AI_EPSILON;

        mKrt = SQR((eta - 1.0f) / (eta + 1.0f));
        mFtOut = 1.0f - AiFresnelWeight(mAxisN, sg->Rd, mKrt);
    }

    AtColor     integrateScatter(AtShaderGlobals *sg, const ShaderData *data)
    {
        if (sg->Rt & AI_RAY_DIFFUSE) {
            // Use OrenNayar as the approximation of diffuse component.
            auto directDiffuse = AI_RGB_BLACK;
            auto diffuseData = AiOrenNayarMISCreateData(sg, 0.0f);
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg)) {
                if (AiLightGetAffectDiffuse(sg->Lp)) {
                    auto irradiance = AiEvaluateLightSample(sg, diffuseData,
                        AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
                    directDiffuse += irradiance * AiLightGetDiffuse(sg->Lp);
                }
            }

            return mBaseColor * directDiffuse;
        }

        int rayType = kSssRayUndefined;
        AiStateGetMsgInt(gRlsRayType, &rayType);

        MsgData *msgData = nullptr;
        AiStateGetMsgPtr(gRlsSssMsgData, reinterpret_cast<void**>(&msgData));

        if (kSssProbeRay == rayType) {
            shadeProbeSample(sg, data, msgData);
            // This is a dummy return value, all the shading results for each
            // probe sample are stored in the msgData.
            return AI_RGB_BLACK;
        } else if (!msgData) {
            // Allocate message data structure for further tracing.
            msgData = (MsgData*)AiShaderGlobalsQuickAlloc(sg, sizeof(MsgData));
            memset(msgData, 0, sizeof(MsgData));
            AiStateSetMsgPtr(gRlsSssMsgData, msgData);
        }

        // Scatter samples around the shading point.

        // The traced probe wound not succeessfully return the intersection
        // within the same triangle of current shading point. Change sg->fi
        // to workaround this problem.
        AtUInt32 oldPrimitiveId = sg->fi;
        sg->fi = UINT_MAX;

        AtRay ray;
        AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, nullptr, AI_BIG, sg);
        auto sampleIdx = 0;
        auto result = AI_RGB_BLACK;

        AtScrSample scrs;
        AtSamplerIterator *sampleIter = data->getSssSamplerIter(sg);
        float   samples[2];

        while (AiSamplerGetSample(sampleIter, samples)) {
            float r = getProbeRay(samples[0], samples[1], sampleIdx++, sg->P, ray);

            AiStateSetMsgInt(gRlsRayType, kSssProbeRay);
            AiStateSetMsgPtr(gRlsSssOrigObj, sg->Op);

            msgData->probeDepth = 0;
            msgData->maxDist = ray.maxdist;
            msgData->Po = sg->P;
            msgData->No = sg->Ns;

            // AiTrace would split the shading tree with ray type "kSssProbeRay".
            if (AiTrace(&ray, &scrs)) {
                // Combine the shading results from each proble samples.
                for (auto i = 0u; i < msgData->probeDepth; i++) {
                    auto &sample = msgData->probeSampleArray[i];

                    if (AiColorIsZero(sample.irradiance)) continue;

                    // Apply MIS techniques to compute the combined pdf and geometry terms.
                    AtVector offset;
                    AiM4VectorByMatrixMult(&offset, mWorldToLocalMat, &sample.disp);
                    offset *= offset;   // component-wise multiplication

                    float rr[3];
                    rr[0] = sqrt(offset[1] + offset[2]);
                    rr[1] = sqrt(offset[0] + offset[2]);
                    rr[2] = sqrt(offset[0] + offset[1]);

                    float pdf = mProfile.getPdf(rr[0]) * ABS(AiV3Dot(mAxisU, sample.N)) * 0.25f
                              + mProfile.getPdf(rr[1]) * ABS(AiV3Dot(mAxisV, sample.N)) * 0.25f
                              + mProfile.getPdf(rr[2]) * ABS(AiV3Dot(mAxisN, sample.N)) * 0.5f;

                    result += sample.irradiance / pdf;
                }
            }
        }

        sg->fi = oldPrimitiveId;
        AiStateSetMsgInt(gRlsRayType, kSssRayUndefined);

        return result * mBaseColor * AiSamplerGetSampleInvCount(sampleIter);
    }

private:
    void    alignDir(const AtVector refDir, AtVector &dir)
    {
        if (AiV3Dot(refDir, dir) < 0.0f) {
            dir *= -1.0f;
        }
    }

    //! Perform shading along probe ray.
    void    shadeProbeSample(AtShaderGlobals *sg, const ShaderData *data, MsgData *msgData)
    {
        if (sg->psg && AiNodeGetNodeEntry(sg->shader) != AiNodeGetNodeEntry(sg->psg->shader)) {
            return;
        }

        const auto probeDepth = msgData->probeDepth;
        if (probeDepth >= kMaxProbeDepth) return;

        void *origObjPtr;
        AiStateGetMsgPtr(gRlsSssOrigObj, &origObjPtr);

        // When probe ray hits other objects, just skip this sample.
        if (origObjPtr != sg->Op) {
            stepProbeRay(sg, msgData);
            return;
        }

        auto &sample = msgData->probeSampleArray[probeDepth];

        sample.disp = sg->P - msgData->Po;   // the offset from current pos of probe ray to the original shading point
        sample.r = AiV3Length(sample.disp);
        msgData->maxDist -= sg->Rl;
        //sample.irradiance = AI_RGB_BLACK;

        if (sample.r > mProfile.maxRadius()) {
            return;
        }

        // Set state to undefined for evaluation of light samples.
        AiStateSetMsgInt(gRlsRayType, kSssRayUndefined);

        // Align normal directions for shading calculation
        AtVector Nref = sg->N;
        alignDir(Nref, sg->Nf);
        alignDir(Nref, sg->Ngf);

        sample.N = sg->Ns;

        auto directResult = evalLightSample(sg, sample.r);
        auto indirectResult = integrateDiffuse(sg, data, sample.r);
        sample.irradiance = directResult + indirectResult;
        msgData->probeDepth++;

        AiStateSetMsgInt(gRlsRayType, kSssProbeRay);
        stepProbeRay(sg, msgData);
    }

    void    stepProbeRay(AtShaderGlobals *sg, MsgData *msgData)
    {
        if (msgData->probeDepth < kMaxProbeDepth && msgData->maxDist > 0.0f) {
            AtRay       ray;
            AtScrSample scrs;
            AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, msgData->maxDist, sg);
            AiTrace(&ray, &scrs);
        }
    }

    //! Evaluate direct lighting.
    AtRGB   evalLightSample(AtShaderGlobals *sg, float r)
    {
        AtRGB result = AI_RGB_BLACK;
        //sg->fhemi = false;

        auto diffuseData = AiOrenNayarMISCreateData(sg, 0.0f);
        AiLightsPrepare(sg);
        while (AiLightsGetSample(sg)) {
            if (AiLightGetAffectDiffuse(sg->Lp)) {
                auto irradiance = AiEvaluateLightSample(sg, diffuseData,
                     AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
                result += irradiance * AiLightGetDiffuse(sg->Lp);
            }
        }

        //AiLightsPrepare(sg);
        //while (AiLightsGetSample(sg)) {
        //    auto d = AiLightGetDiffuse(sg->Lp);
        //    auto NdotL = CLAMP(ABS(AiV3Dot(sg->Ld, sg->Nf)), 0.0f, 1.0f);
        //    AtRGB L = AI_ONEOVERPI * sg->Li * sg->we * NdotL * d;
        //    if (!AiColorIsZero(L)) {
        //        // TODO: Support multiple components
        //        result += mProfile.evalProfile(r) * L;
        //    }
        //}
        //AiLightsResetCache(sg);

        return result * mProfile.evalProfile(r);
    }

    AtRGB   integrateDiffuse(AtShaderGlobals *sg, const ShaderData *data, float r)
    {
        if (!data->shouldTraceDiffuse(sg)) {
            return AI_RGB_BLACK;
        }

        AtRGB result = AI_RGB_BLACK;
        AtRay ray;
        AiMakeRay(&ray, AI_RAY_DIFFUSE, &sg->P, &sg->N, AI_BIG, sg);

        AtScrSample scrs;

        AtSamplerIterator *sampleIter = data->getDiffuseSamplerIter(sg);
        float samples[2];
        while (AiSamplerGetSample(sampleIter, samples)) {
            ray.dir = sampleDiffuseDirection(samples[0], samples[1]);

            if (AiTrace(&ray, &scrs)) {
                // TODO: Support different diffusion profile interfaces such as directional diffusion.
                auto NdotL = CLAMP(AiV3Dot(sg->Nf, ray.dir), 0.0f, 1.0f);
                result += scrs.color * NdotL;
            }
        }

        assert(AiSamplerGetSampleInvCount(sampleIter) == 1);
        return result * AI_ONEOVERPI * mProfile.evalProfile(r);
    }

    //! Return the probe ray to find intersection near x_0.
    //! Based on BSSRDF importance sampling by Arnold.
    float   getProbeRay(float rx, float ry, int idx, const AtVector &origin, AtRay &ray) const
    {
        float r = mProfile.getRadius(rx);
        assert(AiIsFinite(r) && r >= 0.0f);

        float rmax = mProfile.maxRadius();

        float theta = AI_PITIMES2 * ry;
        AtVector offset;
        offset.x = cosf(theta) * r;
        offset.z = sinf(theta) * r;
        offset.y = sqrt(rmax * rmax - r * r);

        // Choose the probe axis from the last two bits.
        AtVector disp;

        FloatInt fi;
        fi.fValue = rx;
        idx = fi.iValue;

        if ((idx & 0x03) < 2) {
            // Use normal as probe direction
            ray.dir = -mAxisN;
            disp = mAxisU * offset.x + mAxisV * offset.z - ray.dir * offset.y;
        } else if ((idx & 0x03) == 2) {
            ray.dir = (idx & 0x04) > 0 ? -mAxisU : mAxisU;
            disp = mAxisV * offset.x + mAxisN * offset.z - ray.dir * offset.y;
        } else {
            ray.dir = (idx & 0x04) > 0 ? -mAxisV : mAxisV;
            disp = mAxisN * offset.x + mAxisU * offset.z - ray.dir * offset.y;
        }

        ray.origin = origin + disp;
        ray.maxdist = offset.y * 2.0f;

        return r;
    }

    //! Sample hemisphere according to cosine-weighted solidangle.
    AtVector    sampleDiffuseDirection(float rx, float ry) const
    {
        AtVector omega = rls::concentricDiskSample(rx, ry);
        omega.z = sqrtf(MAX(0.0f, 1.0f - SQR(omega.x) - SQR(omega.y)));
        AiV3RotateToFrame(omega, mAxisU, mAxisV, mAxisN);
        return omega;
    }

private:
    NDProfile       mProfile;
    //GaussianProfile mProfile;
    AtMatrix    mWorldToLocalMat;
    AtColor     mBaseColor;
    AtVector    mAxisU, mAxisV, mAxisN;
    float       mFtOut, mKrt;
    bool        mEnableFresnel;
};

//-----------------------------------------------------------------------------

node_parameters
{
    AiParameterRGB("color", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("distance_multiplier", 1.0f);
    AiMetaDataSetFlt(mds, "distance_multiplier", "min", 0.0f);
    AiMetaDataSetFlt(mds, "distance_multiplier", "softmax", 5.0f);
    AiParameterVec("scatter_dist", 1.0f, 1.0f, 1.0f);

    AiParameterFLT("ior", 1.0f);
    AiMetaDataSetFlt(mds, "ior", "min", 0.0f);

    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
}

node_initialize
{
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

    if (sg->Rt & AI_RAY_SHADOW) {
        sg->out_opacity *= AiColorClamp(opacity, 0.0f, 1.0f);
        return;
    }

    SssSampler sampler(node, sg);
    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));
    // TODO: Use GGX for glossy component?
    sg->out.RGB = sampler.integrateScatter(sg, data);
    sg->out_opacity = AiColorClamp(opacity, 0.0f, 1.0f);
}