#ifndef RLS_SSS_H
#define RLS_SSS_H

#include <assert.h>
#include <array>

#include <ai.h>

#include "rlUtil.h"

#define STACKLESS_PROBE_TRACING

namespace rls
{

extern const char *gRlsRayType;
extern const char *gRlsSssMsgData;
extern const char *gRlsSssOrigObj;

//! The "Normalized Diffusion" profile from Disney.
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

    inline float   maxRadius() const
    {
        return mMaxRadius;
    }

    void    setDistance(const AtVector &dist, const AtColor &albedo);

    float   getRadius(float rx) const;

    float   getPdf(float r) const;

    AtRGB   evalProfile(float r) const;

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

    void    setDistance(const AtVector &dist, const AtColor &albedo)
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
template <typename Profile>
class   SssSampler
{
private:
    static const int kMaxProbeDepth = 12;

    union FloatParts
    {
        float value;

        struct
        {
            unsigned int sign : 1;
            unsigned int exponent : 8;
            unsigned int mantissa : 23;
        } parts;
    };

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
    SssSampler(AtShaderGlobals *sg, const AtRGB &albedo, const AtVector &dist)
        : mBaseColor(albedo)
    {
        mProfile.setDistance(dist, mBaseColor);

        // We use geometry normal as the up-axis during probe tracing.
        mAxisN = sg->Ns;

        if (!AiV3isZero(sg->dPdu) && AiV3Exists(sg->dPdu)) {
            AtVector U = AiV3Normalize(sg->dPdu);
            mAxisV = AiV3Normalize(AiV3Cross(mAxisN, U));
            mAxisU = AiV3Cross(mAxisV, mAxisN);
        } else {
            AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);
            //AiBuildLocalFrameShirley(&mAxisU, &mAxisV, &mAxisN);
        }

        AiM4Frame(mWorldToLocalMat, &sg->P, &mAxisU, &mAxisV, &mAxisN);

        /*auto eta = AiShaderEvalParamFlt(p_ior);
        mEnableFresnel = eta > AI_EPSILON;

        mKrt = SQR((eta - 1.0f) / (eta + 1.0f));
        mFtOut = 1.0f - AiFresnelWeight(mAxisN, sg->Rd, mKrt);*/
    }

    template <typename ShaderData>
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
        auto result = AI_RGB_BLACK;

        AtSamplerIterator *sampleIter = data->getSssSamplerIter(sg);
        float   samples[2];

        while (AiSamplerGetSample(sampleIter, samples)) {
            float r = getProbeRay(samples[0], samples[1], sg->P, ray);

            AiStateSetMsgInt(gRlsRayType, kSssProbeRay);
            AiStateSetMsgPtr(gRlsSssOrigObj, sg->Op);

            msgData->probeDepth = 0;
            msgData->maxDist = ray.maxdist;
            msgData->Po = sg->P;
            msgData->No = sg->Ns;

#ifdef STACKLESS_PROBE_TRACING
            traceProbe(sg, ray, data, msgData);
#else
            // AiTrace would split the shading tree with ray type "kSssProbeRay".
            AtScrSample scrs;
            AiTrace(&ray, &scrs);
#endif
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

        sg->fi = oldPrimitiveId;
        AiStateSetMsgInt(gRlsRayType, kSssRayUndefined);

        return mBaseColor * result * AiSamplerGetSampleInvCount(sampleIter);
    }

private:
    void    alignDir(const AtVector &refDir, AtVector &dir)
    {
        if (AiV3Dot(refDir, dir) < 0.0f) {
            dir *= -1.0f;
        }
    }

    template <typename ShaderData>
    bool    traceProbe(AtShaderGlobals *sg, AtRay &ray, const ShaderData *data, MsgData *msgData)
    {
        auto sgOld = *sg;

        int trialCount = kMaxProbeDepth;
        AtShaderGlobals sgOut;
        sg->fhemi = false;

        while (trialCount > 0 && AiTraceProbe(&ray, &sgOut)) {
            --trialCount;

            if (sgOut.Op != sg->Op) {
                continue;
            }

            auto dist = AiV3Dist(sg->P, sgOut.P);
            if (dist > AI_EPSILON) {
                sg->P = sgOut.P;
                sg->Po = sgOut.Po;

                sg->N = sgOut.N;
                sg->Nf = sgOut.Nf;
                sg->Ng = sgOut.Ng;
                sg->Ngf = sgOut.Ngf;
                sg->Ns = sgOut.Ns;
                sg->Rl = sgOut.Rl;
                sg->fi = sgOut.fi;

                sg->Rr = sgOut.Rr;
                sg->Rt = sgOut.Rt;
                sg->area = sgOut.area;
                sg->Ro = sgOut.Ro;
                sg->Rd = sgOut.Rd;
                sg->bu = sgOut.bu;
                sg->bv = sgOut.bv;
                sg->dPdu = sgOut.dPdu;
                sg->dPdv = sgOut.dPdv;
                sg->dPdx = sgOut.dPdx;
                sg->dPdy = sgOut.dPdy;
                sg->dNdx = sgOut.dNdx;
                sg->dNdy = sgOut.dNdy;
                shadeProbeSample(sg, data, msgData);
            }

            ray.origin = sgOut.P;
            ray.maxdist = msgData->maxDist;
            sg->fi = UINT_MAX;
            sgOut.fi = UINT_MAX;

            if (msgData->maxDist <= 0.0f) {
                break;
            }
        }

        *sg = sgOld;
        return msgData->probeDepth > 0;
    }

    //! Perform shading along probe ray.
    template <typename ShaderData>
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
        sample.irradiance = AI_RGB_BLACK;

        if (sample.r > mProfile.maxRadius()) {
            return;
        }

        // Set state to undefined for evaluation of light samples.
        AiStateSetMsgInt(gRlsRayType, kSssRayUndefined);

        // Align normal directions for shading calculation
        AtVector Nref = sg->N;
        alignDir(Nref, sg->Nf);
        alignDir(Nref, sg->Ngf);
        alignDir(Nref, sg->Ns);
        alignDir(Nref, sg->Ng);
        sample.N = sg->Ns;

        auto cavityFade = 1.0f;
        if (data->useCavityFade()) {
            auto dispDir = sample.disp / sample.r;
            if (AiV3Dot(msgData->No, dispDir) < 0.0f) {
                // Use cos for gradation.
                auto cosCavityAngle = ABS(AiV3Dot(sample.N, msgData->No));
                cavityFade = sqrt((1.0f + cosCavityAngle) * 0.5f);
            } else {
                // Use cos(theta/2) to approximate the reduce rate of diffusion due to cavity.
                auto cosCavityAngle = CLAMP(AiV3Dot(sample.N, msgData->No), -1.0f, 1.0f);
                cavityFade = sqrt((1.0f + cosCavityAngle) * 0.5f);
            }
        }

        if (cavityFade > AI_EPSILON) {
            auto directResult = evalLightSample(sg, sample.r);
            auto indirectResult = integrateDiffuse(sg, data, sample.r);
            sample.irradiance = (directResult + indirectResult) * cavityFade;
            msgData->probeDepth++;
        }

        AiStateSetMsgInt(gRlsRayType, kSssProbeRay);
        stepProbeRay(sg, msgData);
    }

    void    stepProbeRay(AtShaderGlobals *sg, MsgData *msgData)
    {
#ifndef STACKLESS_PROBE_TRACING
        if (msgData->probeDepth < kMaxProbeDepth && msgData->maxDist > 0.0f) {
            AtRay       ray;
            AtScrSample scrs;
            AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, msgData->maxDist, sg);
            AiTrace(&ray, &scrs);
        }
#endif
    }

    //! Evaluate direct lighting.
    AtRGB   evalLightSample(AtShaderGlobals *sg, float r)
    {
        AtRGB result = AI_RGB_BLACK;

        auto diffuseData = AiOrenNayarMISCreateData(sg, 0.0f);
        AiLightsPrepare(sg);
        while (AiLightsGetSample(sg)) {
            if (AiLightGetAffectDiffuse(sg->Lp)) {
                auto irradiance = AiEvaluateLightSample(sg, diffuseData,
                    AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
                result += irradiance * AiLightGetDiffuse(sg->Lp);
            }
        }

        return result * mProfile.evalProfile(r);
    }

    template <typename ShaderData>
    AtRGB   integrateDiffuse(AtShaderGlobals *sg, const ShaderData *data, float r)
    {
        if (!data->shouldTraceDiffuse(sg)) {
            return AI_RGB_BLACK;
        }

        AtRGB result = AI_RGB_BLACK;
        AtRay ray;
        AiMakeRay(&ray, AI_RAY_DIFFUSE, &sg->P, nullptr, AI_BIG, sg);

        AtScrSample scrs;
        AtSamplerIterator *sampleIter = data->getDiffuseSamplerIter(sg);
        float samples[2];

        while (AiSamplerGetSample(sampleIter, samples)) {
            ray.dir = sampleDiffuseDirection(samples[0], samples[1], sg->Nf);

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
    float   getProbeRay(float rx, float ry, const AtVector &origin, AtRay &ray) const
    {
        auto idx = 0;

        if (rx < 0.5f) {
            idx = 0;
            rx = LINEARSTEP(0.0f, 0.5f, rx);
        } else if (rx < 0.75f) {
            idx = 2;
            rx = LINEARSTEP(0.5f, 0.75f, rx);
        } else {
            idx = 3;
            rx = LINEARSTEP(0.75f, 1.0f, rx);
        }

        float r = mProfile.getRadius(rx);
        assert(AiIsFinite(r) && r >= 0.0f);

        float rmax = mProfile.maxRadius();
        float phi = AI_PITIMES2 * ry;

        AtVector offset;
        offset.x = cosf(phi) * r;
        offset.z = sinf(phi) * r;
        offset.y = sqrt(rmax * rmax - r * r);
        ray.maxdist = offset.y * 2.0f;

        // Choose the probe axis from the last two bits of mantissa
        /*FloatParts fi;
        fi.value = rx;
        idx = fi.parts.mantissa;*/

        if ((idx & 0x03) < 2) {
            // Use normal as probe direction
            ray.dir = -mAxisN;
            AiV3RotateToFrame(offset, mAxisU, -ray.dir, mAxisV);
        } else if ((idx & 0x03) == 2) {
            ray.dir = (idx & 0x04) > 0 ? -mAxisU : mAxisU;
            AiV3RotateToFrame(offset, mAxisV, -ray.dir, mAxisN);
        } else {
            ray.dir = (idx & 0x04) > 0 ? -mAxisV : mAxisV;
            AiV3RotateToFrame(offset, mAxisN, -ray.dir, mAxisU);
        }

        ray.origin = origin + offset;
        return r;
    }

    //! Sample hemisphere according to cosine-weighted solidangle.
    AtVector    sampleDiffuseDirection(float rx, float ry, const AtVector &normal) const
    {
        AtVector omega = rls::concentricDiskSample(rx, ry);
        omega.z = sqrtf(MAX(0.0f, 1.0f - SQR(omega.x) - SQR(omega.y)));

        AtVector u, v;
        AiBuildLocalFramePolar(&u, &v, &normal);
        AiV3RotateToFrame(omega, u, v, normal);
        return omega;
    }

private:
    Profile     mProfile;
    AtMatrix    mWorldToLocalMat;
    AtColor     mBaseColor;
    AtVector    mAxisU, mAxisV, mAxisN;
    float       mFtOut, mKrt;
    bool        mEnableFresnel;
};

} // Namespace rls
#endif