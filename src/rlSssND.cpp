#include <algorithm>
#include <vector>
#include <assert.h>
#include <ai.h>

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

struct  ShaderData
{
    ShaderData() : mSssSampler(nullptr)
    {
    }

    virtual ~ShaderData()
    {
        AiSamplerDestroy(mSssSampler);
    }

    void    update()
    {
        AtNode *options = AiUniverseGetOptions();

        mMaxRayDepth = AiNodeGetInt(options, "GI_total_depth");

        auto sampleNum = AiNodeGetInt(options, "sss_bssrdf_samples");
        AiSamplerDestroy(mSssSampler);
        mSssSampler = AiSampler(sampleNum, 2);
    }

    AtSamplerIterator *getSamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mSssSampler, sg);
    }

    bool    shouldTraceSss(AtShaderGlobals *sg) const
    {
        return sg->Rr < mMaxRayDepth;
    }

private:
    AtSampler   *mSssSampler;
    int         mMaxRayDepth;
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

    float maxDistance() const
    {
        return mMaxDist;
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
        mMaxDist = MAX(mDistance.x, MAX(mDistance.y, mDistance.z)) * 20.0f;

        for (auto i = 0; i < 3; i++) {
            auto d = mDistance[i];
            mC1[i] = 1.0f - exp(-mMaxDist / d);
            mC2[i] = 1.0f - exp(-mMaxDist / d / 3.0f);
        }
    }

    float   getRadius(float rx) const
    {
        if (mMaxDist < AI_EPSILON) {
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
        if (mMaxDist < AI_EPSILON) {
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
        if (mMaxDist < AI_EPSILON) {
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
    float       mMaxDist;
};


class   GaussianProfile
{
public:
    float maxDistance() const
    {
        return mMaxRadius;
    }

    void    setDistance(const AtVector &dist)
    {
        mVariance = dist.x;
        mMaxRadius = sqrt(12.46f * mVariance);
    }

    float   getRadius(float rx) const
    {
        auto x = 1.0f - exp(-SQR(mMaxRadius) / 2.0f / mVariance);
        return sqrt(-2.0f * mVariance * log(1.0f - rx * x));
    }

    float   getPdf(float r) const
    {
        float denom = 1.0f - exp(-SQR(mMaxRadius) / 2.0f / mVariance);
        return evalProfile(r) / denom;
    }

    float   evalProfile(float r) const
    {
        return AI_ONEOVER2PI / mVariance * exp(-r * r / 2.0f / mVariance);
    }

private:
    float   mVariance;
    float   mMaxRadius;
};

//! Sub-surface scattering based on normalized diffusion.
class   SssNDSampler
{
public:
    SssNDSampler(AtNode *node, AtShaderGlobals *sg)
    {
        mBaseColor = AiShaderEvalParamRGB(p_base_color);
        auto distanceScale = AiShaderEvalParamFlt(p_distance_multiplier);
        auto scatterDist = AiShaderEvalParamVec(p_scatter_distance);
        mProfile.setDistance(scatterDist * distanceScale);

        // We use geometry normal as the up-axis during probe tracing.
        mAxisN = sg->Ngf;
        AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);

        auto eta = AiShaderEvalParamFlt(p_ior);
        mEnableFresnel = eta > AI_EPSILON;

        mKrt = SQR((eta - 1.0f) / (eta + 1.0f));
        mFtOut = 1.0f - AiFresnelWeight(mAxisN, sg->Rd, mKrt);
    }

    AtColor     integrateScatter(AtShaderGlobals *sg, const ShaderData *data)
    {
        if (!data->shouldTraceSss(sg)) {
            return AI_RGB_BLACK;
        }

        // The traced probe wound not succeessfully return the intersection
        // within the same triangle of current shading point. We have to
        // switch shading context other than AI_CONTEXT_SURFACE.
        auto oldContext = sg->sc;
        sg->sc = AI_CONTEXT_VOLUME;
        AtRay ray;
        AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, nullptr, mProfile.maxDistance(), sg);

        AtSamplerIterator *sampleIter = data->getSamplerIter(sg);
        float   samples[2];
        AtColor result = AI_RGB_BLACK;
        AtShaderGlobals sgProbe;

        auto origin = sg->P;

        std::vector<AtShaderGlobals> inScatterPoints;
        std::vector<float> radiusArray;

        auto sampleIdx = 0;
        while (AiSamplerGetSample(sampleIter, samples)) {
            float r = getProbeDir(samples[0], samples[1], sampleIdx++, origin, ray);

            if (AiTraceProbe(&ray, &sgProbe)) {
                inScatterPoints.push_back(sgProbe);
                radiusArray.push_back(r);
            }
        }

        sg->sc = oldContext;

        for (auto idx = 0; idx < inScatterPoints.size(); idx++) {
            auto r = radiusArray[idx];
            auto &sgIn = inScatterPoints[idx];
            sg->P = sgIn.P;
            sg->Nf = sgIn.Nf;
            sg->N = sgIn.N;
            sg->Ngf = sgIn.Ngf;
            sg->Ng = sgIn.Ng;
            sg->Ns = sgIn.Ns;

            // Sample lights at each in-scattering point
            AtRGB inScatter = AI_RGB_BLACK;
            auto N = sgIn.Nf;

            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg)) {
                auto L = sg->Ld;
                auto NdotL = CLAMP(AiV3Dot(N, L), 0.0f, 1.0f);
                auto radiance = sg->Li * sg->we * NdotL;

                if (mEnableFresnel) {
                    float FtIn = 1.0f - AiFresnelWeight(N, -L, mKrt);
                    radiance *= FtIn * mFtOut;
                }

                inScatter += radiance;
            }
            AiLightsResetCache(sg);

            // Compute the refract weight for in/out scattering.
            auto combineWeight = ABS(AiV3Dot(mAxisN, N)) * 0.5f
                + ABS(AiV3Dot(mAxisU, N)) * 0.25f
                + ABS(AiV3Dot(mAxisV, N)) * 0.25f;

            auto dist = AiV3Dist(origin, sgIn.P);
            auto profileWeight = mProfile.evalProfile(dist) / mProfile.getPdf(r) / combineWeight;
            result += inScatter * profileWeight;
        }

        auto validSampleCount = static_cast<float>(inScatterPoints.size());

        if (validSampleCount == 0.0f)
            validSampleCount = 1.0f;

        return result * AI_ONEOVERPI * mBaseColor / validSampleCount;
        //return result * AI_ONEOVERPI * mBaseColor * AiSamplerGetSampleInvCount(sampleIter);
    }

private:
    //! Return the probe ray to find intersection near x_0.
    //! Based on BSSRDF importance sampling by Arnold.
    float   getProbeDir(float rx, float ry, int idx, const AtVector &origin, AtRay &ray) const
    {
        float r = mProfile.getRadius(rx);
        assert(r >= 0.0f);

        if (!AiIsFinite(r)) {
            r = 0.0f;
        }

        float theta = AI_PITIMES2 * ry;
        AtVector offset;
        offset.x = cosf(theta) * r;
        offset.z = sinf(theta) * r;
        offset.y = 1.0f;

        // Use idx choose the axis from the last two bits.
        AtVector dir;

        if ((idx & 0x03) < 2) {
            // Use normal as probe direction
            ray.dir = -mAxisN;
            dir = mAxisU * offset.x + mAxisV * offset.z - ray.dir * offset.y;
        } else if ((idx & 0x03) == 2) {
            ray.dir = (idx & 0x04) > 0 ? -mAxisU : mAxisU;
            dir = mAxisV * offset.x + mAxisN * offset.z - ray.dir* offset.y;
        } else {
            ray.dir = (idx & 0x04) > 0 ? -mAxisV : mAxisV;
            dir = mAxisU * offset.x + mAxisN * offset.z - ray.dir* offset.y;
        }

        ray.origin = origin + dir;

        return r;
    }

private:
    NDProfile   mProfile;
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

    SssNDSampler sampler(node, sg);
    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));

    // TODO: Use ggx for reflection?
    sg->out.RGB = sampler.integrateScatter(sg, data);
    sg->out_opacity = AiColorClamp(opacity, 0.0f, 1.0f);
}