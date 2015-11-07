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
#include "rlSss.h"

AI_SHADER_NODE_EXPORT_METHODS(SssNDMethod);

namespace
{

enum   SssParams
{
    p_base_color,
    p_distance_multiplier,
    p_scatter_distance,
    p_ior,
    p_cavity_fadeout,
    p_opacity,
};

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

    void    update(AtNode *node)
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

        mUseCavityFade = AiNodeGetBool(node, AtString("cavity_fadeout"));
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

    bool    useCavityFade() const
    {
        return mUseCavityFade;
    }

private:
    AtSampler   *mDiffuseSampler;
    AtSampler   *mSssSampler;
    int         mMaxRayDepth;
    int         mDiffuseDepth;
    int         mSssDepth;
    bool        mUseCavityFade;
};

}

//-----------------------------------------------------------------------------

node_parameters
{
    AiParameterRGB("color", 1.0f, 1.0f, 1.0f);
    AiMetaDataSetBool(mds, "color", "always_linear", true);

    AiParameterFLT("distance_multiplier", 1.0f);
    AiMetaDataSetFlt(mds, "distance_multiplier", "min", 0.0f);
    AiMetaDataSetFlt(mds, "distance_multiplier", "softmax", 5.0f);
    AiParameterVec("scatter_dist", 1.0f, 1.0f, 1.0f);

    AiParameterFLT("ior", 1.0f);
    AiMetaDataSetFlt(mds, "ior", "min", 0.0f);

    AiParameterBool("cavity_fadeout", true);
    AiMetaDataSetBool(mds, "cavity_fadeout", "linkable", false);

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
        data->update(node);
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

    auto albedo = AiShaderEvalParamRGB(p_base_color);
    auto distanceScale = AiShaderEvalParamFlt(p_distance_multiplier);
    auto scatterDist = AiShaderEvalParamVec(p_scatter_distance) * distanceScale;

    rls::SssSampler<rls::NDProfile> sampler(sg, albedo, scatterDist);
    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));

    // TODO: Use GGX for glossy component?
    sg->out.RGB = sampler.integrateScatter(sg, data);
    sg->out_opacity = AiColorClamp(opacity, 0.0f, 1.0f);
}