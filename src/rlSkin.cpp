#include <ai.h>

#include "rlGgx.h"
#include "rlSss.h"

AI_SHADER_NODE_EXPORT_METHODS(SkinMethod);

namespace
{

enum   SkinParams
{
    p_sss_color,
    p_sss_weight,
    p_distance_multiplier,
    p_scatter_distance,
    p_cavity_fadeout,

    p_specular_color,
    p_specular_weight,
    p_specular_roughness,
    p_specular_ior,

    p_sheen_color,
    p_sheen_weight,
    p_sheen_roughness,
    p_sheen_ior,

    p_opacity,
    p_opacity_color,

    p_aov_sheen,
    p_aov_specular,
    p_aov_sss,
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
    AiParameterRGB("sss_color", 1.0f, 1.0f, 1.0f);
    AiMetaDataSetBool(mds, "sss_color", "always_linear", true);
    AiParameterFLT("sss_weight", 1.0f);
    AiParameterFLT("sss_dist_multiplier", 1.0f);
    AiParameterVec("sss_scatter_dist", 1.0f, 1.0f, 1.0f);

    AiParameterBool("sss_cavity_fadeout", true);
    AiMetaDataSetBool(mds, "sss_cavity_fadeout", "linkable", false);

    AiParameterRGB("specular_color", 1.0f, 1.0f, 1.0f);
    AiMetaDataSetBool(mds, "specular_color", "always_linear", true);
    AiParameterFLT("specular_weight", 0.6f);
    AiParameterFLT("specular_roughness", 0.5f);
    AiParameterFLT("specular_ior", 1.44f);

    AiParameterRGB("sheen_color", 1.0f, 1.0f, 1.0f);
    AiMetaDataSetBool(mds, "sheen_color", "always_linear", true);
    AiParameterFLT("sheen_weight", 0.0f);
    AiParameterFLT("sheen_roughness", 0.35f);
    AiParameterFLT("sheen_ior", 1.44f);

    AiParameterFLT("opacity", 1.0f);
    AiParameterRGB("opacity_color", 1.0f, 1.0f, 1.0f);

    AiParameterSTR("aov_sheen", "sheen");
    AiMetaDataSetInt(mds, "aov_sheen", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_specular", "specular");
    AiMetaDataSetInt(mds, "aov_specular", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_sss", "sss");
    AiMetaDataSetInt(mds, "aov_sss", "aov.type", AI_TYPE_RGB);
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
    AtRGB opacity = AiShaderEvalParamFlt(p_opacity) * AiShaderEvalParamRGB(p_opacity_color);

    if (sg->Rt & AI_RAY_SHADOW) {
        sg->out_opacity *= AiColorClamp(opacity, 0.0f, 1.0f);
        return;
    }

    float sheenFresnel = 0.0f;
    AtRGB sheen = AI_RGB_BLACK;

    float specularFresnel = 0.0f;
    AtRGB specular = AI_RGB_BLACK;

    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));
    AtNode *options = AiUniverseGetOptions();
    auto maxGlossyDepth = AiNodeGetInt(options, "GI_glossy_depth");

    // Evaluate sheen and specular
    if (sg->Rr_gloss < maxGlossyDepth) {
        auto sheenWeight = AiShaderEvalParamFlt(p_sheen_weight);
        auto sheenColor = AiShaderEvalParamRGB(p_sheen_color);
        auto sheenIor = AiShaderEvalParamFlt(p_sheen_ior);
        auto sheenRoughness = AiShaderEvalParamFlt(p_sheen_roughness);
        auto sheenF0 = SQR((sheenIor - 1.0f) / (sheenIor + 1.0f));
        sheenFresnel = AiFresnelWeight(sg->Nf, sg->Rd, sheenF0) * sheenWeight;

        if (sheenFresnel > AI_EPSILON) {
            rls::GgxSampler sheenSampler(sg, sheenColor, sheenIor, sheenRoughness);
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg)) {
                if (AiLightGetAffectSpecular(sg->Lp)) {
                    sheen += sheenSampler.evalLightSample(sg);
                }
            }

            if (sg->Rr == 0) {
                sheen += sheenSampler.integrateGlossy(sg);
            }
        }

        sheen *= sheenWeight;

        auto specularWeight = AiShaderEvalParamFlt(p_specular_weight);
        auto specularColor = AiShaderEvalParamRGB(p_specular_color);
        auto specularIor = AiShaderEvalParamFlt(p_specular_ior);
        auto specularRoughness = AiShaderEvalParamFlt(p_specular_roughness);

        auto specularF0 = SQR((specularIor - 1.0f) / (specularIor + 1.0f));
        specularFresnel = AiFresnelWeight(sg->Nf, sg->Rd, specularF0) * specularWeight;

        if (specularFresnel > AI_EPSILON) {
            rls::GgxSampler specularSampler(sg, specularColor, specularIor, specularRoughness);

            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg)) {
                if (AiLightGetAffectSpecular(sg->Lp)) {
                    specular += specularSampler.evalLightSample(sg);
                }
            }

            if (sg->Rr == 0) {
                specular += specularSampler.integrateGlossy(sg);
            }
        }

        specular *= specularWeight * (1.0f - sheenFresnel);
    }

    auto albedo = AiShaderEvalParamRGB(p_sss_color);
    auto distanceScale = AiShaderEvalParamFlt(p_distance_multiplier);
    auto scatterDist = AiShaderEvalParamVec(p_scatter_distance) * distanceScale;
    auto sssWeight = AiShaderEvalParamFlt(p_sss_weight);
    //sssWeight *= (1.0f - sheenFresnel) * (1.0f - specularFresnel);

    rls::SssSampler<rls::NDProfile> sssSampler(sg, albedo, scatterDist);
    auto sss = AI_RGB_BLACK;
    if (sssWeight > AI_EPSILON) {
        sss = sssSampler.integrateScatter(sg, data) * sssWeight;
    }

    if (sg->Rt & AI_RAY_CAMERA) {
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_sheen), sheen);
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_specular), specular);
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_sss), sss);
    }

    sg->out.RGB = sheen + specular + sss;
    sg->out_opacity = AiColorClamp(opacity, 0.0f, 1.0f);
}