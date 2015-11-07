// Microfacet BSDF based on
// Microfacet Models for Refraction through Rough Surfaces. Bruce Walter et al. EGSR'07
// http://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html
#include <algorithm>
#include <ai.h>

#include "rlGgx.h"

AI_SHADER_NODE_EXPORT_METHODS(GgxMethod);

enum    GgxParams
{
    p_Kd_color,
    p_Kd,
    p_Ks_color,
    p_Ks,
    p_Kt_Color,
    p_Kt,
    p_ior,
    p_roughness,
    p_opacity,
};

namespace
{

struct  ShaderData
{
    ShaderData() : mRefractionSampler(nullptr)
    {
    }

    virtual ~ShaderData()
    {
        AiSamplerDestroy(mRefractionSampler);
    }

    void    update()
    {
        AtNode *options = AiUniverseGetOptions();

        mMaxRayDepth = AiNodeGetInt(options, "GI_total_depth");
        mRefractionDepth = AiNodeGetInt(options, "GI_refraction_depth");

        auto refractionSampleNum = AiNodeGetInt(options, "GI_refraction_samples");
        AiSamplerDestroy(mRefractionSampler);
        mRefractionSampler = AiSampler(refractionSampleNum, 2);
    }

    bool    shouldTraceRefract(AtShaderGlobals *sg) const
    {
        return (sg->Rr_refr < mRefractionDepth && sg->Rr < mMaxRayDepth);
    }

    AtSamplerIterator *getSamplerIter(AtShaderGlobals *sg) const
    {
        return AiSamplerIterator(mRefractionSampler, sg);
    }

private:
    AtSampler   *mRefractionSampler;
    int         mMaxRayDepth;
    int         mRefractionDepth;
};

}


node_parameters
{
    AiParameterRGB("KdColor", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("Kd", 0.5f);
    AiMetaDataSetFlt(mds, "Kd", "min", 0.0f);
    AiMetaDataSetFlt(mds, "Kd", "max", 1.0f);

    AiParameterRGB("KsColor", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("Ks", 0.5f);
    AiMetaDataSetFlt(mds, "Ks", "min", 0.0f);
    AiMetaDataSetFlt(mds, "Ks", "max", 1.0f);

    AiParameterRGB("KtColor", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("Kt", 0.0f);
    AiMetaDataSetFlt(mds, "Kt", "min", 0.0f);
    AiMetaDataSetFlt(mds, "Kt", "max", 1.0f);


    AiParameterFLT("ior", 1.0f);
    AiMetaDataSetFlt(mds, "ior", "min", 0.0f);

    AiParameterFLT("roughness", 0.0f);
    AiMetaDataSetFlt(mds, "roughness", "min", 0.0f);
    AiMetaDataSetFlt(mds, "roughness", "max", 1.0f);

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

    if (AiShaderGlobalsApplyOpacity(sg, opacity)) {
        return;
    }

    auto specColor = AiShaderEvalParamRGB(p_Ks_color);
    auto ior = AiShaderEvalParamFlt(p_ior);
    auto roughness = AiShaderEvalParamFlt(p_roughness);

    rls::GgxSampler sampler(sg, specColor, ior, roughness);
    ShaderData *data = static_cast<ShaderData *>(AiNodeGetLocalData(node));

    if (sg->Rt & AI_RAY_SHADOW) {
        auto shadowColor = AiShaderEvalParamRGB(p_Kt_Color) * AiShaderEvalParamFlt(p_Kt);
        //auto shadowColor = sampler.traceShadow(sg, data);
        sg->out_opacity = AI_RGB_WHITE - shadowColor;
        return;
    }

    AtNode *options = AiUniverseGetOptions();
    auto maxDiffuseDepth = AiNodeGetInt(options, "GI_diffuse_depth");
    auto maxGlossyDepth = AiNodeGetInt(options, "GI_glossy_depth");

    auto diffData = AiOrenNayarMISCreateData(sg, roughness);

    auto diffuseWeight = AiShaderEvalParamFlt(p_Kd);
    auto diffuseColor = AiShaderEvalParamRGB(p_Kd_color) * diffuseWeight;
    bool sampleDiffuse = !AiColorIsSmall(diffuseColor) && sg->Rr_diff < maxDiffuseDepth;

    AtColor diffuse = AI_RGB_BLACK;
    AtColor specular = AI_RGB_BLACK;

    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg)) {
        if (AiLightGetAffectDiffuse(sg->Lp) && sampleDiffuse) {
            diffuse += AiEvaluateLightSample(sg, diffData,
                            AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
        }

        if (AiLightGetAffectSpecular(sg->Lp) && sg->Rr_gloss < maxGlossyDepth) {
            //specular += sampler.evalBrdf(&sampler, &sg->Ld) * sg->Li * sg->we;
            specular += sampler.evalLightSample(sg);
        }
    }

    auto specularWeight = AiShaderEvalParamFlt(p_Ks);

    auto transmissionWeight = AiShaderEvalParamFlt(p_Kt);
    auto ktColor = AiShaderEvalParamRGB(p_Kt_Color) * transmissionWeight;
    AtColor transmission = AiColorIsSmall(ktColor) ? AI_RGB_BLACK : sampler.integrateRefract(sg, data) * ktColor;

    AtColor result = diffuseColor * diffuse + specularWeight * specular + transmission;

    if (sg->Rr > 0) {
        sg->out.RGB = result;
        return;
    }

    auto indirectDiffuse = sampleDiffuse ? AiOrenNayarIntegrate(&sg->Nf, sg, roughness) : AI_RGB_BLACK;
    auto indirectGlossy = sampler.integrateGlossy(sg) * specularWeight;

    result += diffuseColor * indirectDiffuse + indirectGlossy;
    sg->out.RGB = result;
}