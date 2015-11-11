#include <algorithm>
#include <ai.h>

#include "rlGgx.h"

AI_SHADER_NODE_EXPORT_METHODS(GgxMethod);

//-----------------------------------------------------------------------------

AtVector2   rls::VNDFKernel::sampleSlope(float theta, float rx, float ry) const
{
    AtVector2 slope;

    auto uniformSample = [](float rx, float ry) {
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

AtVector rls::VNDFKernel::evalSample(float rx, float ry) const
{
    auto V = mViewDir;

    // Transform to local frame
    float cosThetaV = CLAMP(AiV3Dot(mBasis.N, V), -1.0f, 1.0f);
    //assert(cosThetaV >= 0.0f);

    float phiV = atan2f(AiV3Dot(mBasis.V, V), AiV3Dot(mBasis.U, V));
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

    auto slope = sampleSlope(theta, rx, ry);

    // Rotate & unstretch
    float cosPhi = cosf(phi);
    float sinPhi = sinf(phi);
    AtVector omega;
    omega.x = -(cosPhi * slope.x - sinPhi * slope.y) * mAlphaX;
    omega.y = -(sinPhi * slope.x + cosPhi * slope.y) * mAlphaY;
    omega.z = 1.0f;

    AiV3RotateToFrame(omega, mBasis.U, mBasis.V, mBasis.N);
    return AiV3Normalize(omega);
}

//-----------------------------------------------------------------------------

namespace
{

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
    p_anisotropic,
    p_opacity,
    p_opacity_color,

    p_aov_direct_diffuse,
    p_aov_direct_specular,
    p_aov_refract,
    p_aov_indirect_diffuse,
    p_aov_indirect_specular,
};

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

    AiParameterRGB("KsColor", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("Ks", 0.5f);

    AiParameterRGB("KtColor", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("Kt", 0.0f);
    AiParameterFLT("ior", 1.0f);

    AiParameterFLT("roughness", 0.0f);
    AiParameterFLT("anisotropic", 0.0f);

    AiParameterFLT("opacity", 1.0f);
    AiParameterRGB("opacity_color", 1.0f, 1.0f, 1.0f);

    AiParameterSTR("aov_direct_diffuse", "direct_diffuse");
    AiMetaDataSetInt(mds, "aov_direct_diffuse", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_direct_specular", "direct_specular");
    AiMetaDataSetInt(mds, "aov_direct_specular", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_refract", "refraction");
    AiMetaDataSetInt(mds, "aov_refract", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_indirect_diffuse", "indirect_diffuse");
    AiMetaDataSetInt(mds, "aov_indirect_diffuse", "aov.type", AI_TYPE_RGB);
    AiParameterSTR("aov_indirect_specular", "indirect_specular");
    AiMetaDataSetInt(mds, "aov_indirect_specular", "aov.type", AI_TYPE_RGB);
}

node_initialize
{
    // Dump samples
    //auto sg = AiShaderGlobals();
    //sg->Rd = -rls::sphericalDirection(cosf(AI_PIOVER2 * 0.95f), AI_PIOVER2 * 0.0f);
    ////sg->Rd = -rls::sphericalDirection(cosf(AI_PIOVER2 * 0.5f), AI_PIOVER2 * 0.0f);
    //AiV3Create(sg->Nf, 0.0f, 0.0f, 1.0f);

    //auto specColor = AiShaderEvalParamRGB(p_Ks_color);
    //auto ior = AiShaderEvalParamFlt(p_ior);
    //auto roughness = AiShaderEvalParamFlt(p_roughness);
    //auto anisotropic = AiShaderEvalParamFlt(p_anisotropic);

    //rls::GgxSampler brdf(sg, specColor, ior, roughness, anisotropic);
    //rls::SampleWriter samplerWiter(512, 256);

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
    AtRGB opacity = AiShaderEvalParamFlt(p_opacity) * AiShaderEvalParamRGB(p_opacity_color);

    if (AiShaderGlobalsApplyOpacity(sg, opacity)) {
        return;
    }

    auto specColor = AiShaderEvalParamRGB(p_Ks_color);
    auto ior = AiShaderEvalParamFlt(p_ior);
    auto roughness = AiShaderEvalParamFlt(p_roughness);
    auto anisotropic = AiShaderEvalParamFlt(p_anisotropic);

    rls::GgxSampler sampler(sg, specColor, ior, roughness, anisotropic);
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
    bool sampleDiffuse = !AiColorIsSmall(diffuseColor) && sg->Rr_diff <= maxDiffuseDepth;

    AtColor diffuse = AI_RGB_BLACK;
    AtColor specular = AI_RGB_BLACK;

    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg)) {
        if (AiLightGetAffectDiffuse(sg->Lp) && sampleDiffuse) {
            diffuse += AiEvaluateLightSample(sg, diffData,
                            AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
        }

        if (AiLightGetAffectSpecular(sg->Lp) && sg->Rr_gloss <= maxGlossyDepth) {
            specular += sampler.evalLightSample(sg);
        }
    }

    auto specularWeight = AiShaderEvalParamFlt(p_Ks);
    diffuse *= diffuseColor;
    specular *= specularWeight;

    auto transmissionWeight = AiShaderEvalParamFlt(p_Kt);
    auto ktColor = AiShaderEvalParamRGB(p_Kt_Color) * transmissionWeight;
    AtColor transmission = AiColorIsSmall(ktColor) ? AI_RGB_BLACK : sampler.integrateRefract(sg, data) * ktColor;

    AtColor result = diffuse + specular + transmission;

    if (sg->Rt & AI_RAY_CAMERA) {
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_direct_diffuse), diffuse);
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_direct_specular), specular);
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_refract), transmission);

        auto indirectDiffuse = sampleDiffuse ? AiOrenNayarIntegrate(&sg->Nf, sg, roughness) : AI_RGB_BLACK;
        auto indirectGlossy = sampler.integrateGlossy(sg) * specularWeight;

        result += diffuseColor * indirectDiffuse + indirectGlossy;
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_indirect_diffuse), indirectDiffuse);
        AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_indirect_specular), indirectGlossy);
    }

    sg->out.RGB = result;
    sg->out_opacity = AiColorClamp(opacity, 0.0f, 1.0f);
}