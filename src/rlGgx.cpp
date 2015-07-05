// Microfacet BSDF based on
// Microfacet Models for Refraction through Rough Surfaces. Bruce Walter et al. EGSR'07
// http://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.html
#include <algorithm>
#include <ai.h>

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
    GgxSampler(AtNode *node, AtShaderGlobals *sg)
    {
        mSpecularColor = AiShaderEvalParamRGB(p_Ks_color);
        
        bool isEntering = AiV3Dot(sg->N, sg->Rd) < AI_EPSILON;
        mIorIn = 1.0f;
        mIorOut = AiShaderEvalParamFlt(p_ior);
        if (!isEntering) {
            std::swap(mIorIn, mIorOut);
        }

        mViewDir = -sg->Rd;
        mAxisN = sg->Nf;        
        AiBuildLocalFramePolar(&mAxisU, &mAxisV, &mAxisN);

        // Remap the roughness value to approximate the highlight range in aiStandard.
        mRoughness = MAX(1e-5f, SQR(AiShaderEvalParamFlt(p_roughness)));
        //mRoughness = (1.2f - 0.2f * sqrtf(ABS(AiV3Dot(mViewDir, Nf)))) * mRoughness;
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
        auto thetaM = atanf((mRoughness) * sqrtf(rx / (1.0f - rx)));
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
            // Replace division with multiplication, since we only needs the sign.
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

    GgxSampler sampler(node, sg);
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
   
    auto roughness = AiShaderEvalParamFlt(p_roughness);
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