#include "rlSss.h"
#include "rlUtil.h"
using namespace rls;

const char *rls::gRlsRayType = "rls_raytype";
const char *rls::gRlsSssMsgData = "rls_sss_msg_data";
const char *rls::gRlsSssOrigObj = "rls_sss_orig_obj";

//-----------------------------------------------------------------------------
// Normalized Diffusion.
// [1]  Approximate Reflectance Profiles for Efficient Subsurface Scattering.
//      Per H. Christensen, Brent Burley.
//      http://graphics.pixar.com/library/ApproxBSSRDF/index.html
// [2]  Extending the Disney BRDF to a BSDF with Integrated Subsurface Scattering.
//      Brent Burley. SIG'15 course notes.
//      http://blog.selfshadow.com/publications/s2015-shading-course/burley/s2015_pbs_disney_bsdf_notes.pdf
// [3]  BSSRDF Importance Sampling. Alan King et al. SIG'13 Talks.
//      https://www.solidangle.com/research/s2013_bssrdf_slides.pdf

void    NDProfile::setDistance(const AtVector &dist, const AtColor &albedo)
{
    auto intensity = rls::colorToLuminance(albedo);
    auto s = 1.85f - intensity + 7.0f * powf(intensity - 0.8f, 3.0f);
    mDistance = dist;
    // Scale the max tracing distance empirically.
    mMaxRadius = MAX(mDistance.x, MAX(mDistance.y, mDistance.z)) * 3.0f;
    //mMaxRadius = MAX(dist.x, MAX(dist.y, dist.z));

    for (auto i = 0; i < 3; i++) {
        auto d = mDistance[i];
        mC1[i] = 1.0f - exp(-mMaxRadius / d);
        mC2[i] = 1.0f - exp(-mMaxRadius / d / 3.0f);
    }
}

float   NDProfile::getRadius(float rx) const
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

float   NDProfile::getPdf(float r) const
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

AtRGB   NDProfile::evalProfile(float r) const
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