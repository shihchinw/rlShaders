#include "rlUtil.h"

AtVector   rls::concentricDiskSample(float rx, float ry)
{
    rx = rx * 2.0f - 1.0f;
    ry = ry * 2.0f - 1.0f;

    AtVector result;
    if (rx == 0.0f && ry == 0.0f) {
        // Handle degenerate case
        result.x = result.y = 0.0f;
        return result;
    }

    float r, phi;
    if (ABS(rx) > ABS(ry)) {
        r = rx;
        phi = AI_PIOVER2 * 0.5f * ry / rx;
    } else {
        r = ry;
        phi = AI_PIOVER2 * (1.0f - 0.5f * rx / ry);
    }

    result.x = r * cosf(phi);
    result.y = r * sinf(phi);
    return result;
}