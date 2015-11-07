#ifndef RLS_UTIL_H
#define RLS_UTIL_H

#include <vector>
#include <functional>
#include <memory>

#include <ai.h>

#include "ext/tinyexr.h"

namespace rls
{

inline AtVector sphericalDirection(float cosTheta, float phi)
{
    AtVector omega;
    omega.z = cosTheta;
    float r = sqrtf(1.0f - SQR(omega.z));
    omega.x = r * cosf(phi);
    omega.y = r * sinf(phi);
    return omega;
}

inline float colorToLuminance(const AtRGB &col)
{
    return col.r * 0.212671f + col.g * 0.715160f + col.b * 0.072169f;
}

AtVector   concentricDiskSample(float rx, float ry);

//! Image writer for sampling pattern.
class   SampleWriter
{
public:
    static const int kChannelNum = 3;

public:
    SampleWriter(int w, int h, const std::string &outpath = "")
        : mWidth(w), mHeight(h), mImagePath(outpath)
    {
        mData.resize(w * h * kChannelNum);
        mStride = w * h;

        mImage.width = w;
        mImage.height = h;
        mImage.num_channels = kChannelNum;

        static const char *channel_names[] = { "B", "G", "R" };
        mImage.channel_names = channel_names;

        static int pixel_types[kChannelNum] = { TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT };
        static int request_pixel_types[kChannelNum] = { TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF };

        mImage.pixel_types = pixel_types;
        mImage.requested_pixel_types = request_pixel_types;

        if (mImagePath.empty()) {
            // Set the default output path to the same directory of beauty pass.
            auto driver = AiNodeLookUpByName("defaultArnoldDriver@driver_exr.RGBA");
            if (driver) {
                mImagePath = AiNodeGetStr(driver, "filename");
                auto pos = mImagePath.find_last_of("/");
                mImagePath.erase(pos);
                mImagePath.insert(pos, "/rls_sampling_pattern.exr");
            }
        }
    }

    ~SampleWriter()
    {
        auto data = std::make_unique<unsigned char *[]>(3);
        data[0] = reinterpret_cast<unsigned char *>(&mData[0]);
        data[1] = reinterpret_cast<unsigned char *>(&mData[0 + mStride]);
        data[2] = reinterpret_cast<unsigned char *>(&mData[0 + mStride * 2]);
        mImage.images = data.get();

        const char *err;
        if (SaveMultiChannelEXRToFile(&mImage, mImagePath.c_str(), &err) == 0) {
            AiMsgInfo("[RLS] Write sampling pattern to %s", mImagePath.c_str());
        } else {
            AiMsgWarning("[RLS] Failed to dump the image of sampling pattern! %s", err);
        }
    }

public:
    template <typename Brdf>
    void    writeRadiance(AtShaderGlobals *sg, Brdf &brdf)
    {
        using namespace std::placeholders;
        auto evalBrdf = std::bind(Brdf::evalBrdf, &brdf, _1);

        for (auto j = 0; j < mHeight; j++) {
            auto theta = AI_PIOVER2 * j / mHeight;

            for (auto i = 0; i < mWidth; i++) {
                auto phi = AI_PITIMES2 * i / mWidth;
                AtVector dir = sphericalDirection(cosf(theta), phi);
                auto color = evalBrdf(&dir);
                writePixel(i, j, color);
            }
        }
    }

    template <typename Brdf>
    void    writeSample(AtShaderGlobals *sg, Brdf &brdf, int sampleCount)
    {
        using namespace std::placeholders;
        auto evalSample = std::bind(Brdf::evalSample, &brdf, _1, _2);

        auto sampler = AiSampler(sampleCount, 2);
        auto sampleIter = AiSamplerIterator(sampler, sg);
        float samples[2];

        int missingCount = 0;

        while (AiSamplerGetSample(sampleIter, samples)) {
            auto dir = evalSample(samples[0], samples[1]);
            if (AiV3IsZero(dir)) {
                continue;
            }

            float theta = acosf(dir.z);
            float phi = atan2f(dir.y, dir.x);

            if (phi < 0.0f) {
                phi += AI_PITIMES2;
            }

            int i = CLAMP(static_cast<int>(phi * AI_ONEOVER2PI * mWidth), 0, mWidth - 1);
            int j = CLAMP(static_cast<int>(theta / AI_PIOVER2 * mHeight), 0, mHeight - 1);

            if (theta > AI_PIOVER2) {
                writePixel(i, j, AI_RGB_RED);
                AiMsgInfo("[RLS] Theta diff %f!", (theta - AI_PIOVER2) * AI_ONEOVERPI * 180.0f);
                missingCount++;
            } else {
                writePixel(i, j, AI_RGB_GREEN);
            }
        }

        auto totalCount = AiSamplerGetSampleCount(sampleIter);
        AiMsgInfo("[RLS] Missing %d ray samples out of %d!", missingCount, totalCount);
        AiSamplerDestroy(sampler);
    }

    inline void writePixel(int x, int y, AtColor color)
    {
        mData[x + y * mWidth] = color.b;
        mData[x + y * mWidth + mStride] = color.g;
        mData[x + y * mWidth + mStride * 2] = color.r;
    }

private:
    EXRImage            mImage;
    std::vector<float>  mData;
    std::string         mImagePath;
    int                 mWidth, mHeight;
    int                 mStride;
};

}
#endif
