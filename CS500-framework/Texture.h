#pragma once

#include <stb_image.h>

class _Image {
public:
    int height, width, channels;
    uint8_t* img;

    _Image(const char* filename) {
        stbi_set_flip_vertically_on_load(true);
        img = stbi_load(filename, &width, &height, &channels, 3);
        if (img == NULL) {
            std::cout << "Image failed to load" << std::endl;
            exit(1);
        }
    }

    vec3 getUV(float u, float v) {
        // Intentionally incomplete, not rounding up or linearly interpolating color
        u = clamp(u, 0.f, 1.f);
        int x = (width-1) * u;
        v = clamp(v, 0.f, 1.f);
        int y = (height-1) * v;

        int index = channels * (x + y * width);
        return vec3(img[index], img[index + 1], img[index + 2]) / 256.f;
    }

    vec3 getUV(vec2 uv) { return getUV(uv.x, uv.y); }
};



#if false

class Texture
{
public:
    unsigned int textureId;
    int width, height, depth;
    unsigned char* image;
    Texture();
    Texture(const std::string& filename);

    void BindTexture(const int unit, const int programId, const std::string& name);
    void UnbindTexture(const int unit);
};

#endif
