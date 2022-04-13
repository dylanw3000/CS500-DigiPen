#pragma once

#include <stb_image.h>
#include <glm/glm.hpp>
#include <iostream>

using namespace glm;

class _Image {
public:
    int height, width, channels;
    uint8_t* img;

    _Image(const char* filename);

    vec3 getUV(float u, float v);

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
