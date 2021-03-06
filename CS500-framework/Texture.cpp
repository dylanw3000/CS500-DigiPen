
#include "Texture.h"

_Image::_Image(const char* filename) {
    // stbi_set_flip_vertically_on_load(true);
    img = stbi_load(filename, &width, &height, &channels, 3);
    if (img == NULL) {
        std::cout << "Image failed to load" << std::endl;
        exit(1);
    }
}

vec3 _Image::getUV(float u, float v) {
    u = clamp(u, 0.f, 1.f);
    int x = (width - 1) * u;
    v = clamp(v, 0.f, 1.f);
    int y = (height - 1) * v;

    int index = channels * (x + y * width);
    if (index < 0) index = 0;
    // if (index >= (height-1) * width) index = 0;

    vec3 out(1.f);
    out.x = img[index];
    out.y = img[index+1];
    out.z = img[index+2];
    return out / 255.f;
}

#if false

#include <glbinding/gl/gl.h>
#include <glbinding/Binding.h>
using namespace gl;

// #define GLM_FORCE_RADIANS
// #define GLM_SWIZZLE
#include <glm/glm.hpp>

#include "texture.h"

// #define STB_IMAGE_IMPLEMENTATION
// #define STBI_FAILURE_USERMSG
#include "stb_image.h"

#include <glu.h>                // For gluErrorString
#define CHECKERROR {GLenum err = glGetError(); if (err != GL_NO_ERROR) { fprintf(stderr, "OpenGL error (at line texture.cpp:%d): %s\n", __LINE__, gluErrorString(err)); exit(-1);} }

Texture::Texture(const std::string& path) : textureId(0)
{
    stbi_set_flip_vertically_on_load(true);
    image = stbi_load(path.c_str(), &width, &height, &depth, 4);
    printf("%d %d %d %s\n", depth, width, height, path.c_str());
    if (!image) {
        printf("\nRead error on file %s:\n  %s\n\n", path.c_str(), stbi_failure_reason());
        exit(-1);
    }

    glGenTextures(1, &textureId);   // Get an integer id for this texture from OpenGL
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D, 0, (GLint)GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 10);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (int)GL_LINEAR);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (int)GL_LINEAR_MIPMAP_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    stbi_image_free(image);
}

#endif
