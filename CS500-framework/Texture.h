#pragma once

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
