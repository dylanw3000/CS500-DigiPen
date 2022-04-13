///////////////////////////////////////////////////////////////////////
// A framework for a raytracer.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include "Camera.h"

#include "Shape.h"
#include "Ray.h"
#include "acceleration.h"

#include "Texture.h"

const float PI = 3.14159f;
const float Radians = PI / 180.0f;    // Convert degrees to radians

////////////////////////////////////////////////////////////////////////
// Material: encapsulates a BRDF and communication with a shader.
////////////////////////////////////////////////////////////////////////
class Material
{
public:
    vec3 Kd, Ks, Kt;
    float alpha;

    float Pd, Pr, Pt;
    float IoR = 1.f;

    virtual bool isLight() { return false; }

    bool isTexture = false;
    _Image* tex;

    Material() : Kd(vec3(1.0, 0.5, 0.0)), Ks(vec3(1, 1, 1)), alpha(1.0) { Pd = length(Kd) / (length(Kd) + length(Ks)); Pr = length(Ks) / (length(Kd) + length(Ks)); }
    Material(const vec3 Kd_, const vec3 Ks_, const float a_)
        : Kd(Kd_), Ks(Ks_), Kt(vec3(0.f)), alpha(a_) {
        // s = length(Kd) + length(Ks) + length(Kt); 
        Pd = length(Kd) / (length(Kd) + length(Ks)); 
        Pr = length(Ks) / (length(Kd) + length(Ks));
        Pt = 0.f;
        // std::cout << "Kd: " << length(Kd) << ", Ks: " << length(Ks) << ", Pd: " << Pd << std::endl;
    }
    Material(Material& o) { 
        Kd = o.Kd;  
        Ks = o.Ks;  
        alpha = o.alpha;  
        Pd = length(Kd) / (length(Kd) + length(Ks)); 
        Pr = length(Ks) / (length(Kd) + length(Ks)); 
        // std::cout << Pd << std::endl; 
    }

    // transmissive
    Material(const vec3 Kd_, const vec3 Ks_, const float a_, const vec3 Kt_, const float IoR_)
        : Kd(Kd_), Ks(Ks_), Kt(Kt_), alpha(a_), IoR(IoR_)
    {
        float s = length(Kd) + length(Ks) + length(Kt); 
        Pd = length(Kd) / s;
        Pr = length(Ks) / s;
        Pt = length(Kt) / s;
        // std::cout << Pd << ", " << Pr << ", " << Pt << std::endl;
    }

    //virtual void apply(const unsigned int program);
};

////////////////////////////////////////////////////////////////////////
// Data structures for storing meshes -- mostly used for model files
// read in via ASSIMP.
//
// A MeshData holds two lists (stl::vector) one for vertices
// (VertexData: consisting of point, normal, texture, and tangent
// vectors), and one for triangles (ivec3: consisting of three
// indices into the vertex array).

class VertexData
{
public:
    vec3 pnt;
    vec3 nrm;
    vec2 tex;
    vec3 tan;
    VertexData(const vec3& p, const vec3& n, const vec2& t, const vec3& a)
        : pnt(p), nrm(n), tex(t), tan(a)
    {}
};

struct MeshData
{
    std::vector<VertexData> vertices;
    std::vector<ivec3> triangles;
    Material* mat;
};

////////////////////////////////////////////////////////////////////////
// Light: encapsulates a light and communiction with a shader.
////////////////////////////////////////////////////////////////////////
class Light : public Material
{
public:

    Light(const vec3 e) : Material() { Kd = e; }
    virtual bool isLight() { return true; }
    //virtual void apply(const unsigned int program);
};

class LightTexture : public Light {
public:
    LightTexture(std::string path) : Light(vec3(1.f)) {
        isTexture = true;
        tex = new _Image(path.c_str());
    }
};

////////////////////////////////////////////////////////////////////////////////
// Scene

class Scene {
public:
    int width, height;
    Material* currentMat;

    Camera* camera;
    std::vector<Shape*> vectorOfShapes;
    std::vector<Shape*> vectorOfLights;

    vec3 ambient;

    AccelerationBvh* bvh;

    Scene();
    void Finit();

    // The scene reader-parser will call the Command method with the
    // contents of each line in the scene file.
    void Command(const std::vector<std::string>& strings,
        const std::vector<float>& f);

    // To read a model file into the scene via ASSIMP, call ReadAssimpFile.  
    void ReadAssimpFile(const std::string& path, const mat4& M);

    // Once ReadAssimpFile parses the information from the model file,
    // it will call:
    void triangleMesh(MeshData* mesh);

    // The main program will call the TraceImage method to generate
    // and return the texImg.  This is the Ray Tracer!
    void TraceImage(Color* image, const int pass);
    vec3 TracePath(Ray ray);

    // Intersection SampleSphere(vec3 C, float R);
    Intersection SampleLight();
    float PdfLight(Intersection Q);

    _Image* texImg;
};
