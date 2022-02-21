//////////////////////////////////////////////////////////////////////
// Provides the framework for a raytracer.
////////////////////////////////////////////////////////////////////////

#include <vector>

#ifdef _WIN32
    // Includes for Windows
#include <windows.h>
#include <cstdlib>
#include <limits>
#include <crtdbg.h>
#else
    // Includes for Linux
#endif

#include "geom.h"
#include "raytrace.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Random.h"

Scene::Scene()
{
    camera = new Camera();
}

void Scene::Finit()
{
    std::cout << "vectorOfShapes size: " << vectorOfShapes.size() << std::endl;

    for (auto shape : vectorOfShapes) {
        if (shape->mat->isLight()) {
            vectorOfLights.push_back(shape);
        }
    }
    std::cout << "vectorOfLights size: " << vectorOfLights.size() << std::endl;

    bvh = new AccelerationBvh(vectorOfShapes);
}

void Scene::triangleMesh(MeshData* mesh)
{
    // mesh->vertices[0].pnt;
    for (auto tri : mesh->triangles) {
        Triangle* t = new Triangle(
            mesh->vertices[tri[0]].pnt,
            mesh->vertices[tri[1]].pnt,
            mesh->vertices[tri[2]].pnt,
            mesh->vertices[tri[0]].nrm,
            mesh->vertices[tri[1]].nrm,
            mesh->vertices[tri[2]].nrm,
            currentMat
        );
        vectorOfShapes.push_back(t);
    }
}

quat Orientation(int i,
    const std::vector<std::string>& strings,
    const std::vector<float>& f)
{
    quat q(1, 0, 0, 0); // Unit quaternion
    while (i < strings.size()) {
        std::string c = strings[i++];
        if (c == "x")
            q *= angleAxis(f[i++] * Radians, Xaxis());
        else if (c == "y")
            q *= angleAxis(f[i++] * Radians, Yaxis());
        else if (c == "z")
            q *= angleAxis(f[i++] * Radians, Zaxis());
        else if (c == "q") {
            q *= quat(f[i + 0], f[i + 1], f[i + 2], f[i + 3]);
            i += 4;
        }
        else if (c == "a") {
            q *= angleAxis(f[i + 0] * Radians, normalize(vec3(f[i + 1], f[i + 2], f[i + 3])));
            i += 4;
        }
    }
    return q;
}

void Scene::Command(const std::vector<std::string>& strings,
    const std::vector<float>& f)
{
    if (strings.size() == 0) return;
    std::string c = strings[0];

    if (c == "screen") {
        // syntax: screen width height
        width = int(f[1]);
        height = int(f[2]);
    }

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
        camera->eye = vec3(f[1], f[2], f[3]);
        camera->orient = Orientation(5, strings, f);
        camera->ry = f[4];
    }

    else if (c == "ambient") {
        // syntax: ambient radius g b
        // Sets the ambient color.  Note: This parameter is temporary.
        // It will be ignored once your raytracer becomes capable of
        // accurately *calculating* the true ambient light.
        ambient = vec3(f[1], f[2], f[3]);
    }

    else if (c == "brdf") {
        // syntax: brdf  radius g b   radius g b  alpha
        // later:  brdf  radius g b   radius g b  alpha  radius g b ior
        // First rgb is Diffuse reflection, second is specular reflection.
        // third is beer's law transmission followed by index of refraction.
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Material(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), f[7]);
    }

    else if (c == "light") {
        // syntax: light  radius g b   
        // The rgb is the emission of the light
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Light(vec3(f[1], f[2], f[3]));
    }

    else if (c == "sphere") {
        // syntax: sphere x y z   radius
        // Creates a Shape instance for a sphere defined by a center and radius

        Sphere* s = new Sphere(vec3(f[1], f[2], f[3]), f[4], currentMat);
        vectorOfShapes.push_back(s);
    }

    else if (c == "box") {
        // syntax: box bx by bz   dx dy dz
        // Creates a Shape instance for a box defined by a corner point and diagonal vector

        Box* b = new Box(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), currentMat);
        vectorOfShapes.push_back(b);
    }

    else if (c == "cylinder") {
        // syntax: cylinder bx by bz   ax ay az  radius
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius

        Cylinder* c = new Cylinder(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), f[7], currentMat);
        vectorOfShapes.push_back(c);
    }


    else if (c == "mesh") {
        // syntax: mesh   filename   tx ty tz   s   <orientation>
        // Creates many Shape instances (one per triangle) by reading
        // model(s) from filename. All triangles are rotated by a
        // quaternion (qw qx qy qz), uniformly scaled by s, and
        // translated by (tx ty tz) .
        mat4 modelTr = translate(vec3(f[2], f[3], f[4]))
            * scale(vec3(f[5], f[5], f[5]))
            * toMat4(Orientation(6, strings, f));
        ReadAssimpFile(strings[1], modelTr);
    }


    else {
        fprintf(stderr, "\n*********************************************\n");
        fprintf(stderr, "* Unknown command: %s\n", c.c_str());
        fprintf(stderr, "*********************************************\n\n");
    }
}



Intersection Sphere::SampleSphere(vec3 C, float R) {
    float m1 = myrandom(RNGen);
    float m2 = myrandom(RNGen);
    float z = 2.f * m1 - 1.f;
    float r = sqrtf(1.f - z * z);
    float a = 2 * PI * m2;

    Intersection out;
    vec3 N = normalize(vec3(r * cosf(a), r * sinf(a), z));
    out.addIntersect(this, 0.f, C + R * N, N);
    return out;
}

Intersection Scene::SampleLight() {
    Intersection Q;
    if (vectorOfLights.size() == 0) return Q;

    // TODO select light randomly (uniformly) ; currently hardcoded
    Shape* shape = vectorOfLights[0];
    

    if(dynamic_cast<Sphere*>(shape) != nullptr) {
        Sphere* s = static_cast<Sphere*>(shape);
        Q = s->SampleSphere(s->pos, s->radius);
    }

    return Q;
}

vec3 SampleLobe(vec3 A, float c, float phi) {
    float s = sqrtf(1.f - c * c);

    vec3 K(s * cosf(phi), s * sinf(phi), c);
    if (fabs(A.z - 1.f) < 10E-3) return K;
    if (fabs(A.z + 1.f) < 10E-3) return K;

    A = normalize(A);
    vec3 B = normalize(vec3(-A.y, A.x, 0.f));
    vec3 C = cross(A, B);

    return K.x*B + K.y*C + K.z*A;
}

vec3 SampleBrdf(vec3 N) {
    float m1 = myrandom(RNGen);
    float m2 = myrandom(RNGen);
    return SampleLobe(N, sqrtf(m1), 2.f * PI * m2);
}

vec3 Intersection::EvalScattering(vec3 N, vec3 wi) {
    return fabs(dot(N, wi)) * object->mat->Kd / PI;
}

float PdfBrdf(vec3 N, vec3 wi) {
    return length(dot(N, wi)) / PI;
}

float Scene::PdfLight(Intersection L) {
    if (dynamic_cast<Sphere*>(L.object) != nullptr) {
        Sphere* s = static_cast<Sphere*>(L.object);
        return 1.f / (4.f * PI * powf(s->radius, 2) * vectorOfLights.size());
    }
    return 0.f;
}

vec3 EvalRadiance(Intersection Q) {
    if(Q.object->mat->isLight()) return Q.object->mat->Kd;
    return vec3(0.f);
}

float GeometryFactor(Intersection A, Intersection B) {
    vec3 D = A.P - B.P;
    return fabsf( dot(A.N, D) * dot(B.N, D) / pow2(dot(D, D)));
}



vec3 Scene::TracePath(Ray ray) {
    vec3 C(0, 0, 0);
    vec3 W(1, 1, 1);

    Intersection P = bvh->intersect(ray);
    vec3 N = P.N;
    
    if (!P.collision) return C;
    if (P.object->mat->isLight()) return P.object->mat->Kd;     // TODO evalRadiance(P)
    
    // while russian roulette
    float russianRoulette = 0.8f;

    while (myrandom(RNGen) <= russianRoulette) {
        N = P.N;
        vec3 wi = SampleBrdf(N);


        Intersection L = SampleLight();
        float p = PdfLight(L) / GeometryFactor(P, L);
        Intersection I = bvh->intersect(Ray(P.P, normalize(L.P - P.P)));
        if (p > 0.f && I.collision && distance(I.P, L.P) < 10E-3) {
            vec3 f = P.EvalScattering(N, wi);
            C += 0.5f * W * f / p * EvalRadiance(L);
        }
        

        
        Ray t(P.P, wi);
        Intersection Q = bvh->intersect(t);
        if (!Q.collision) break;
        
        vec3 f = P.EvalScattering(N, wi);
        p = PdfBrdf(N, wi) * russianRoulette;
        if (p < 10E-6) break;

        W *= f / p;

        if (Q.object->mat->isLight()) {
            C += 0.5f * W * EvalRadiance(Q);
            break;
        }

        P = Q;
    }

    return C;
}



void Scene::TraceImage(Color* image, const int pass)
{
    float rx = camera->ry * float(width) / float(height);

    glm::vec3 X, Y, Z;
    X = rx * transformVector(camera->orient, Xaxis());
    Y = camera->ry * transformVector(camera->orient, Yaxis());
    Z = transformVector(camera->orient, Zaxis());

    // for(int p=0; p<pass; p++){
#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
        for (int y = 0; y < height; y++) {

            fprintf(stderr, "Rendering %4d %4d\r", y, pass);
            for (int x = 0; x < width; x++) {
                Color color = Color(0, 0, 0);

                float dx = 2.f * (float(x) + myrandom(RNGen)) / float(width) - 1;
                float dy = 2.f * (float(y) + myrandom(RNGen)) / float(height) - 1;
                // myrandom(RNGen)

                Intersection front, current;
                front.t = INFINITY;
                Ray ray(camera->eye, glm::normalize(dx * X + dy * Y - Z));

                image[y * width + x] += TracePath(ray);

#if false
                /*
                for (auto shape : vectorOfShapes) {
                    current = shape->Intersect(ray);
                    if (current.collision && (!front.collision || current.t < front.t)) {
                        front = current;
                    }
                }
                */

                front = bvh->intersect(ray);

                if (front.object != nullptr) {
                    /*
                    float t = (front.t - 5.f) / 4.f;
                    color = Color(t, t, t);
                    color = Color(front.object->mat->Kd);
                    color = Color(abs(front.N));
                    */

                    vec3 lightPos(1.9, 5, 2);
                    vec3 lightColor(3, 3, 3);

                    vec3 N = normalize(front.N);
                    vec3 E = normalize(camera->eye - front.P);

                    vec3 L = normalize(lightPos - front.P);
                    vec3 H = normalize(L + E);

                    float NL = max(0.1f, dot(N, L));
                    float HN = max(0.f, dot(H, N));

                    float alpha = front.object->mat->alpha;
                    // color = Color(lightColor * (front.object->mat->Kd * NL / 3.14f + alpha / 6.28f * front.object->mat->Ks * pow(HN, alpha)));
                    color = Color(lightColor * (front.object->mat->Kd * NL / 3.14f));

                }

                image[y * width + x] += color;
#endif
            }
        }
    // }
    // fprintf(stderr, "\n");
}
