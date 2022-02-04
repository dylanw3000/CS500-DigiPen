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

// A good quality *thread-safe* Mersenne Twister random number generator.
#include <random>
std::random_device device;
std::mt19937_64 RNGen(device());
std::uniform_real_distribution<> myrandom(0.0, 1.0);
// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].

Scene::Scene()
{
    camera = new Camera();
}

void Scene::Finit()
{
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

void Scene::TraceImage(Color* image, const int pass)
{
    // realtime->run();                          // Remove this (realtime stuff)

    std::cout << "vectorOfShapes size: " << vectorOfShapes.size() << std::endl;

    float rx = camera->ry * float(width) / float(height);

    glm::vec3 X, Y, Z;
    X = rx * transformVector(camera->orient, Xaxis());
    Y = camera->ry * transformVector(camera->orient, Yaxis());
    Z = transformVector(camera->orient, Zaxis());

#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
    for (int y = 0; y < height; y++) {

        fprintf(stderr, "Rendering %4d\r", y);
        for (int x = 0; x < width; x++) {
            Color color = Color(0, 0, 0);

            float dx = 2.f * (float(x) + 0.5f) / float(width) - 1;
            float dy = 2.f * (float(y) + 0.5f) / float(height) - 1;

            Intersection front, current;
            front.t = INFINITY;
            Ray ray(camera->eye, glm::normalize(dx * X + dy * Y - Z));

            /*
            for (auto shape : vectorOfShapes) {
                current = shape->Intersect(ray);
                if (current.collision && (!front.collision || current.t < front.t)) {
                    front = current;
                }
            }
            */

            front = bvh->intersect(ray);


            /*
            if ((x-width/2)*(x-width/2)+(y-height/2)*(y-height/2) < 100*100)
                color = Color(myrandom(RNGen), myrandom(RNGen), myrandom(RNGen));
            else if (abs(x-width/2)<4 || abs(y-height/2)<4)
                color = Color(0.0, 0.0, 0.0);
            else
                color = Color(1.0, 1.0, 1.0);
            */

            if (front.object != nullptr) {
                float t = (front.t - 5.f) / 4.f;
                color = Color(t, t, t);
                color = Color(front.object->mat->Kd);

                
                color = Color(abs(front.N));
                vec3 lightPos(1.9, 5, 2);
                vec3 lightColor(3, 3, 3);

                vec3 N = normalize(front.N);
                vec3 E = normalize(camera->eye - front.P);

                vec3 L = normalize(lightPos - front.P);
                vec3 H = normalize(L + E);

                float NL = max(0.1f, dot(N, L));
                float HN = max(0.f, dot(H, N));

                float alpha = front.object->mat->alpha;
                color = Color(lightColor * (front.object->mat->Kd * NL / 3.14f + alpha / 6.28f * front.object->mat->Ks * pow(HN, alpha)));
                
            }

            image[y * width + x] = color;
        }
    }
    fprintf(stderr, "\n");
}
