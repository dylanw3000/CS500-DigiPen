//////////////////////////////////////////////////////////////////////
// Provides the framework for a raytracer.
////////////////////////////////////////////////////////////////////////

#include <vector>
#define TRANSDEBUG false

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
        if (f.size() >= 11) {
            currentMat = new Material(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), f[7], vec3(f[8], f[9], f[10]), f[11]);
        }
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
    // int index = int(myrandom(RNGen) * vectorOfLights.size());
    Shape* shape = vectorOfLights[0];
    

    if(dynamic_cast<Sphere*>(shape) != nullptr) {
        Sphere* s = static_cast<Sphere*>(shape);
        Q = s->SampleSphere(s->center, s->radius);
    }

    return Q;
}

vec3 SampleLobe(vec3 A, float c, float phi) {
    float s = sqrtf(1.f - c * c);

    vec3 K(s * cosf(phi), s * sinf(phi), c);
    if (fabs(A.z - 1.f) < 10E-3) return K;
    if (fabs(A.z + 1.f) < 10E-3) return vec3(K.x, -K.y, -K.z);

    A = normalize(A);
    vec3 B = normalize(vec3(-A.y, A.x, 0.f));
    vec3 C = cross(A, B);

    return normalize(K.x*B + K.y*C + K.z*A);
}

vec3 Intersection::SampleBrdf(vec3 wo, vec3 N) {
    vec3 Kd = object->mat->Kd;
    vec3 Ks = object->mat->Ks;
    vec3 Kt = object->mat->Kt;
    float pd = object->mat->Pd;
    float pr = object->mat->Pr;
    float pt = object->mat->Pt;

    float E1 = myrandom(RNGen);
    float E2 = myrandom(RNGen);

    float E = myrandom(RNGen);
    vec3 wi; // output
    if (E < pd)
    { // diffuse 
        wi = SampleLobe(N, sqrtf(E1), 2 * PI * E2);

    }
    else
    { // reflection
        float cosThetam = cos(atan((object->mat->alpha * sqrtf(E1)) / sqrtf(1 - E1)));
        vec3 m = SampleLobe(N, cosThetam, 2 * PI * E2);
        float w0dotm = dot(wo, m) > 1.f ? 1.f : dot(wo, m);
        wi = normalize(2 * abs(w0dotm) * m - wo);

        if (E >= pd + pr)
        { // transmission


            float w0dotm = dot(wo, m) > 1.f ? 1.f : dot(wo, m);
            float w0dotN = dot(wo, N) > 1.f ? 1.f : dot(wo, N);

            // Calculate n
            float n = 0.f;
            if (w0dotN > 0)
            { // Path leaves an object, passing from inside to outside
                float ni = 1.f;
                float n0 = object->mat->IoR;
                n = ni / n0;
            }
            else
            { // Path enters an object, passing from outside to inside
                float ni = object->mat->IoR;
                float n0 = 1.f;
                n = ni / n0;

            }
            // Test radicand:
            float r = 1 - pow(n, 2) * (1 - pow(w0dotm, 2));
            if (r >= 0)
                wi = (n * w0dotm - sign(w0dotN) * sqrtf(r)) * m - n * wo;

        }

    }

    return wi;
}



vec3 Intersection::F(float d) {
    return object->mat->Ks + (vec3(1.f) - object->mat->Ks) * powf(1 - abs(d), 5);
}

float Xp(float d) {
    if (d > 0) return 1.f;
    return 0.f;
}

float Intersection::D(vec3 m) {
    float mdotN = dot(m, N);
    float Ndotm = dot(N, m);
    if (mdotN > 1.f)
        mdotN = 1.f;
    if (Ndotm > 1.f)
        Ndotm = 1.f;
    float tanThetaM = sqrtf(1.0f - pow(mdotN, 2)) / mdotN;
    float denominator = PI * pow(Ndotm, 4) * pow(pow(object->mat->alpha, 2) + pow(tanThetaM, 2), 2);
    float out = Xp(mdotN) * (pow(object->mat->alpha, 2) / denominator);
    if (out > 1.f)
        out = 1.f;
    return out;
}

float Intersection::G1(vec3 v, vec3 m) {
    float vdotN = dot(v, N);
    if (vdotN > 1.f)
        return 1.0f;
    float vdotm = dot(v, m);
    if (vdotm > 1.f)
        vdotm = 1.f;

    float tanThetaV = sqrtf(1.0f - pow(vdotN, 2)) / vdotN;
    if (abs(tanThetaV) < 0.0001f)
        return 1.0f;
    float denominator = 1 + sqrtf(1 + pow(object->mat->alpha, 2) * pow(tanThetaV, 2));
    return Xp(vdotm / vdotN) * (2 / denominator);
}

float Intersection::G(vec3 wi, vec3 wo, vec3 m) {
    return G1(wi, m) * G1(wo, m);
}

vec3 Intersection::EvalScattering(vec3 wo, vec3 N, vec3 wi) {
    float pd = object->mat->Pd;
    float pr = object->mat->Pr;
    float pt = object->mat->Pt;


    vec3 Ed = vec3(0.f);
    vec3 Er = vec3(0.f);
    vec3 Et = vec3(0.f);
    vec3 m = normalize(wo + wi);
    float widotm = dot(wi, m) > 1.f ? 1.f : dot(wi, m);
    float widotN = dot(wi, N) > 1.f ? 1.f : dot(wi, N);
    float w0dotN = dot(wo, N) > 1.f ? 1.f : dot(wo, N);
    float Ndotwi = dot(N, wi) > 1.f ? 1.f : dot(N, wi);


    // Calculate Ed
    if (pd > 10E-4)
        Ed = object->mat->Kd / PI;


    // Calculate Er
    float denominator = 4 * abs(widotN) * abs(w0dotN);
    if (pr > 10E-4)
        Er = (D(m) * G(wi, wo, m) * F(widotm)) / denominator;

    float Di, Gi;
    vec3 Fi;


    vec3 numerator_0;
    float denominator_0;
    float numerator_1;
    float denominator_1;
    float n = 0.f, n0 = 0.f, ni = 0.f, Pt = 0.f;
    float w0dotm;

    // Calculate Et
    if (pt > 10E-4)
    {
        if (w0dotN > 0)
        { // Path leaves an object, passing form inside to outside
            ni = 1.f;
            n0 = object->mat->IoR;
            m = -(n0 * wi + ni * wo);
            n = ni / n0;
        }
        else
        { // Path enters an object, passing from outside to inside
            ni = object->mat->IoR;
            n0 = 1.f;
            m = -(n0 * wi + ni * wo);
            n = ni / n0;

        }

        m = normalize(m);
        widotm = dot(wi, m) > 1.f ? 1.f : dot(wi, m);
        w0dotm = dot(wo, m) > 1.f ? 1.f : dot(wo, m);


        // Calculate A(t) for Beer's Law
        vec3 At = vec3(1.f);
        if (w0dotN < 0)
            At = glm::exp(t * glm::log(object->mat->Kt));



        // Test radicand:
        float r = 1 - pow(n, 2) * (1 - pow(w0dotm, 2));
        if (r >= 0)
        {
            Di = D(m);
            Gi = G(wi, wo, m);
            Fi = F(widotm);
            numerator_0 = (D(m) * G(wi, wo, m) * (vec3(1.f) - F(widotm)));
            denominator_0 = abs(widotN) * abs(w0dotN);

            numerator_1 = abs(widotm) * abs(w0dotm) * pow(n0, 2);
            denominator_1 = pow(n0 * widotm + ni * w0dotm, 2);

            Et = At * (numerator_0 / denominator_0) * (numerator_1 / denominator_1);
        }
        else
            Et = At * Er;
    }

    vec3 out = abs(Ndotwi) * (Ed + Er + Et);
    return out;
}



float Intersection::PdfBrdf(vec3 wo, vec3 N, vec3 wi) {
    vec3 m = (wo + wi) / glm::length(wo + wi);
    float pd = object->mat->Pd;
    float pr = object->mat->Pr;
    float pt = object->mat->Pt;
    float Pd = 0.f, Pr = 0.f, Pt = 0.f;

    float ior = object->mat->IoR;

    if (pd > 0.0001f)
        Pd = abs(dot(wi, N)) / PI;

    if (pr > 0.0001f)
        Pr = D(m) * abs(dot(m, N)) * (1 / (4 * abs(dot(wi, m))));


    if (pt > 0.0001f)
    {
        // update m for Pt calculation
        float ni, n0, n;
        if (dot(wo, N) > 0)
        { // Path leaves an object, passing form inside to outside
            ni = 1.f;
            n0 = ior;
            m = -(n0 * wi + ni * wo);
            n = ni / n0;
        }
        else
        { // Path enters an object, passing from outside to inside
            ni = ior;
            n0 = 1.f;
            m = -(n0 * wi + ni * wo);
            n = ni / n0;

        }
        m = normalize(m);


        // Test radicand:
        float r = 1.f - pow(n, 2) * (1.f - pow(dot(wo, m), 2));
        if (r >= 0)
        {
            Pt = D(m) * abs(dot(m, N)) * ((pow(n0, 2) * abs(dot(wi, m))) / pow(n0 * dot(wi, m) + ni * dot(wo, m), 2));
        }
        else
        {
            Pt = Pr;
        }
    }


    return pd * Pd + pr * Pr + pt * Pt;
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
    if (P.object->mat->isLight()) return EvalRadiance(P);

    vec3 wo = -ray.dir;
    vec3 wi = normalize(ray.dir);
    
    // while russian roulette
    float russianRoulette = 0.8f;

    while (myrandom(RNGen) <= russianRoulette) {

        // Explicit Light connection
        Intersection L = SampleLight();
        float p = PdfLight(L) / GeometryFactor(P, L);
        
        wi = normalize(L.P - P.P);
        Intersection I = bvh->intersect(Ray(P.P, wi));
        if (p > 0.f && I.collision && distance(I.P, L.P) < 10E-3) {
            float q = P.PdfBrdf(wo, N, wi) * russianRoulette;
            float wmis = (p * p) / ((p * p) + (q * q));
            // std::cout << p << ", " << q << " = " << wmis << std::endl;
            // wmis = .5;

            vec3 f = P.EvalScattering(wo, N, wi);
            C += W * wmis * f / p * EvalRadiance(L);
        }
        

        // Extend Path
        wi = P.SampleBrdf(wo, N);
        Ray t(P.P, wi);
        Intersection Q = bvh->intersect(t);
        if (!Q.collision) break;

        // std::cout << Q.object->mat->Kd.x << std::endl;
        
        vec3 f = P.EvalScattering(wo, N, wi);
        p = P.PdfBrdf(wo, N, wi) * russianRoulette;
        if (p < 10E-6) break;
        
        W *= f / p;


        // Implicit Light Connection
        if (Q.object->mat->isLight()) {
            float q = PdfLight(Q) / GeometryFactor(P, Q);
            float wmis = (p * p) / ((p * p) + (q * q));
            // wmis = .5;
            C += abs(W * wmis * EvalRadiance(Q));
            break;
        }


        // Step Forward
        P = Q;
        N = P.N;
        wo = -wi;
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

#if TRANSDEBUG
    Ray r(camera->eye, vec3(-0.838028, -0.47094, -0.275543));
    TracePath(r);
    return;
#endif
    

    // for(int p=0; p<pass; p++){
#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
        for (int y = 0; y < height; y++) {

            fprintf(stderr, "Rendering %4d %4d\r", y, pass);
            // fprintf(stderr, "Rendering %4d %4d\r", y, pass);
            for (int x = 0; x < width; x++) {
                Color color = Color(0, 0, 0);

                float dx = 2.f * (float(x) + myrandom(RNGen)) / float(width) - 1;
                float dy = 2.f * (float(y) + myrandom(RNGen)) / float(height) - 1;
                // myrandom(RNGen)

                Intersection front, current;
                front.t = INFINITY;
                Ray ray(camera->eye, glm::normalize(dx * X + dy * Y - Z));

                Color C = TracePath(ray);
                if (all(isnan(C)) || all(isinf(C))) continue;
                image[y * width + x] += C;
            }
        }
    // }
    // fprintf(stderr, "\n");
}
