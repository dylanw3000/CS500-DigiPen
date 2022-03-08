
#if false

#include "Intersection.h"
#include "Shape.h"
#include "raytrace.h"
#include "Random.h"

using namespace glm;

#define PI 3.14159f

vec3 SampleLobe(vec3 A, float c, float phi) {
    float s = sqrtf(1.f - c * c);

    vec3 K(s * cosf(phi), s * sinf(phi), c);
    if (fabs(A.z - 1.f) < 10E-3) return K;
    if (fabs(A.z + 1.f) < 10E-3) return K;

    A = normalize(A);
    vec3 B = normalize(vec3(-A.y, A.x, 0.f));
    vec3 C = cross(A, B);

    return K.x * B + K.y * C + K.z * A;
}



vec3 Intersection::SampleBrdf(vec3 wo, vec3 N) {
    float m1 = myrandom(RNGen);
    float m2 = myrandom(RNGen);
    if (myrandom(RNGen) < object->mat->Pd) {
        return SampleLobe(N, sqrtf(m1), 2.f * PI * m2);
    }
    vec3 m = SampleLobe(N, cos(powf(m1, 1.f / (object->mat->alpha + 1))), 2.f * PI * m2);
    return 2.f * fabs(dot(wo, m)) * m - wo;
}



vec3 Intersection::F(float d) {
    return object->mat->Ks + (vec3(1.f) - object->mat->Ks) * powf(1 - abs(d), 5);
}

float Xp(float d) {
    if (d > 0) return 1;
    return 0;
}

float Intersection::D(vec3 m) {
    return Xp(dot(m, N)) * (object->mat->alpha + 2) / (2 * PI) * powf(dot(m, N), object->mat->alpha);
}

float Intersection::G1(vec3 v, vec3 m) {
    float d = dot(v, N);
    if (d > 1.f) return 1;

    float theta = sqrtf(1 - pow(d, 2)) / d;
    if (fabs(theta) == 0) return 1;

    float x = Xp(dot(v, m) / dot(v, N));
    float a = sqrtf(object->mat->alpha / 2 + 1) / theta;
    if (a < 1.6) {
        return x * (3.535 * a + 2.181 * a * a) / (1 + 2.276 * a + 2.577 * a * a);
    }
    return x;
}

float Intersection::G(vec3 wi, vec3 wo, vec3 m) {
    return G1(wi, m) * G1(wo, m);
}

vec3 Intersection::EvalScattering(vec3 wo, vec3 N, vec3 wi) {
    vec3 m = normalize(wo + wi);
    vec3 Ed = object->mat->Kd / PI;
    vec3 Er = D(m) * G(wi, wo, m) * F(dot(wi, m)) / (4 * abs(dot(wi, N)) * abs(pow(wo, N)));
    return abs(dot(N, wi)) * (Ed + Er);
}



float Intersection::PdfBrdf(vec3 wo, vec3 N, vec3 wi) {
    vec3 m = normalize(wo + wi);
    float pd = fabs(dot(wi, N)) / PI;
    float pr = D(m) * fabs(dot(m, N)) / (4 * abs(dot(wi, m)));
    return pd * object->mat->Pd + pr * object->mat->Pr;
    // return length(dot(N, wi)) / PI;
}

#endif
