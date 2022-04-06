#pragma once

#include <glm/glm.hpp>

class Shape;


class Intersection {
public:
	Intersection() : t(INFINITY), object(nullptr) {}

	void addIntersect(Shape* _o, float _t, glm::vec3 _P, glm::vec3 _N, glm::vec2 _uv = glm::vec2(0,0)) {
		collision = true;
		t = _t;
		object = _o;
		P = _P;
		N = _N;
		uv = _uv;
	}

	float distance() const { return t; }

	float t;
	Shape* object;
	bool collision = false;
	glm::vec3 P;
	glm::vec3 N;
	bool isLight = false;

	vec3 EvalScattering(vec3 wo, vec3 N, vec3 wi);
	vec3 SampleBrdf(vec3 wo, vec3 N);
	float PdfBrdf(vec3 wo, vec3 N, vec3 wi);

	vec3 F(float d);
	float D(vec3 m);
	float Intersection::G(vec3 wi, vec3 wo, vec3 m);
	float Intersection::G1(vec3 v, vec3 m);

	glm::vec2 uv;
};
