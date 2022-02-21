#pragma once

#include <glm/glm.hpp>

class Shape;


class Intersection {
public:
	Intersection() : t(INFINITY), object(nullptr) {}

	void addIntersect(Shape* _o, float _t, glm::vec3 _P, glm::vec3 _N) {
		collision = true;
		t = _t;
		object = _o;
		P = _P;
		N = _N;
	}

	float distance() const { return t; }

	float t;
	Shape* object;
	bool collision = false;
	glm::vec3 P;
	glm::vec3 N;
	bool isLight = false;

	vec3 EvalScattering(vec3 N, vec3 wi);
};
