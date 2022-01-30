#pragma once

#include <glm/glm.hpp>
using namespace glm;

class Material;
class Shape;


class Interval {
public:
	Interval() : t0(0.f), t1(INFINITY) {}
	Interval(float _t0, float _t1, glm::vec3 _n0, glm::vec3 _n1) {
		if (_t0 <= _t1) {
			t0 = _t0;
			t1 = _t1;
			N0 = _n0;
			N1 = _n1;
		}
		else {
			t0 = _t1;
			t1 = _t0;
			N0 = _n1;
			N1 = _n0;
		}
	}

public:
	void empty() {
		t0 = 0.f;
		t1 = -1.f;
	}

	void intersect() {
		//
	}

public:
	float t0, t1;
	glm::vec3 N0, N1;
};

class Intersection {
public:
	Intersection() : t(INFINITY), object(nullptr) {}

	void addIntersect(Shape* _o, float _t, vec3 _P, vec3 _N) {
		collision = true;
		t = _t;
		object = _o;
		P = _P;
		N = _N;
	}

	float t;
	Shape* object;
	bool collision = false;
	glm::vec3 P;
	glm::vec3 N;
};

#include "Ray.h"

/*
struct Intersection {
	float t;
	Shape* object;
	glm::vec3 P;
	glm::vec3 N;
};

// #include "Ray.h"
*/


class Shape {
public:
	Shape() {}
	virtual Intersection Intersect(csRay ray) = 0;
	Material* mat;
};



class Sphere : public Shape {
public:
	Sphere(glm::vec3 _pos, float _r, Material* _mat) : pos(_pos), r(_r) {
		mat = _mat;
	}
	// Sphere();

	Intersection Intersect(csRay ray);

	glm::vec3 pos;
	float r;
	// Material* mat;
};



class Box : public Shape {
public:
	Box(glm::vec3 _v1, glm::vec3 _v2, Material* _mat) : b(_v1), d(_v2) {
		mat = _mat;
	}
	Intersection Intersect(csRay ray);

	glm::vec3 b;
	glm::vec3 d;
	// Material* mat;

};



class Cylinder : public Shape {
public:
	Cylinder(glm::vec3 _pos, glm::vec3 _ang, float _r, Material* _mat) : base(_pos), ang(_ang), radius(_r) {
		mat = _mat;
	}
	Intersection Intersect(csRay ray);

	glm::vec3 base;
	glm::vec3 ang;
	float radius;
	// Material* mat;
};



class Triangle : public Shape {
public:
	Triangle(glm::vec3 _v0, glm::vec3 _v1, glm::vec3 _v2, Material* _mat) : v0(_v0), v1(_v1), v2(_v2) {
		mat = _mat;
	}

	Triangle(vec3 _v0, vec3 _v1, vec3 _v2, vec3 _n0, vec3 _n1, vec3 _n2, Material* _mat) : v0(_v0), v1(_v1), v2(_v2), n0(_n0), n1(_n1), n2(_n2) {
		mat = _mat;
	}

public:
	Intersection Intersect(csRay ray);

public:
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;

	glm::vec3 n0;
	glm::vec3 n1;
	glm::vec3 n2;
};









