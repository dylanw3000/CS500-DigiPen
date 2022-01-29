#pragma once

#include <glm/glm.hpp>

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
	Intersection() : t(INFINITY), object(nullptr), collision(false) {}

	void addIntersect(float _t, Shape* _o) {
		collision = true;
		t = _t;
		object = _o;
	}

	float t;
	Shape* object;
	bool collision;
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

	Intersection Intersect(csRay ray) { 
		Intersection out;

		vec3 Q = ray.pos - pos;
		float a = dot(ray.dir, ray.dir);
		float b = 2.f * dot(Q, ray.dir);
		float c = dot(Q, Q) - (r * r);

		float discriminant = b * b - 4 * a * c;
		if (discriminant >= 0) {
			discriminant = sqrt(discriminant / 4.f);
			float left = -b / 2.f;
			out.addIntersect(std::min(left-discriminant, left+discriminant), this);
		}

		return out;
		/*
		vec3 oc = ray.pos - pos;
		float a = dot(ray.dir, ray.dir);
		float b = 2.0 * dot(oc, ray.dir);
		float c = dot(oc, oc) - radius * radius;
		float discriminant = b * b - 4 * a * c;

		Intersection out;
		if (discriminant >= 0) {
			out.addIntersect((-b - sqrt(discriminant)) / (2.0 * a), this);
		}
		return out;
		*/
	}

	glm::vec3 pos;
	float r;
	// Material* mat;
};

class Box : public Shape {
public:
	Box(glm::vec3 _v1, glm::vec3 _v2, Material* _mat) : b(_v1), d(_v2) {
		mat = _mat;
	}
	Intersection Intersect(csRay ray) { 
		Intersection out;

		Interval i1, i2, i3;

		{
			vec3 N(1, 0, 0);
			float t0 = -(-b.x + dot(N, ray.pos)) / dot(N, ray.dir);
			float t1 = -((-b.x + -d.x) + dot(N, ray.pos)) / dot(N, ray.dir);
			i1 = Interval(t0, t1, N, -N);
		}

		{
			vec3 N(0, 1, 0);
			float t0 = -(-b.y + dot(N, ray.pos)) / dot(N, ray.dir);
			float t1 = -((-b.y + -d.y) + dot(N, ray.pos)) / dot(N, ray.dir);
			i2 = Interval(t0, t1, N, -N);
		}

		{
			vec3 N(0, 0, 1);
			float t0 = -(-b.z + dot(N, ray.pos)) / dot(N, ray.dir);
			float t1 = -((-b.z + -d.z) + dot(N, ray.pos)) / dot(N, ray.dir);
			i3 = Interval(t0, t1, N, -N);
		}

		float tmin = std::max(i1.t0, std::max(i2.t0, i3.t0));
		float tmax = std::min(i1.t1, std::min(i2.t1, i3.t1));

		if (tmin < tmax) {
			if (tmin > 0.f) {
				out.addIntersect(tmin, this);
			}
			else if (tmax > 0.f) {
				out.addIntersect(tmax, this);
			}
		}

		return out;
	}

	glm::vec3 b;
	glm::vec3 d;
	// Material* mat;

};

class Cylinder : public Shape {
public:
	Cylinder(glm::vec3 _pos, glm::vec3 _ang, float _r, Material* _mat) : base(_pos), ang(_ang), radius(_r) {
		mat = _mat;
	}
	Intersection Intersect(csRay ray) {
		// Intersection tmp;  return tmp;
		Intersection out;
		return out;

		vec3 A = normalize(ang);
		vec3 B = normalize(cross({ 0,0,1 }, A));
		vec3 C = cross(A, B);

		mat3 Rinv(B, C, A);
		mat3 R(glm::transpose(Rinv));

		vec3 Q = ray.pos - base;
		vec3 D = R * ray.dir;

		Interval i1;
		{
			vec3 N(0, 0, 1);
			float d0 = 0;
			float d1 = -length(A);

			float t0 = -(d0 + dot(N, ray.pos)) / dot(N, ray.dir);
			float t1 = -(d1 + dot(N, ray.pos)) / dot(N, ray.dir);
			i1 = Interval(t0, t1, N, -N);
		}

		float a = (D.x * D.x + D.y + D.y);
		float b = 2.f * (D.x * Q.x + D.y + Q.y);
		float c = Q.x * Q.x + Q.y * Q.y - (radius * radius);

		float discriminant = b * b - 4 * a * c;
		if (discriminant < 0) { return out; }

		Interval i2;
		{
			float t0 = (-b - sqrtf(discriminant)) / (2 * a);
			float t1 = (-b + sqrtf(discriminant)) / (2 * a);
			Interval i2(t0, t1, { 0,0,1 }, { 0,0,-1 });
		}

		float tmin = std::max(i1.t0, i2.t0);
		float tmax = std::min(i1.t1, i2.t1);

		if (tmin < tmax) {
			if (tmin > 0) {
				out.addIntersect(tmin, this);
			}
			else if (tmax > 0) {
				out.addIntersect(tmax, this);
			}
		}

		return out;
		
		/*
		vec3 pos = ray.pos;
		vec3 dir = ray.dir;
		vec3 center = base;
		float radius = radius;

		float a = (dir.x * dir.x) + (dir.z * dir.z);
		float b = 2 * (dir.x * (pos.x - center.x) + dir.z * (pos.z - center.z));
		float c = (pos.x - center.x) * (pos.x - center.x) + (pos.z - center.z) * (pos.z - center.z) - (radius * radius);

		float delta = b * b - 4 * (a * c);
		if (fabs(delta) < 0.001) return out;
		if (delta < 0.0) return out;

		float t1 = (-b - sqrt(delta)) / (2 * a);
		float t2 = (-b + sqrt(delta)) / (2 * a);
		float t;

		if (t1 > t2) t = t2;
		else t = t1;

		float r = pos.y + t * dir.y;

		if ((r >= center.y) && (r <= center.y + length(ang))) out.addIntersect(t, this);
		
		return out;
		*/
		
	}

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

public:
	Intersection Intersect(csRay ray) {
		Intersection out;

		
		glm::vec3 e1 = v1 - v0;
		glm::vec3 e2 = v2 - v0;
		glm::vec3 p = glm::cross(ray.dir, e2);
		float d = glm::dot(p, e1);

		if (d == 0) { return out; }
		
		glm::vec3 s = ray.pos - v0;
		float u = glm::dot(p, s) / d;

		if (u < 0.f || u > 1.f) { return out; }

		vec3 q = cross(s, e1);
		float v = dot(ray.dir, q) / d;

		if (v < 0.f || u + v > 1.f) { return out; }

		float t = dot(e2, q) / d;
		if (t < 0.f) { return out; }

		out.addIntersect(t, this);
		return out;
		

		/*
		vec3 v0v1 = v1 - v0;
		vec3 v0v2 = v2 - v0;
		vec3 pvec = cross(ray.dir, v0v2); // dir.crossProduct(v0v2);
		float det = dot(v0v1, pvec);

		// ray and triangle are parallel if det is close to 0
		if (fabs(det) <= 0.f) return out;
 
		float invDet = 1 / det;

		vec3 tvec = ray.pos - v0;
		float u = dot(tvec, pvec) * invDet;
		if (u < 0.f || u > 1.f) return out;

		vec3 qvec = cross(tvec, v0v1);
		float v = dot(ray.dir, qvec) * invDet;
		if (v < 0 || u + v > 1) return out;

		float t = dot(v0v2, qvec) * invDet;
		// out.addIntersect(t, this);
		return out;
		*/
	}

public:
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;

	glm::vec3 n0;
	glm::vec3 n1;
	glm::vec3 n2;
};









