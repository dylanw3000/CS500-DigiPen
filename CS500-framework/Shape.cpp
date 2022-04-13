#include "Shape.h"

#define BOUNDINGBOX true
#define PI 3.14159f

Intersection Sphere::Intersect(Ray ray) {
	Intersection out;

	vec3 Q = ray.origin - center;

	float t0 = -dot(Q, ray.dir) + pow(pow(dot(Q, ray.dir), 2.f) - dot(Q, Q) + pow(radius, 2.f), 0.5f);
	float t1 = -dot(Q, ray.dir) - pow(pow(dot(Q, ray.dir), 2.f) - dot(Q, Q) + pow(radius, 2.f), 0.5f);

	float t = -1.f;
	if (t0 >= 10E-3 && t1 >= 10E-3) {
		t = min(t0, t1);
	}
	else if (t0 >= 10E-3) {
		t = t0;
	}
	else if (t1 > 10E-3) {
		t = t1;
	}
	else {
		return out;
	}

	vec3 pos = ray.eval(t);
	vec3 normal = glm::normalize(pos - center);
	float theta = atan2(normal.y, normal.x);
	float psi = acos(normal.z);
	vec2 uv = vec2(theta / (2 * PI), psi / PI);

	out.addIntersect(this, t, pos, normal, uv);
	return out;
}





Interval slabInterval(Ray ray, float d0, float d1, vec3 N) {
	float denom = dot(N, ray.dir);
	if (abs(denom) < 0.001f) {
		Interval i;
		i.empty();
		return i;
	}
	float t0 = -(d0 + dot(N, ray.origin)) / denom;
	float t1 = -(d1 + dot(N, ray.origin)) / denom;
	return Interval(t0, t1, -N, N);
}

Intersection Box::Intersect(Ray ray) {
	Intersection out;

	if (false) {
		out.addIntersect(this, 2, ray.eval(2), vec3(0, 1, 0));
		return out;
	}

	Interval i1, i2, i3;

	i1 = slabInterval(ray, -b.x, -b.x - d.x, { 1,0,0 });
	i2 = slabInterval(ray, -b.y, -b.y - d.y, { 0,1,0 });
	i3 = slabInterval(ray, -b.z, -b.z - d.z, { 0,0,1 });

	if (i1.t0 > i1.t1 || i2.t0 > i2.t1 || i3.t0 > i3.t1) {
		return out;
	}

	Interval* tmin = &i1;
	if (i2.t0 > tmin->t0) tmin = &i2;
	if (i3.t0 > tmin->t0) tmin = &i3;

	Interval* tmax = &i1;
	if (i2.t1 < tmax->t1) tmax = &i2;
	if (i3.t1 < tmax->t1) tmax = &i3;

	/*
	float tmin = max(i1.t0, max(i2.t0, i3.t0));
	float tmax = min(i1.t1, min(i2.t1, i3.t1));
	*/

	if (tmin->t0 < tmax->t1) {
		if (tmin->t0 > 0.f) {
			out.addIntersect(this, tmin->t0, ray.eval(tmin->t0), tmin->N0);
		}
		else if (tmax->t1 > 0.f) {
			out.addIntersect(this, tmax->t1, ray.eval(tmax->t1), tmax->N1);
		}
	}

	return out;
}



Intersection Cylinder::Intersect(Ray ray) {
	Intersection out;
	
	if (false) {
		out.addIntersect(this, 2.f, ray.eval(2.f), vec3(0, 1, 0));
		return out;
	}

	vec3 A = normalize(ang);
	bool zAxis = (A == vec3(0, 0, 1));
	vec3 B = normalize(cross({ 0,0,1 }, A));
	vec3 C = cross(A, B);

	mat3 Rinv(B, C, A);
	mat3 R = glm::transpose(Rinv);

	vec3 Q = R * (ray.origin - base);
	vec3 D = R * ray.dir;

	// zAxis = false;
	if (zAxis) {
		Q = ray.origin - base;
		D = ray.dir;
	}

	Ray tRay(Q, D);

	Interval i1 = slabInterval(tRay, 0, -length(ang), { 0,0,1 });
	if (i1.t0 > i1.t1) {
		return out;
	}

	/*
	{
		vec3 N(0, 0, 1);
		float d0 = 0;
		float d1 = -length(ang);

		float t0 = -(d0 + dot(N, ray.origin)) / dot(N, ray.dir);
		float t1 = -(d1 + dot(N, ray.origin)) / dot(N, ray.dir);
		i1 = Interval(t0, t1, N, -N);
	}
	*/

	float a = (D.x*D.x + D.y*D.y);
	float b = 2.f * (D.x*Q.x + D.y*Q.y);
	float c = Q.x*Q.x + Q.y*Q.y - radius*radius;

	float discriminant = b*b - 4*a*c;
	if (discriminant < 0) { return out; }

	Interval i2;
	{
		float t0 = (-b + sqrtf(discriminant)) / (2 * a);
		float t1 = (-b - sqrtf(discriminant)) / (2 * a);

		float t = t0<t1 ? t0 : t1;
		vec3 N = normalize(vec3(Q.x + t * D.x, Q.y + t * D.y, 0));
		i2 = Interval(t0, t1, -N, N);
	}

	float tmin = max(i1.t0, i2.t0);
	float tmax = min(i1.t1, i2.t1);

	if (tmin < tmax) {
		if (tmin > 0) {
			Interval* imin = i1.t0 < i1.t0 ? &i1 : &i2;
			out.addIntersect(this, tmin, ray.eval(tmin), zAxis ? imin->N0 : Rinv*imin->N0);
		}
		else if (tmax > 0) {
			vec3 N = vec3(Q.x + tmax * D.x, Q.y + tmax * D.y, 0);
			out.addIntersect(this, tmax, ray.eval(tmax), zAxis ? N : Rinv*N);
		}
	}

	return out;
}



Intersection Triangle::Intersect(Ray ray) {
	Intersection out;


	glm::vec3 e1 = v1 - v0;
	glm::vec3 e2 = v2 - v0;
	glm::vec3 p = glm::cross(ray.dir, e2);
	float d = glm::dot(p, e1);

	if (d == 0) { return out; }

	glm::vec3 s = ray.origin - v0;
	float u = glm::dot(p, s) / d;

	if (u < 0.f || u > 1.f) { return out; }

	vec3 q = cross(s, e1);
	float v = dot(ray.dir, q) / d;

	if (v < 0.f || u + v > 1.f) { return out; }

	float t = dot(e2, q) / d;
	if (t < 0.f) { return out; }

	out.addIntersect(this, t, ray.eval(t), (1.f-u-v)*n0 + u*n1 + v*n2);
	return out;
}
