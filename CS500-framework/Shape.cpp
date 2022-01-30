#include "Shape.h"


Intersection Sphere::Intersect(csRay ray) {
	Intersection out;

	vec3 Q = ray.pos - pos;
	float a = dot(ray.dir, ray.dir);
	float b = 2.f * dot(Q, ray.dir);
	float c = dot(Q, Q) - (r * r);

	float discriminant = b * b - 4 * a * c;
	if (discriminant >= 0) {
		discriminant = sqrtf(discriminant / 4.f);
		float left = -b / 2.f;
		float t = min(left - discriminant, left + discriminant);
		out.addIntersect(this, t, ray.eval(t), normalize(ray.eval(t)-pos));
	}

	return out;
}



Interval slabInterval(csRay ray, float d0, float d1, vec3 N) {
	float t0 = -(d0 + dot(N, ray.pos)) / dot(N, ray.dir);
	float t1 = -(d1 + dot(N, ray.pos)) / dot(N, ray.dir);
	return Interval(t0, t1, -N, N);
}

Intersection Box::Intersect(csRay ray) {
	Intersection out;

	Interval i1, i2, i3;

	i1 = slabInterval(ray, -b.x, -b.x - d.x, { 1,0,0 });
	i2 = slabInterval(ray, -b.y, -b.y - d.y, { 0,1,0 });
	i3 = slabInterval(ray, -b.z, -b.z - d.z, { 0,0,1 });

	/*
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
	*/

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



Intersection Cylinder::Intersect(csRay ray) {
	// Intersection tmp;  return tmp;
	Intersection out;
	// return out;

	vec3 A = normalize(ang);
	bool zAxis = (A == vec3(0, 0, 1));
	vec3 B = normalize(cross({ 0,0,1 }, A));
	vec3 C = cross(A, B);

	mat3 Rinv(B, C, A);
	mat3 R = glm::transpose(Rinv);

	vec3 Q = R * (ray.pos - base);
	vec3 D = R * ray.dir;

	// zAxis = false;
	if (zAxis) {
		Q = ray.pos - base;
		D = ray.dir;
	}

	csRay tRay(Q, D);

	Interval i1 = slabInterval(tRay, 0, -length(ang), { 0,0,1 });

	/*
	{
		vec3 N(0, 0, 1);
		float d0 = 0;
		float d1 = -length(ang);

		float t0 = -(d0 + dot(N, ray.pos)) / dot(N, ray.dir);
		float t1 = -(d1 + dot(N, ray.pos)) / dot(N, ray.dir);
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



Intersection Triangle::Intersect(csRay ray) {
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

	out.addIntersect(this, t, ray.eval(t), (1.f-u-v)*n0 + u*n1 + v*n2);
	return out;
}
