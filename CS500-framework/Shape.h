#pragma once

#include <glm/glm.hpp>

class Shape {
public:
	Shape() {}
};

class Sphere : public Shape {
public:
	Sphere();

	glm::vec3 pos;
	float r;
};

class Box : public Shape {
public:
	Box();

	glm::vec3 pos;
};

class Cylinder : public Shape {
public:
	Cylinder();

	glm::vec3 pos;
};

class Triangle : public Shape {
public:
	Triangle();

	glm::vec3 v1;
	glm::vec3 v2;
	glm::vec3 v3;
};
