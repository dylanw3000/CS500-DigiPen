#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

class Ray {
public:
	Ray(glm::vec3 _pos, glm::vec3 _dir) : origin(_pos), dir(_dir) {}

	glm::vec3 origin;	// Q
	glm::vec3 dir;	// D

	glm::vec3 eval(float t) { return origin + t * dir; }
	/*
	glm::quat a = glm::angleAxis(3.14f, glm::vec3(1, 0, 0));
	glm::vec3 b = glm::eulerAngles(a);
	*/
};
