#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

class csRay {
public:
	csRay(glm::vec3 _pos, glm::vec3 _dir) : pos(_pos), dir(_dir) {}

	glm::vec3 pos;	// Q
	glm::vec3 dir;	// D

	glm::vec3 eval(float t) { return pos + t * dir; }
	/*
	glm::quat a = glm::angleAxis(3.14f, glm::vec3(1, 0, 0));
	glm::vec3 b = glm::eulerAngles(a);
	*/
};
