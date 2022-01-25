#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

class Camera {
public:
    Camera() {}

public:
	glm::vec3 eye;
    glm::quat orient;   // Represents rotation of -Z to view direction
    float ry;
    float front, back;
    float spin, tilt;
    float cDist;              // Distance from eye to center of scene

};
