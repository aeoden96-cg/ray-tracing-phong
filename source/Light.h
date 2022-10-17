#ifndef LIGHT_H
#define LIGHT_H

#include "utils.h"

class Light
{
public:
    Light(const glm::vec3 &p, const glm::vec3 &i) : position(p), intensity(i) {}
    glm::vec3 position;
    glm::vec3 intensity;
};


#endif // LIGHT_H