#ifndef HITTABLE_H
#define HITTABLE_H

#include "utils.h"

class Hittable
{
 public:
    Hittable() :
        materialType(DIFFUSE_AND_GLOSSY),
        ior(1.3), Kd(0.8), Ks(0.2), diffuseColor(0.2), specularExponent(25) {}
    virtual ~Hittable() {}
    virtual bool intersect(const glm::vec3 &, const glm::vec3 &, float &, uint32_t &, glm::vec2 &) const = 0;
    virtual void getSurfaceProperties(const glm::vec3 &, const glm::vec3 &, const uint32_t &, const glm::vec2 &, glm::vec3 &, glm::vec2 &) const = 0;
    virtual glm::vec3 evalDiffuseColor(const glm::vec2 &) const { return diffuseColor; }
    // material properties
    MaterialType materialType;
    float ior;
    float Kd, Ks;
    glm::vec3 diffuseColor;
    float specularExponent;
};







#endif // HITTABLE_H