#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

const float kInfinity = std::numeric_limits<float>::max();

enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION };

/*struct Material
{
    MaterialType type;
    glm::vec3 color;
    float ior; // index of refraction
    float diffuse;
    float specular;
    float roughness;
};*/


struct Options
{
    uint32_t width;
    uint32_t height;
    float fov;
    float imageAspectRatio;
    uint8_t maxDepth;
    glm::vec3 backgroundColor;
    float bias;
};

float clamp_(const float &lo, const float &hi, const float &v);

glm::vec3 mix(const glm::vec3 &a, const glm::vec3& b, const float &mixValue);

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1);

bool rayTriangleIntersect(
    const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2,
    const glm::vec3 &orig, const glm::vec3 &dir,
    float &tnear, float &u, float &v);


// Compute reflection direction
glm::vec3 reflect(const glm::vec3 &I, const glm::vec3 &N);

// Compute refraction direction using Snell's law
//
// We need to handle with care the two possible situations:
//
//    - When the ray is inside the Hittable
//
//    - When the ray is outside.
//
// If the ray is outside, you need to make cosi positive cosi = -N.I
//
// If the ray is inside, you need to invert the refractive indices and negate the normal N
glm::vec3 refract(const glm::vec3 &I, const glm::vec3 &N, const float &ior);


// Compute Fresnel equation
//
// \param I is the incident view direction
//
// \param N is the normal at the intersection point
//
// \param ior is the mateural refractive index
//
// \param[out] kr is the amount of light reflected
float fresnel(const glm::vec3 &I, const glm::vec3 &N, const float &ior);

#endif // UTILS_H