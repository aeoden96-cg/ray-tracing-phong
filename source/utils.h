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

enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION , METAL};

//    typedef std::vector<PointList> PatchList;
typedef std::vector<glm::vec3> PointList;
typedef std::vector<unsigned > IndexList;
typedef std::vector<glm::vec2> TexSTList;

typedef std::unique_ptr<std::vector<glm::vec3>> UniqPointList;
typedef std::unique_ptr<std::vector<unsigned>> UniqIndexList;
typedef std::unique_ptr<std::vector<glm::vec2>> UniqTexList;

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
    glm::mat4 cameraToWorld;
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

void multVecMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4& x);

void multDirMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4& x);

// Compute the position of a point along a Bézier curve at t [0:1]
glm::vec3 evalBezierCurve(const std::vector<glm::vec3> &P, const float &t, int offset);

glm::vec3 evalBezierPatch(const std::vector<glm::vec3> &controlPoints, const float &u, const float &v);
glm::vec3 derivBezier(const std::vector<glm::vec3>&P, const float &t);

// Compute the derivative of a point on Bezier patch along the u parametric direction
glm::vec3 dUBezier(const std::vector<glm::vec3>& controlPoints, const float &u, const float &v);

// Compute the derivative of a point on Bezier patch along the v parametric direction
glm::vec3 dVBezier(const std::vector<glm::vec3>& controlPoints, const float &u, const float &v);




#endif // UTILS_H