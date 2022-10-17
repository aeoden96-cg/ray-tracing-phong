#include "utils.h"



float clamp_(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

glm::vec3 mix(const glm::vec3 &a, const glm::vec3& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0) return false;
    else if (discr == 0) x0 = x1 = - 0.5 * b / a;
    else {
        float q = (b > 0) ?
            -0.5 * (b + sqrt(discr)) :
            -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);
    return true;
}

bool rayTriangleIntersect(
    const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2,
    const glm::vec3 &orig, const glm::vec3 &dir,
    float &tnear, float &u, float &v)
{
    glm::vec3 edge1 = v1 - v0;
    glm::vec3 edge2 = v2 - v0;
    glm::vec3 pvec = glm::cross(dir, edge2);
    float det = glm::dot(edge1, pvec);
    if (det == 0 || det < 0) return false;

    glm::vec3 tvec = orig - v0;
    u = glm::dot(tvec, pvec);
    if (u < 0 || u > det) return false;

    glm::vec3 qvec = glm::cross(tvec, edge1);
    v = glm::dot(dir, qvec);
    if (v < 0 || u + v > det) return false;

    float invDet = 1 / det;
    
    tnear = glm::dot(edge2, qvec) * invDet;
    u *= invDet;
    v *= invDet;

    return true;
}

// Compute reflection direction
glm::vec3 reflect(const glm::vec3 &I, const glm::vec3 &N)
{
    return I - 2 * glm::dot(I, N) * N;
}


glm::vec3 refract(const glm::vec3 &I, const glm::vec3 &N, const float &ior)
{
    float cosi = clamp_(-1, 1, glm::dot(I, N));
    float etai = 1, etat = ior;
    glm::vec3 n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? glm::vec3(0,0,0) : eta * I + (eta * cosi - sqrtf(k)) * n;
}

// Compute Fresnel equation
//
// \param I is the incident view direction
//
// \param N is the normal at the intersection point
//
// \param ior is the mateural refractive index
//
// \param[out] kr is the amount of light reflected
void fresnel(const glm::vec3 &I, const glm::vec3 &N, const float &ior, float &kr)
{
    float cosi = clamp_(-1, 1, glm::dot(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}