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

