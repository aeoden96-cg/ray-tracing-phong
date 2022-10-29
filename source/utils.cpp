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
float fresnel(const glm::vec3 &I, const glm::vec3 &N, const float &ior)
{
    float kr = 0;
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
    return kr;
}



glm::vec3 evalBezierCurve(const std::vector<glm::vec3> &P, const float &t, int offset)
{
    float b0 = (1 - t) * (1 - t) * (1 - t);
    float b1 = 3 * t * (1 - t) * (1 - t);
    float b2 = 3 * t * t * (1 - t);
    float b3 = t * t * t;

    return  P[offset] * b0 +
            P[offset + 1] * b1 +
            P[offset + 2] * b2 +
            P[offset + 3] * b3;
}

glm::vec3 dVBezier(const std::vector<glm::vec3> &controlPoints, const float &u, const float &v)
{
    std::vector<glm::vec3> uCurve(4);
    for (int i = 0; i < 4; ++i) {
        uCurve[i] = evalBezierCurve(controlPoints, u, i * 4);
    }

    return derivBezier(uCurve, v);
}

glm::vec3 derivBezier(const std::vector<glm::vec3>& P, const float &t)
{
    return -3 * (1 - t) * (1 - t) * P[0] +
           (3 * (1 - t) * (1 - t) - 6 * t * (1 - t)) * P[1] +
           (6 * t * (1 - t) - 3 * t * t) * P[2] +
           3 * t * t * P[3];
}

void multDirMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4& x)
{
    float a, b, c;

    a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0];
    b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1];
    c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2];

    dst.x = a;
    dst.y = b;
    dst.z = c;
}

void multVecMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4& x)
{
    float a, b, c, w;

    a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
    b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
    c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
    w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];

    dst.x = a / w;
    dst.y = b / w;
    dst.z = c / w;
}



glm::vec3 evalBezierPatch(const std::vector<glm::vec3>& controlPoints, const float &u, const float &v)
{
    std::vector<glm::vec3> uCurve(4);
    for (int i = 0; i < 4; ++i)
        uCurve[i] = evalBezierCurve(controlPoints, u, i * 4);

    return evalBezierCurve(uCurve, v, 0);
}


glm::vec3 dUBezier(const std::vector<glm::vec3> &controlPoints, const float &u, const float &v)
{
    std::vector<glm::vec3> P(4);
    std::vector<glm::vec3> vCurve(4);
    for (int i = 0; i < 4; ++i) {
        P[0] = controlPoints[i];
        P[1] = controlPoints[4 + i];
        P[2] = controlPoints[8 + i];
        P[3] = controlPoints[12 + i];
        vCurve[i] = evalBezierCurve(P, v, 0);
    }

    return derivBezier(vCurve, u);
}