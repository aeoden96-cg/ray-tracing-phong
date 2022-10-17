#ifndef SPHERE_H
#define SPHERE_H

#include "utils.h"
#include "Hittable.h"


class Sphere : public Hittable{
public:
    Sphere(const glm::vec3 &c, const float &r) : center(c), radius(r), radius2(r * r) {}
    bool intersect(const glm::vec3 &orig, const glm::vec3 &dir, float &tnear, uint32_t &index, glm::vec2 &uv,hit_record& rec) const
    {
        // analytic solution
        glm::vec3 L = orig - center;
        float a = glm::dot(dir, dir);
        float b = 2 * glm::dot(dir, L);
        float c = glm::dot(L, L) - radius2;
        float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) return false;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        tnear = t0;

        return true;
    }

    void getSurfaceProperties(const glm::vec3 &I, hit_record& rec) const
    { rec.normal = normalize(rec.p - center); }

    glm::vec3 center;
    float radius, radius2;
};







#endif // SPHERE_H