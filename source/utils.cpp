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