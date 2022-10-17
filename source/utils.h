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



float clamp_(const float &lo, const float &hi, const float &v);

glm::vec3 mix(const glm::vec3 &a, const glm::vec3& b, const float &mixValue);

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1);







#endif // UTILS_H