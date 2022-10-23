#include "MeshTriangle.h"

bool MeshTriangle::intersect(
        const glm::vec3 &orig,
        const glm::vec3 &dir,
        float &tnear,
        uint32_t &index,
        glm::vec2 &uv,
        hit_record& rec
) const
{
    int k = 0;
    bool intersect = false;
    for (auto triangle : vertexIndex){
        const glm::vec3 & v0 = vertices[triangle.x];
        const glm::vec3 & v1 = vertices[triangle.y];
        const glm::vec3 & v2 = vertices[triangle.z];
        float t, u, v;
        if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
            tnear = t;
            uv.x = u;
            uv.y = v;
            index = k;
            intersect |= true;
        }
        k++;
    }

    return intersect;
}


void MeshTriangle::getSurfaceProperties(
        const glm::vec3 &I,
        hit_record& rec) const
{
    //first we need to get the triangle from the index

    auto& triangle = vertexIndex[rec.triIndex];

    //--------------------normal--------------------

    //now we need to get the vertices of the triangle

    const glm::vec3 & v0 = vertices[triangle.x];
    const glm::vec3 & v1 = vertices[triangle.y];
    const glm::vec3 & v2 = vertices[triangle.z];

    //now we need to get the st coordinates of the triangle
    //st coordinates are the texture coordinates
    //we need to interpolate the st coordinates of the triangle
    //based on the uv coordinates of the intersection point

    //e0 is the vector from v0 to v1
    glm::vec3 e0 = normalize(v1 - v0);
    //e1 is the vector from v0 to v2
    glm::vec3 e1 = normalize(v2 - v1);

    //now we need to get the normal of the triangle using the cross product of e0 and e1
    //this is the same as the normal of the plane that the triangle lies on
    rec.normal = normalize(glm::cross(e0, e1));

    //--------------------st coordinates-------------------

    const glm::vec2 &st0 = stCoordinates[triangle.x];
    const glm::vec2 &st1 = stCoordinates[triangle.y];
    const glm::vec2 &st2 = stCoordinates[triangle.z];

    //now we need to interpolate the st coordinates of the triangle
    //based on the uv coordinates of the intersection point
    //we can do this by using the barycentric coordinates of the intersection point
    //we can get the barycentric coordinates of the intersection point by using the uv coordinate
    rec.st = st0 * (1 - rec.uv.x - rec.uv.y) + st1 * rec.uv.x + st2 * rec.uv.y;

}

glm::vec3 MeshTriangle::evalDiffuseColor(const glm::vec2 &st) const
{

    return checkerPattern(st);
//        return gaussianNoise(st);

}


glm::vec3 MeshTriangle::gaussianNoise(const glm::vec2 &st) const{
    const double mean = 0.0;
    const double stddev = 0.1;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);
    return glm::vec3(dist(generator), dist(generator), dist(generator));
}

glm::vec3 MeshTriangle::checkerPattern(const glm::vec2 &st) const
{
    float scale =5;
//        float sines = sin(10 * st.x) * sin(10 * st.y);
//        if (sines < 0) {
//            return glm::vec3(0.0f);
//        }
//        else {
//            return glm::vec3(1.0f);
//        }
    //pattern is a checkerboard pattern, so we need to get the integer part of the st coordinates
    //and then we need to check if the sum of the integer parts is even or odd
    //if the sum is even, then we return white, otherwise we return black
    //we need integer part of the st coordinates because we want the pattern to repeat
    //every integer part of the st coordinates
    int s = (int)floor(st.x * scale);
    int t = (int)floor(st.y * scale);
    if ((s + t) % 2 == 0) {
        return glm::vec3(0.815, 0.235, 0.031);
    }
    else {
        return glm::vec3(0.937, 0.937, 0.231);
    }

    //        float pattern = (fmodf(st.x * scale, 1) > 0.5f) ^ (fmodf(st.y * scale, 1) > 0.5);
//        return mix(
//                glm::vec3(0.815, 0.235, 0.031),
//                glm::vec3(0.937, 0.937, 0.231),
//                pattern);
}
