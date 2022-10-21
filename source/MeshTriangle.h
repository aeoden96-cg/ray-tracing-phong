#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include "utils.h"
#include "Hittable.h"

class MeshTriangle : public Hittable
{
public:
    MeshTriangle(
        const std::vector<glm::vec3>& verts,
        const std::vector<glm::uvec3>&  vertInd,
        const std::vector<glm::uvec2>&st): vertices(verts), vertexIndex(vertInd), stCoordinates(st)
    {
        numTriangles = vertInd.size();

    }

    bool intersect(
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

    void getSurfaceProperties(
        const glm::vec3 &I, 
        hit_record& rec) const
    {

        const glm::vec3 & v0 = vertices[vertexIndex[rec.triIndex].x];
        const glm::vec3 & v1 = vertices[vertexIndex[rec.triIndex].y];
        const glm::vec3 & v2 = vertices[vertexIndex[rec.triIndex].z];

        glm::vec3 e0 = normalize(v1 - v0);
        glm::vec3 e1 = normalize(v2 - v1);
        rec.normal = normalize(glm::cross(e0, e1));

        const glm::vec2 &st0 = stCoordinates[vertexIndex[rec.triIndex].x];
        const glm::vec2 &st1 = stCoordinates[vertexIndex[rec.triIndex].y];
        const glm::vec2 &st2 = stCoordinates[vertexIndex[rec.triIndex].z];
      
        rec.st = st0 * (1 - rec.uv.x - rec.uv.y) + st1 * rec.uv.x + st2 * rec.uv.y;

        return ;
    }

    glm::vec3 evalDiffuseColor(const glm::vec2 &st) const
    {
        float scale = 5;
        float pattern = (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
        return mix(glm::vec3(0.815, 0.235, 0.031), glm::vec3(0.937, 0.937, 0.231), pattern);
    }

    const std::vector<glm::vec3> vertices;
    const std::vector<glm::uvec3> vertexIndex;
    uint32_t numTriangles;
    const std::vector<glm::uvec2> stCoordinates;
};




#endif // MESHTRIANGLE_H