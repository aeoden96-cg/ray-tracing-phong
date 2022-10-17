#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include "utils.h"
#include "Hittable.h"

class MeshTriangle : public Hittable
{
public:
    MeshTriangle(
        const glm::vec3 *verts,
        const uint32_t *vertsIndex,
        const uint32_t &numTris,
        const glm::vec2 *st)
    {
        uint32_t maxIndex = 0;
        for (uint32_t i = 0; i < numTris * 3; ++i)
            if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
        maxIndex += 1;
        vertices = std::unique_ptr<glm::vec3[]>(new glm::vec3[maxIndex]);
        memcpy(vertices.get(), verts, sizeof(glm::vec3) * maxIndex);
        vertexIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
        memcpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
        numTriangles = numTris;
        stCoordinates = std::unique_ptr<glm::vec2[]>(new glm::vec2[maxIndex]);
        memcpy(stCoordinates.get(), st, sizeof(glm::vec2) * maxIndex);
    }

    bool intersect(const glm::vec3 &orig, const glm::vec3 &dir, float &tnear, uint32_t &index, glm::vec2 &uv) const
    {
        bool intersect = false;
        for (uint32_t k = 0; k < numTriangles; ++k) {
            const glm::vec3 & v0 = vertices[vertexIndex[k * 3]];
            const glm::vec3 & v1 = vertices[vertexIndex[k * 3 + 1]];
            const glm::vec3 & v2 = vertices[vertexIndex[k * 3 + 2]];
            float t, u, v;
            if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
                tnear = t;
                uv.x = u;
                uv.y = v;
                index = k;
                intersect |= true;
            }
        }

        return intersect;
    }

    void getSurfaceProperties(const glm::vec3 &P, const glm::vec3 &I, const uint32_t &index, const glm::vec2 &uv, glm::vec3 &N, glm::vec2 &st) const
    {
        const glm::vec3 &v0 = vertices[vertexIndex[index * 3]];
        const glm::vec3 &v1 = vertices[vertexIndex[index * 3 + 1]];
        const glm::vec3 &v2 = vertices[vertexIndex[index * 3 + 2]];
        glm::vec3 e0 = normalize(v1 - v0);
        glm::vec3 e1 = normalize(v2 - v1);
        N = normalize(glm::cross(e0, e1));
        const glm::vec2 &st0 = stCoordinates[vertexIndex[index * 3]];
        const glm::vec2 &st1 = stCoordinates[vertexIndex[index * 3 + 1]];
        const glm::vec2 &st2 = stCoordinates[vertexIndex[index * 3 + 2]];
        st = st0 * (1 - uv.x - uv.y) + st1 * uv.x + st2 * uv.y;
    }

    glm::vec3 evalDiffuseColor(const glm::vec2 &st) const
    {
        float scale = 5;
        float pattern = (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
        return mix(glm::vec3(0.815, 0.235, 0.031), glm::vec3(0.937, 0.937, 0.231), pattern);
    }

    std::unique_ptr<glm::vec3[]> vertices;
    uint32_t numTriangles;
    std::unique_ptr<uint32_t[]> vertexIndex;
    std::unique_ptr<glm::vec2[]> stCoordinates;
};




#endif // MESHTRIANGLE_H