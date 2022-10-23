#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include <random>
#include "utils.h"
#include "Hittable.h"

class MeshTriangle : public Hittable
{
public:
    MeshTriangle(
        const std::vector<glm::vec3>& verts,
        const std::vector<glm::uvec3>&  vertInd,
        const std::vector<glm::uvec2>&st):
        vertices(verts),
        vertexIndex(vertInd),
        stCoordinates(st){}

//
//    MeshTriangle(
//
//            const glm::mat4 &o2w,
//            const uint32_t nfaces,
//            const std::unique_ptr<uint32_t []> &faceIndex,
//            const std::unique_ptr<uint32_t []> &vertsIndex,
//            const std::unique_ptr<glm::vec3 []> &verts,
//            std::unique_ptr<glm::vec3 []> &normals,
//            std::unique_ptr<glm::vec2 []> &st,
//            bool singleVertAttr = true) :
//            Hittable(o2w),
//            numTris(0),
//            isSingleVertAttr(singleVertAttr)
//    {
//        uint32_t k = 0, maxVertIndex = 0;
//        // find out how many triangles we need to create for this mesh
//        for (uint32_t i = 0; i < nfaces; ++i) {
//            numTris += faceIndex[i] - 2;
//            for (uint32_t j = 0; j < faceIndex[i]; ++j)
//                if (vertsIndex[k + j] > maxVertIndex)
//                    maxVertIndex = vertsIndex[k + j];
//            k += faceIndex[i];
//        }
//        maxVertIndex += 1;
//
//        // allocate memory to store the position of the mesh vertices
//        vertices = std::unique_ptr<glm::vec3 []>(new glm::vec3[maxVertIndex]);
//
//
//        for (uint32_t i = 0; i < maxVertIndex; ++i) {
//            multVecMatrix(verts[i], vertices[i],objectToWorld);
//        }
//
//        // allocate memory to store triangle indices
//        vertexIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
//        glm::mat4 transformNormals = glm::transpose(worldToObject);
//        // [comment]
//        // Sometimes we have 1 vertex attribute per vertex per face. So for example of you have 2
//        // quads this would be defefined by 6 vertices but 2 * 4 vertex attribute values for
//        // each vertex attribute (normal, tex. coordinates, etc.). But in some cases you may
//        // want to have 1 single value per vertex. So in the quad example this would be 6 vertices
//        // and 6 vertex attributes values per attribute. We need to provide both option to users.
//        // [/comment]
//        if (isSingleVertAttr) {
//            N = std::unique_ptr<glm::vec3 []>(new glm::vec3[maxVertIndex]);
//            stCoordinates = std::unique_ptr<glm::vec2 []>(new glm::vec2[maxVertIndex]);
//            for (uint32_t i = 0; i < maxVertIndex; ++i) {
//                stCoordinates[i] = st[i];
//                multDirMatrix(normals[i], N[i],transformNormals);
//            }
//        }
//        else {
//            N = std::unique_ptr<glm::vec3 []>(new glm::vec3[numTris * 3]);
//            stCoordinates = std::unique_ptr<glm::vec2 []>(new glm::vec2[numTris * 3]);
//            for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
//                for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) {
//                    multDirMatrix(normals[k], N[l],transformNormals);
//                    multDirMatrix(normals[k + j + 1], N[l + 1],transformNormals);
//                    multDirMatrix(normals[k + j + 2], N[l + 2],transformNormals);
//                    N[l] = glm::normalize(N[l]);
//                    N[l + 1] = glm::normalize(N[l + 1]);
//                    N[l + 2] = glm::normalize(N[l + 2]);
//                    stCoordinates[l] = st[k];
//                    stCoordinates[l + 1] = st[k + j + 1];
//                    stCoordinates[l + 2] = st[k + j + 2];
//                }
//                k += faceIndex[i];
//            }
//        }
//
//        // generate the triangle index array and set normals and st coordinates
//        for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
//            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) { // for each triangle in the face
//                vertexIndex[l] = vertsIndex[k];
//                vertexIndex[l + 1] = vertsIndex[k + j + 1];
//                vertexIndex[l + 2] = vertsIndex[k + j + 2];
//                l += 3;
//            }
//            k += faceIndex[i];
//        }
//    }

    bool intersect(
        const glm::vec3 &orig,
        const glm::vec3 &dir,
        float &tnear, 
        uint32_t &index, 
        glm::vec2 &uv,
        hit_record& rec
        ) const override;

    void getSurfaceProperties(
        const glm::vec3 &I, 
        hit_record& rec) const override;

    glm::vec3 evalDiffuseColor(const glm::vec2 &st) const;

    const std::vector<glm::vec3> vertices;
    const std::vector<glm::uvec3> vertexIndex;
    const std::vector<glm::uvec2> stCoordinates;

    float numTris;
    std::unique_ptr<glm::vec3 []> N;              // triangles vertex normals


    bool smoothShading = true;       // smooth shading by default
    bool isSingleVertAttr = true;   // single vertex attribute by default


private:
    glm::vec3 gaussianNoise(const glm::vec2 &st) const;
    glm::vec3 checkerPattern(const glm::vec2 &st) const;
};




#endif // MESHTRIANGLE_H