//
// Created by mateo on 23.10.22..
//

#ifndef PROJECT_PATCHMESH_H
#define PROJECT_PATCHMESH_H


#include "utils.h"
#include "Hittable.h"


class PatchMesh : public Hittable
{
public:
    // Build a triangle mesh from a face index array and a vertex index array
    PatchMesh(
            const glm::mat4 &o2w,
            const uint32_t nfaces,
            const std::unique_ptr<uint32_t []> &faceIndex,
            const std::unique_ptr<uint32_t []> &vertsIndex,
            const std::unique_ptr<glm::vec3 []> &verts,
            std::unique_ptr<glm::vec3 []> &normals,
            std::unique_ptr<glm::vec2 []> &st,
            bool singleVertAttr = true) :
            Hittable(o2w),
            numTris(0),
            isSingleVertAttr(singleVertAttr)
    {
        uint32_t k = 0, maxVertIndex = 0;
        // find out how many triangles we need to create for this mesh
        for (uint32_t i = 0; i < nfaces; ++i) {
            numTris += faceIndex[i] - 2;
            for (uint32_t j = 0; j < faceIndex[i]; ++j)
                if (vertsIndex[k + j] > maxVertIndex)
                    maxVertIndex = vertsIndex[k + j];
            k += faceIndex[i];
        }
        maxVertIndex += 1;

        // allocate memory to store the position of the mesh vertices
        P = std::unique_ptr<glm::vec3 []>(new glm::vec3[maxVertIndex]);


        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            multVecMatrix(verts[i], P[i],objectToWorld);
        }

        // allocate memory to store triangle indices
        trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
        glm::mat4 transformNormals = glm::transpose(worldToObject);
        // [comment]
        // Sometimes we have 1 vertex attribute per vertex per face. So for example of you have 2
        // quads this would be defefined by 6 vertices but 2 * 4 vertex attribute values for
        // each vertex attribute (normal, tex. coordinates, etc.). But in some cases you may
        // want to have 1 single value per vertex. So in the quad example this would be 6 vertices
        // and 6 vertex attributes values per attribute. We need to provide both option to users.
        // [/comment]
        if (isSingleVertAttr) {
            N = std::unique_ptr<glm::vec3 []>(new glm::vec3[maxVertIndex]);
            texCoordinates = std::unique_ptr<glm::vec2 []>(new glm::vec2[maxVertIndex]);
            for (uint32_t i = 0; i < maxVertIndex; ++i) {
                texCoordinates[i] = st[i];
                multDirMatrix(normals[i], N[i],transformNormals);
            }
        }
        else {
            N = std::unique_ptr<glm::vec3 []>(new glm::vec3[numTris * 3]);
            texCoordinates = std::unique_ptr<glm::vec2 []>(new glm::vec2[numTris * 3]);
            for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
                for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) {
                    multDirMatrix(normals[k], N[l],transformNormals);
                    multDirMatrix(normals[k + j + 1], N[l + 1],transformNormals);
                    multDirMatrix(normals[k + j + 2], N[l + 2],transformNormals);
                    N[l] = glm::normalize(N[l]);
                    N[l + 1] = glm::normalize(N[l + 1]);
                    N[l + 2] = glm::normalize(N[l + 2]);
                    texCoordinates[l] = st[k];
                    texCoordinates[l + 1] = st[k + j + 1];
                    texCoordinates[l + 2] = st[k + j + 2];
                }
                k += faceIndex[i];
            }
        }

        // generate the triangle index array and set normals and st coordinates
        for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) { // for each triangle in the face
                trisIndex[l] = vertsIndex[k];
                trisIndex[l + 1] = vertsIndex[k + j + 1];
                trisIndex[l + 2] = vertsIndex[k + j + 2];
                l += 3;
            }
            k += faceIndex[i];
        }
    }
    // Test if the ray interesests this triangle mesh

    bool intersect(
            const glm::vec3 &orig,
            const glm::vec3 &dir,
            float &tNear,
            uint32_t& triIndex,
            glm::vec2 &uv,
            hit_record& rec) const override
    {
        uint32_t j = 0;
        bool isect = false;
        for (uint32_t i = 0; i < numTris; ++i) {
            const glm::vec3 &v0 = P[trisIndex[j]];
            const glm::vec3 &v1 = P[trisIndex[j + 1]];
            const glm::vec3 &v2 = P[trisIndex[j + 2]];
            float t = kInfinity, u, v;
            if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                triIndex = i;
                isect = true;
            }
            j += 3;
        }

        return isect;
    }

    void getSurfaceProperties(
            const glm::vec3 &viewDirection,
            hit_record& rec) const override
    {
        rec.material = MaterialType::DIFFUSE_AND_GLOSSY;
        uint32_t vai[3]; // vertex attr index
        if (isSingleVertAttr) {
            vai[0] = trisIndex[rec.triIndex * 3];
            vai[1] = trisIndex[rec.triIndex * 3 + 1];
            vai[2] = trisIndex[rec.triIndex * 3 + 2];
        }
        else {
            vai[0] = rec.triIndex * 3;
            vai[1] = rec.triIndex * 3 + 1;
            vai[2] = rec.triIndex * 3 + 2;
        }
        if (smoothShading) {
            // vertex normal
            const glm::vec3 &n0 = N[vai[0]];
            const glm::vec3 &n1 = N[vai[1]];
            const glm::vec3 &n2 = N[vai[2]];
            rec.normal = (1 - rec.uv.x - rec.uv.y) * n0 + rec.uv.x * n1 + rec.uv.y * n2;
        }
        else {
            // face normal
            const glm::vec3 &v0 = P[trisIndex[rec.triIndex * 3]];
            const glm::vec3 &v1 = P[trisIndex[rec.triIndex * 3 + 1]];
            const glm::vec3 &v2 = P[trisIndex[rec.triIndex * 3 + 2]];
            rec.normal = glm::cross((v1 - v0),(v2 - v0));
        }

        // doesn't need to be normalized as the N's are normalized but just for safety
        rec.normal = glm::normalize(rec.normal);

        // texture coordinates
        const glm::vec2 &st0 = texCoordinates[vai[0]];
        const glm::vec2 &st1 = texCoordinates[vai[1]];
        const glm::vec2 &st2 = texCoordinates[vai[2]];
        rec.st = (1 - rec.uv.x - rec.uv.y) * st0 + rec.uv.x * st1 + rec.uv.y * st2;
    }

    glm::vec3 checkerPattern(const glm::vec2 &st) const
    {
        float scale = 20;

        int s = (int)floor(st.x * scale);
        int t = (int)floor(st.y * scale);
        if ((s + t) % 2 == 0) {
            return glm::vec3(0.815, 0.235, 0.031);
        }
        else {
            return glm::vec3(0.937, 0.937, 0.231);
        }
    }

    glm::vec3 evalDiffuseColor(const glm::vec2 &st) const override
    {
        return checkerPattern(st);
    }


    void displayInfo() const
    {
        std::cerr << "Number of triangles in this mesh: " << numTris << std::endl;
        std::cerr << BBox[0][0] << " " << BBox[0][1] << " " << BBox[0][2] << std::endl;
        std::cerr << BBox[1][0] << " " << BBox[1][1] << " " << BBox[1][2] << std::endl;
    }
    // member variables
    uint32_t numTris;                         // number of triangles
    std::unique_ptr<glm::vec3 []> P;              // triangles vertex position
    std::unique_ptr<uint32_t []> trisIndex;   // vertex index array
    std::unique_ptr<glm::vec3 []> N;              // triangles vertex normals
    std::unique_ptr<glm::vec2 []> texCoordinates; // triangles texture coordinates
    bool smoothShading = true;                // smooth shading by default
    bool isSingleVertAttr = true;
};



#endif //PROJECT_PATCHMESH_H
