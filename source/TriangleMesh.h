#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include "utils.h"
#include "Hittable.h"


class TriangleMesh : public Hittable 
{ 
public: 
    // Build a triangle mesh from a face index array and a vertex index array
    TriangleMesh(
        const glm::mat4 &o2w,
        const uint32_t nfaces,
        std::unique_ptr<std::vector<uint32_t>> &faceIndex,
        std::unique_ptr<std::vector<uint32_t>> &vertexIndex,
        std::unique_ptr<std::vector<glm::vec3>> &vertices,
        std::unique_ptr<std::vector<glm::vec3>> &normals,
        std::unique_ptr<std::vector<glm::vec2>> &st
    ){

    }
    
    TriangleMesh( 
        const glm::mat4 &o2w, 
        const uint32_t nfaces, 
        const std::unique_ptr<uint32_t []> &faceIndex, 
        const std::unique_ptr<uint32_t []> &vertsIndex, 
        const std::unique_ptr<glm::vec3 []> &verts, 
        std::unique_ptr<glm::vec3 []> &normals, 
        std::unique_ptr<glm::vec3 []> &st, 
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
            for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) {  //for each  face 
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
        for (uint32_t i = 0, k = 0, l = 0; i < nfaces; ++i) {  //for each  face 
            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) {  //for each triangle in the face 
                trisIndex[l] = vertsIndex[k]; 
                trisIndex[l + 1] = vertsIndex[k + j + 1]; 
                trisIndex[l + 2] = vertsIndex[k + j + 2]; 
                l += 3; 
            } 
            k += faceIndex[i]; 
        } 
    } 

    void multDirMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4 &x) const 
    { 
        float a, b, c; 
 
        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0]; 
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1]; 
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2]; 
 
        dst.x = a; 
        dst.y = b; 
        dst.z = c; 
    } 

    void multVecMatrix(const glm::vec3 &src, glm::vec3 &dst,glm::mat4 &x ) const
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

    // Test if the ray interesests this triangle mesh
    bool intersect(const glm::vec3 &orig, const glm::vec3 &dir, float &tNear, uint32_t &triIndex, glm::vec2 &uv) const 
    { 
        uint32_t j = 0; 
        bool isect = false; 
        for (uint32_t i = 0; i < numTris; ++i) { 
            const glm::vec3 &v0 = P[trisIndex[j]]; 
            const glm::vec3 &v1 = P[trisIndex[j + 1]]; 
            const glm::vec3 &v2 = P[trisIndex[j + 2]]; 
            float t = kInfinity, u, v; 
            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v) && t < tNear) { 
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
        const glm::vec3 &hitPoint, 
        const glm::vec3 &viewDirection, 
        const uint32_t &triIndex, 
        const glm::vec2 &uv, 
        glm::vec3 &hitNormal, 
        glm::vec2 &hitTextureCoordinates) const 
    { 
        uint32_t vai[3];  //vertex attr index 
        if (isSingleVertAttr) { 
            vai[0] = trisIndex[triIndex * 3]; 
            vai[1] = trisIndex[triIndex * 3 + 1]; 
            vai[2] = trisIndex[triIndex * 3 + 2]; 
        } 
        else { 
            vai[0] = triIndex * 3; 
            vai[1] = triIndex * 3 + 1; 
            vai[2] = triIndex * 3 + 2; 
        } 
        if (smoothShading) { 
            // vertex normal
            const glm::vec3 &n0 = N[vai[0]]; 
            const glm::vec3 &n1 = N[vai[1]]; 
            const glm::vec3 &n2 = N[vai[2]]; 
            hitNormal = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2; 
        } 
        else { 
            // face normal
            const glm::vec3 &v0 = P[trisIndex[triIndex * 3]]; 
            const glm::vec3 &v1 = P[trisIndex[triIndex * 3 + 1]]; 
            const glm::vec3 &v2 = P[trisIndex[triIndex * 3 + 2]]; 
            hitNormal = glm::cross(v1 - v0,v2 - v0); 
        } 
 
        // doesn't need to be normalized as the N's are normalized but just for safety
        hitNormal = glm::normalize(hitNormal);
 
        // texture coordinates
        const glm::vec2 &st0 = texCoordinates[vai[0]]; 
        const glm::vec2 &st1 = texCoordinates[vai[1]]; 
        const glm::vec2 &st2 = texCoordinates[vai[2]]; 
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2; 
    } 
    void displayInfo() const 
    { 
        std::cerr << "Number of triangles in this mesh: " << numTris << std::endl; 
        //std::cerr << BBox[0] << ", " << BBox[1] << std::endl; 
    } 
    // member variables
    uint32_t numTris;                          //number of triangles 
    std::unique_ptr<glm::vec3 []> P;               //triangles vertex position 
    std::unique_ptr<uint32_t []> trisIndex;    //vertex index array 
    std::unique_ptr<glm::vec3 []> N;               //triangles vertex normals 
    std::unique_ptr<glm::vec2 []> texCoordinates;  //triangles texture coordinates 
    bool smoothShading = true;                 //smooth shading by default 
    bool isSingleVertAttr = true; 
}; 

#endif // MESHTRIANGLE_H