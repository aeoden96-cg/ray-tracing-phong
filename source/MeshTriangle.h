#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include <random>
#include <utility>
#include "utils.h"
#include "Hittable.h"

enum class MeshType {
    QUAD,
    FILE,
    POT
};

class MeshTriangle : public Hittable
{
public:
    MeshTriangle(
        std::vector<glm::vec3>& vertices_in,
        std::vector<glm::ivec3>&  vertIndices_in,
        std::vector<glm::vec2>& st_in):
        smoothShading(false),
        meshType(MeshType::QUAD)
        {

        this->vertices = std::make_unique<std::vector<glm::vec3>>(vertices_in);
        this->vertIndices = std::make_unique<std::vector<glm::ivec3>>(vertIndices_in);
        this->st = std::make_unique<std::vector<glm::vec2>>(st_in);

    }

    MeshTriangle(
            std::unique_ptr<std::vector<glm::vec3>>& vertices_in,
            std::unique_ptr<std::vector<glm::ivec3>>&  vertIndices_in,
            std::unique_ptr<std::vector<glm::vec2>>& st_in,
            std::unique_ptr<std::vector<glm::vec3>>& normals_in):
            smoothShading(false),
            meshType(MeshType::FILE)
    {

        this->vertices = std::move(vertices_in);
        this->vertIndices = std::move(vertIndices_in);
        this->st = std::move(st_in);
        this->N = std::move(normals_in);

    }


    MeshTriangle(
            std::unique_ptr<std::vector<glm::vec3>>& vertices_in,
            std::unique_ptr<std::vector<unsigned int>>& vertIndices_in,
            std::unique_ptr<std::vector<glm::vec2>>& st_in,
            std::unique_ptr<std::vector<unsigned int>>& faceIndices_in,
            std::unique_ptr<std::vector<glm::vec3>>& normals_in,
            const glm::mat4 &o2w,
            bool singleVertAttr = true);

    MeshTriangle(
            const uint32_t numFaces,
            const std::unique_ptr<std::vector<uint32_t>> &faceIndex,
            const std::unique_ptr<std::vector<uint32_t>> &vertsIndex,
            const std::unique_ptr<std::vector<glm::vec3>> &verts,
            std::unique_ptr<std::vector<glm::vec3>> &normals,
            std::unique_ptr<std::vector<glm::vec2>> &st);

    MeshTriangle(
            const uint32_t numFaces,
            const std::unique_ptr<std::vector<uint32_t>> &faceIndex,
            const std::unique_ptr<std::vector<uint32_t>> &vertsIndex,
            const std::unique_ptr<std::vector<glm::vec3>> &verts,
            std::unique_ptr<std::vector<glm::vec3>> &normals);

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

    glm::vec3 evalDiffuseColor(const glm::vec2 &st) const override;

    std::unique_ptr<std::vector<glm::vec3>> vertices;
    std::unique_ptr<std::vector<glm::ivec3>> vertIndices;
    std::unique_ptr<std::vector<glm::vec2>> st;
    std::unique_ptr<std::vector<glm::vec3>> N;              // triangles vertex normals

    MeshType meshType;
    bool smoothShading = true;       // smooth shading by default
    bool isSingleVertAttr = true;   // single vertex attribute by default


private:
    static glm::vec3 gaussianNoise(const glm::vec2 &st) ;
    static glm::vec3 checkerPattern(const glm::vec2 &st) ;

    void calcNormal(hit_record& rec) const;
    void calcST(hit_record& rec) const;

    void initVertices(UniqPointList& vertices_in,
                      unsigned maxVertIndex);
    void initVertIndices();
};




#endif // MESHTRIANGLE_H