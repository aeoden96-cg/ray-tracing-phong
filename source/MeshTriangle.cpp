#include "MeshTriangle.h"

#include <memory>


MeshTriangle::MeshTriangle(
    std::unique_ptr<std::vector<glm::vec3>>& vertices_in,
    std::unique_ptr<std::vector<unsigned int>>& vertIndices_in,
    std::unique_ptr<std::vector<glm::uvec2>>& st_in,
    std::unique_ptr<std::vector<unsigned int>>& faceIndices_in,
    std::unique_ptr<std::vector<glm::vec3>>& normals_in,
    const glm::mat4 &o2w,
    bool singleVertAttr) :
        Hittable(o2w),
        isSingleVertAttr(singleVertAttr)
{
    unsigned numTris = 0;

    //k is the index of the face, maxVertIndex is the number of vertices
    unsigned int k = 0;
    unsigned int maxVertIndex = 0;

    // find out how many triangles we need to create for this mesh
    for (auto& faceInd : *faceIndices_in) {
        numTris += faceInd - 2;
        for (uint32_t j = 0; j < faceInd; ++j)
            if (vertIndices_in->at(k + j) > maxVertIndex)
                maxVertIndex = vertIndices_in->at(k + j);
        k += faceInd;
    }
    maxVertIndex += 1;

    // allocate memory to store the position of the mesh vertices
//    std::cout << "maxVertIndex: " << maxVertIndex << std::endl;
//    std::cout << "numTris: " << numTris << std::endl;
//    std::cout << vertices.size() << std::endl;

    this->vertices = std::make_unique<std::vector<glm::vec3>>(maxVertIndex);


    for (uint32_t i = 0; i < maxVertIndex; ++i) {
        multVecMatrix(vertices_in->at(i), this->vertices->at(i), objectToWorld);
    }

    // allocate memory to store triangle indices
    this->vertIndices = std::make_unique<std::vector<glm::uvec3>>(numTris);

    glm::mat4 transformNormals = glm::transpose(worldToObject);
    // [comment]
    // Sometimes we have 1 vertex attribute per vertex per face. So for example of you have 2
    // quads this would be defefined by 6 vertices but 2 * 4 vertex attribute values for
    // each vertex attribute (normal, tex. coordinates, etc.). But in some cases you may
    // want to have 1 single value per vertex. So in the quad example this would be 6 vertices
    // and 6 vertex attributes values per attribute. We need to provide both option to users.
    // [/comment]
    unsigned nfaces = faceIndices_in->size();

    if (isSingleVertAttr) {
        this->N = std::make_unique<std::vector<glm::vec3>>(maxVertIndex);
        this->st = std::make_unique<std::vector<glm::uvec2>>(maxVertIndex);

        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            this->st->at(i) = st_in->at(i);
            multDirMatrix(normals_in->at(i), N->at(i), transformNormals);
        }
    }
    else {
        this->N = std::make_unique<std::vector<glm::vec3>>(numTris);
        this->st = std::make_unique<std::vector<glm::uvec2>>(numTris);



        for (unsigned i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
            for (uint32_t j = 0; j < faceIndices_in->at(i) - 2; ++j) {
                multDirMatrix(normals_in->at(k), N->at(l), transformNormals);
                multDirMatrix(normals_in->at(k + j + 1), N->at(l+1), transformNormals);
                multDirMatrix(normals_in->at(k + j + 2), N->at(l + 2), transformNormals);
                this->N->at(l) = glm::normalize(N->at(l));
                this->N->at(l + 1) = glm::normalize(N->at(l + 1));
                this->N->at(l + 2) = glm::normalize(N->at(l + 2));
                this->st->at(l) = st_in->at(k);
                this->st->at(l + 1)  = st_in->at(k + j + 1);
                this->st->at(l + 2)  = st_in->at(k + j + 2);
            }
            k += faceIndices_in->at(i);
        }
    }

    // generate the triangle index array and set normals_in and st_in coordinates
    for (unsigned int i = 0, k = 0, l = 0; i < nfaces; ++i) { // for each  face
        for (unsigned int j = 0; j < faceIndices_in->at(i) - 2; ++j) { // for each triangle in the face
            this->vertIndices->at(l) = glm::uvec3(
                    vertIndices_in->at(k),
                    vertIndices_in->at(k + j + 1),
                    vertIndices_in->at(k + j + 2)
                    );
            l += 1;
        }
        k += faceIndices_in->at(i);
    }
}

bool MeshTriangle::intersect(
        const glm::vec3 &orig,
        const glm::vec3 &dir,
        float &tnear,
        unsigned &index,
        glm::vec2 &uv,
        hit_record& rec
) const
{
    int k = 0;
    bool intersect = false;

    for (const glm::uvec3& triangle : *this->vertIndices) {
        const glm::vec3 & v0 = this->vertices->at(triangle.x);
        const glm::vec3 & v1 = this->vertices->at(triangle.y);
        const glm::vec3 & v2 = this->vertices->at(triangle.z);
        float t, u, v;
        if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
            tnear = t;
            uv.x = u;
            uv.y = v;
            index = k;
            intersect = true;
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
    auto triangle = vertIndices->at(rec.triIndex);


    //--------------------normal--------------------
    //vai is vertex attribute index, which is the index of the vertex in the vertices array
    //triIndex is the index of the triangle in the vertIndices array

    //isSingleVertAttr is a boolean that tells us if we have 1 vertex attribute per vertex per face
    //or 1 vertex attribute per vertex

    if (isSingleVertAttr) {
        triangle = vertIndices->at(rec.triIndex);
    }
    else {
        triangle = glm::vec3(
                rec.triIndex,
                rec.triIndex,
                rec.triIndex);
    }

    if (smoothShading) {
        // vertex normal
        const glm::vec3 &n0 = N->at(triangle.x);
        const glm::vec3 &n1 = N->at(triangle.y);
        const glm::vec3 &n2 = N->at(triangle.z);
        rec.normal = (1 - rec.uv.x - rec.uv.y) * n0 + rec.uv.x * n1 + rec.uv.y * n2;
    }
    else {
        // face normal

        //now we need to get the vertices of the triangle
        const glm::vec3 & v0 = vertices->at(triangle.x);
        const glm::vec3 & v1 = vertices->at(triangle.y);
        const glm::vec3 & v2 = vertices->at(triangle.z);
        rec.normal = glm::cross((v1 - v0),(v2 - v0));
    }


    // doesn't need to be normalized as the N's are normalized but just for safety
    rec.normal = glm::normalize(rec.normal);

    //now we need to get the st coordinates of the triangle
    //st coordinates are the texture coordinates
    //we need to interpolate the st coordinates of the triangle
    //based on the uv coordinates of the intersection point

    //now we need to get the normal of the triangle using the cross product of e0 and e1
    //this is the same as the normal of the plane that the triangle lies on


    //--------------------st coordinates-------------------

    const glm::vec2 &st0 = st->at(triangle.x);
    const glm::vec2 &st1 = st->at(triangle.y);
    const glm::vec2 &st2 = st->at(triangle.z);

    //now we need to interpolate the st coordinates of the triangle
    //based on the uv coordinates of the intersection point
    //we can do this by using the barycentric coordinates of the intersection point
    //we can get the barycentric coordinates of the intersection point by using the uv coordinate
    rec.st = st0 * (1 - rec.uv.x - rec.uv.y) + st1 * rec.uv.x + st2 * rec.uv.y;

}

glm::vec3 MeshTriangle::evalDiffuseColor(const glm::vec2 &stCoords) const
{

    return checkerPattern(stCoords);
//        return gaussianNoise(st);

}


glm::vec3 MeshTriangle::gaussianNoise(const glm::vec2 &st) {
    const double mean = 0.0;
    const double stddev = 0.1;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);
    return glm::vec3(dist(generator), dist(generator), dist(generator));
}

glm::vec3 MeshTriangle::checkerPattern(const glm::vec2 &st)
{
    float scale =10;
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
