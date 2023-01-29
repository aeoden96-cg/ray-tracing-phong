#include "MeshTriangle.h"

#include <memory>


void MeshTriangle::initVertices(
        UniqPointList& vertices_in,
        unsigned maxVertIndex
        ) {
    this->vertices = std::make_unique<PointList>(maxVertIndex);

    for (uint32_t i = 0; i < maxVertIndex; ++i) {
        multVecMatrix(vertices_in->at(i), this->vertices->at(i), objectToWorld);
    }

}

MeshTriangle::MeshTriangle(
        const uint32_t numFaces,
        const std::unique_ptr<std::vector<uint32_t>> &faceIndex,
        const std::unique_ptr<std::vector<uint32_t>> &vertsIndex,
        const std::unique_ptr<std::vector<glm::vec3>> &verts,
        std::unique_ptr<std::vector<glm::vec3>> &normals):
        meshType(MeshType::FILE)
{
    int numTris = 0;
    uint32_t k = 0, maxVertIndex = 0;
    // find out how many triangles we need to create for this mesh
    for (uint32_t i = 0; i < numFaces; ++i) {
        numTris += faceIndex->at(i) - 2;
        for (uint32_t j = 0; j < faceIndex->at(i); ++j)
            if (vertsIndex->at(k + j) > maxVertIndex)
                maxVertIndex = vertsIndex->at(k + j);
        k += faceIndex->at(i);
    }
    maxVertIndex += 1;

    // allocate memory to store the position of the mesh vertices
//        P = std::unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
    vertices = std::make_unique<std::vector<glm::vec3>>(maxVertIndex);
    for (uint32_t i = 0; i < maxVertIndex; ++i) {
        vertices->at(i) = verts->at(i);
    }

    // allocate memory to store triangle indices
//        trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
//    trisIndex = std::make_unique<std::vector<uint32_t>>(numTris * 3);
    vertIndices = std::make_unique<std::vector<glm::ivec3>>(numTris);
    uint32_t l = 0;
    // [comment]
    // Generate the triangle index array
    // Keep in mind that there is generally 1 vertex attribute for each vertex of each face.
    // So for example if you have 2 quads, you only have 6 vertices but you have 2 * 4
    // vertex attributes (that is 8 normals, 8 texture coordinates, etc.). So the easiest
    // lazziest method in our triangle mesh, is to create a new array for each supported
    // vertex attribute (st, normals, etc.) whose size is equal to the number of triangles
    // multiplied by 3, and then set the value of the vertex attribute at each vertex
    // of each triangle using the input array (normals[], st[], etc.)
    // [/comment]
//        N = std::unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
    N = std::make_unique<std::vector<glm::vec3>>(numTris * 3);
//      texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
    this->st = std::make_unique<std::vector<glm::vec2>>(numTris * 3);



    for (uint32_t i = 0, k = 0; i < numFaces; ++i) { // for each  face
        for (uint32_t j = 0; j < faceIndex->at(i) - 2; ++j) { // for each triangle in the face
//            trisIndex->at(l) = vertsIndex->at(k);
//            trisIndex->at(l + 1) = vertsIndex->at(k + j + 1);
//            trisIndex->at(l + 2) = vertsIndex->at(k + j + 2);
            vertIndices->at(l / 3) = glm::ivec3(vertsIndex->at(k), vertsIndex->at(k + j + 1), vertsIndex->at(k + j + 2));

            N->at(l) = normals->at(k);
            N->at(l + 1) = normals->at(k + j + 1);
            N->at(l + 2) = normals->at(k + j + 2);
            this->st->at(l) = st->at(k);
            this->st->at(l + 1) = st->at(k + j + 1);
            this->st->at(l + 2) = st->at(k + j + 2);
            l += 3;
        }
        k += faceIndex->at(i);
    }
    // you can use move if the input geometry is already triangulated
    //N = std::move(normals); // transfer ownership
    //sts = std::move(st); // transfer ownership


    /*std::cout << "Loaded: -------------" << std::endl;
    std::cout << "  " <<  vertices->size() << " vertices" << std::endl;
    std::cout << "  " <<  vertIndices->size() << " triangle indices" << std::endl;
    std::cout << "  " <<  N->size() << " normals" << std::endl;
    std::cout << "  " <<  this->st->size() << " texture coordinates" << std::endl;

    std:: cout <<  "  "<<"Mesh type: " << int(meshType) << std::endl;
    std:: cout <<  "  "<<"smoothShading: " << smoothShading << std::endl;
    std:: cout <<  "  "<<"isSingleVertAttr: " << isSingleVertAttr << std::endl;

    std::cout << "-------------------" << std::endl;*/


}

MeshTriangle::MeshTriangle(
        const uint32_t numFaces,
        const std::unique_ptr<std::vector<uint32_t>> &faceIndex,
        const std::unique_ptr<std::vector<uint32_t>> &vertsIndex,
        const std::unique_ptr<std::vector<glm::vec3>> &verts,
        std::unique_ptr<std::vector<glm::vec3>> &normals,
        std::unique_ptr<std::vector<glm::vec2>> &st):
        meshType(MeshType::FILE)
{
    int numTris = 0;
    uint32_t k = 0, maxVertIndex = 0;
    // find out how many triangles we need to create for this mesh
    for (uint32_t i = 0; i < numFaces; ++i) {
        numTris += faceIndex->at(i) - 2;
        for (uint32_t j = 0; j < faceIndex->at(i); ++j)
            if (vertsIndex->at(k + j) > maxVertIndex)
                maxVertIndex = vertsIndex->at(k + j);
        k += faceIndex->at(i);
    }
    maxVertIndex += 1;

    // allocate memory to store the position of the mesh vertices
//        P = std::unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
    vertices = std::make_unique<std::vector<glm::vec3>>(maxVertIndex);
    for (uint32_t i = 0; i < maxVertIndex; ++i) {
        vertices->at(i) = verts->at(i);
    }

    // allocate memory to store triangle indices
//        trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
//    trisIndex = std::make_unique<std::vector<uint32_t>>(numTris * 3);
    vertIndices = std::make_unique<std::vector<glm::ivec3>>(numTris);
    uint32_t l = 0;
    // [comment]
    // Generate the triangle index array
    // Keep in mind that there is generally 1 vertex attribute for each vertex of each face.
    // So for example if you have 2 quads, you only have 6 vertices but you have 2 * 4
    // vertex attributes (that is 8 normals, 8 texture coordinates, etc.). So the easiest
    // lazziest method in our triangle mesh, is to create a new array for each supported
    // vertex attribute (st, normals, etc.) whose size is equal to the number of triangles
    // multiplied by 3, and then set the value of the vertex attribute at each vertex
    // of each triangle using the input array (normals[], st[], etc.)
    // [/comment]
//        N = std::unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
    N = std::make_unique<std::vector<glm::vec3>>(numTris * 3);
//      texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
    this->st = std::make_unique<std::vector<glm::vec2>>(numTris * 3);



    for (uint32_t i = 0, k = 0; i < numFaces; ++i) { // for each  face
        for (uint32_t j = 0; j < faceIndex->at(i) - 2; ++j) { // for each triangle in the face
//            trisIndex->at(l) = vertsIndex->at(k);
//            trisIndex->at(l + 1) = vertsIndex->at(k + j + 1);
//            trisIndex->at(l + 2) = vertsIndex->at(k + j + 2);
            vertIndices->at(l / 3) = glm::ivec3(vertsIndex->at(k), vertsIndex->at(k + j + 1), vertsIndex->at(k + j + 2));

            N->at(l) = normals->at(k);
            N->at(l + 1) = normals->at(k + j + 1);
            N->at(l + 2) = normals->at(k + j + 2);
            this->st->at(l) = st->at(k);
            this->st->at(l + 1) = st->at(k + j + 1);
            this->st->at(l + 2) = st->at(k + j + 2);
            l += 3;
        }
        k += faceIndex->at(i);
    }
    // you can use move if the input geometry is already triangulated
    //N = std::move(normals); // transfer ownership
    //sts = std::move(st); // transfer ownership


    /*std::cout << "Loaded: -------------" << std::endl;
    std::cout << "  " <<  vertices->size() << " vertices" << std::endl;
    std::cout << "  " <<  vertIndices->size() << " triangle indices" << std::endl;
    std::cout << "  " <<  N->size() << " normals" << std::endl;
    std::cout << "  " <<  this->st->size() << " texture coordinates" << std::endl;

    std:: cout <<  "  "<<"Mesh type: " << int(meshType) << std::endl;
    std:: cout <<  "  "<<"smoothShading: " << smoothShading << std::endl;
    std:: cout <<  "  "<<"isSingleVertAttr: " << isSingleVertAttr << std::endl;

    std::cout << "-------------------" << std::endl;*/

}

MeshTriangle::MeshTriangle(
        UniqPointList& vertices_in,
        UniqIndexList& vertIndices_in,
        UniqTexList& st_in,
        UniqIndexList& faceIndices_in,
        UniqPointList& normals_in,
        const glm::mat4 &o2w,
        bool singleVertAttr) :
        Hittable(o2w),
        isSingleVertAttr(singleVertAttr),
        meshType(MeshType::POT)
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

    this->initVertices(vertices_in, maxVertIndex);


    // allocate memory to store triangle indices
    this->vertIndices = std::make_unique<std::vector<glm::ivec3>>(numTris);

    glm::mat4 transformNormals = glm::transpose(worldToObject);
    // [comment]
    // Sometimes we have 1 vertex attribute per vertex per face. So for example of you have 2
    // quads this would be defefined by 6 vertices but 2 * 4 vertex attribute values for
    // each vertex attribute (normal, tex. coordinates, etc.). But in some cases you may
    // want to have 1 single value per vertex. So in the quad example this would be 6 vertices
    // and 6 vertex attributes values per attribute. We need to provide both option to users.
    // [/comment]
    unsigned nfaces = faceIndices_in->size();


    //isSingleVertAttr is true by default
    //if isSingleVertAttr is true, then we have 1 vertex attribute per vertex per face
    //if isSingleVertAttr is false, then we have 1 vertex attribute per vertex
    if (isSingleVertAttr) {
        this->N = std::make_unique<std::vector<glm::vec3>>(maxVertIndex);
        this->st = std::make_unique<std::vector<glm::vec2>>(maxVertIndex);

        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            this->st->at(i) = st_in->at(i);
            multDirMatrix(normals_in->at(i), N->at(i), transformNormals);
        }
    }
    else {
        this->N = std::make_unique<std::vector<glm::vec3>>(numTris);
        this->st = std::make_unique<std::vector<glm::vec2>>(numTris);



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
            this->vertIndices->at(l) = glm::vec3(
                    vertIndices_in->at(k),
                    vertIndices_in->at(k + j + 1),
                    vertIndices_in->at(k + j + 2)
            );
            l += 1;
        }
        k += faceIndices_in->at(i);
    }

    /*std::cout << "Loaded: -------------" << std::endl;
    std::cout << "  " <<  vertices->size() << " vertices" << std::endl;
    std::cout << "  " <<  vertIndices->size() << " triangle indices" << std::endl;
    std::cout << "  " <<  N->size() << " normals" << std::endl;
    std::cout << "  " <<  this->st->size() << " texture coordinates" << std::endl;

    std:: cout <<  "  "<<"Mesh type: " << int(meshType) << std::endl;
    std:: cout <<  "  "<<"smoothShading: " << smoothShading << std::endl;
    std:: cout <<  "  "<<"isSingleVertAttr: " << isSingleVertAttr << std::endl;

    std::cout << "-------------------" << std::endl;*/
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
    bool intersect = false;
    uint32_t j = 0;
    for (uint32_t i = 0; i < vertIndices->size(); ++i) {
        glm::vec3 v0 = vertices->at(vertIndices->at(i).x);
        glm::vec3 v1 = vertices->at(vertIndices->at(i).y);
        glm::vec3 v2 = vertices->at(vertIndices->at(i).z);
        float t, u, v;
        if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
            tnear = t;
            uv.x = u;
            uv.y = v;
            index = i;
            intersect = true;
        }
        j += 3;
    }
    return intersect;
}


void MeshTriangle::calcNormal(hit_record& rec) const{

    glm::ivec3 triangle;

    if (isSingleVertAttr) {
        triangle = vertIndices->at(rec.triIndex);
    }
    else {
        triangle = glm::vec3(rec.triIndex,rec.triIndex,rec.triIndex);
    }

    if (smoothShading) {
        const glm::vec3 &n0 = N->at(triangle.x);
        const glm::vec3 &n1 = N->at(triangle.y);
        const glm::vec3 &n2 = N->at(triangle.z);
        rec.normal = (1 - rec.uv.x - rec.uv.y) * n0 + rec.uv.x * n1 + rec.uv.y * n2;
    }
    else {

        const glm::vec3 & v0 = vertices->at(triangle.x);
        const glm::vec3 & v1 = vertices->at(triangle.y);
        const glm::vec3 & v2 = vertices->at(triangle.z);
        rec.normal = glm::cross((v1 - v0),(v2 - v0));
    }

    rec.normal = glm::normalize(rec.normal);

}

void MeshTriangle::calcST(hit_record& rec) const{

    if (meshType == MeshType::QUAD) {

        auto triangle = vertIndices->at(rec.triIndex);

        const glm::vec2 &st0 = st->at(triangle.x);
        const glm::vec2 &st1 = st->at(triangle.y);
        const glm::vec2 &st2 = st->at(triangle.z);

        rec.st = st0 * (1 - rec.uv.x - rec.uv.y) + st1 * rec.uv.x + st2 * rec.uv.y;
    }
    else if( meshType == MeshType::FILE)
    {
        auto triangle = vertIndices->at(rec.triIndex);

        const glm::vec2 &st01 = st->at(triangle.x);
        const glm::vec2 &st11 = st->at(triangle.y);
        const glm::vec2 &st21 = st->at(triangle.z);


        rec.st = (1 - rec.uv.x - rec.uv.y) * st01 + rec.uv.x * st11 + rec.uv.y * st21;
//        rec.st = {triangle.x,triangle.y};
    }


}



void MeshTriangle::getSurfaceProperties(
        const glm::vec3 &I,
        hit_record& rec) const
{

    rec.material = MaterialType::DIFFUSE_AND_GLOSSY;
    this->calcNormal(rec);
    this->calcST(rec);

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
    float scale =5;

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

//    float pattern = (fmodf(st.x * scale, 1) > 0.5f) ^ (fmodf(st.y * scale, 1) > 0.5);
//    return mix(
//        glm::vec3(0.815, 0.235, 0.031),
//        glm::vec3(0.937, 0.937, 0.231),
//        pattern);
}

