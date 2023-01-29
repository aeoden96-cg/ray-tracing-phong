//[header]
// A simple program to demonstrate how to implement Whitted-style ray-tracing
//[/header]
//[compile]
// Download the whitted.cpp file to a folder.
// Open a shell/terminal, and run the following command where the files is saved:
//
// c++ -o whitted whitted.cpp -std=c++11 -O3
//
// Run with: ./whitted. Open the file ./out.png in Photoshop or any program
// reading PPM files.
//[/compile]
//[ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//[/ignore]


#include "utils.h"
#include "Hittable.h"
#include "Sphere.h"
#include "MeshTriangle.h"
#include "Light.h"
#include "teapotdata.h"
#include <boost/range/combine.hpp>
#include <memory>

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>

#include "yaml-cpp/yaml.h"





void createPolyTeapot(const glm::mat4& o2w, std::vector<std::unique_ptr<Hittable>> &objects)
{
    uint32_t divs = 8;

    std::unique_ptr<PointList> verts = std::make_unique<PointList>((divs + 1) * (divs + 1));
    std::unique_ptr<IndexList> faceIndices = std::make_unique<IndexList>(divs * divs);
    std::unique_ptr<IndexList> vertIndices = std::make_unique<IndexList>(divs * divs * 4);
    std::unique_ptr<PointList> normals = std::make_unique<PointList>((divs + 1) * (divs + 1));
    std::unique_ptr<TexSTList> st = std::make_unique<TexSTList>((divs + 1) * (divs + 1));

    // face connectivity - all patches are subdivided the same way so there
    // share the same topology and uvs
    for (unsigned j = 0, k = 0; j < divs; ++j) {
        for (unsigned i = 0; i < divs; ++i, ++k) {
            faceIndices->at(k) = 4;
            vertIndices->at(k * 4) = (divs + 1) * j + i;
            vertIndices->at(k * 4 + 1)  = (divs + 1) * j + i + 1;
            vertIndices->at(k * 4 + 2)  = (divs + 1) * (j + 1) + i + 1;
            vertIndices->at(k * 4 + 3)  = (divs + 1) * (j + 1) + i;
        }
    }

    std::vector<glm::vec3> controlPoints(16);
    for (auto & teapotPatch : teapotPatches) {  //kTeapotNumPatches
        // set the control points for the current patch
        for (auto tup : boost::combine(controlPoints, teapotPatch)) {
            //boost:combine is a zip function, it combines two vectors into a tuple, which is then unpacked
            glm::vec3 &cp = tup.get<0>();
            unsigned &idx = tup.get<1>();
            cp = teapotVertices[idx - 1];
        }


        // generate grid
        for (unsigned j = 0, k = 0; j <= divs; ++j) {
            float v = (float)j / (float)divs;
            for (unsigned i = 0; i <= divs; ++i, ++k) {
                float u = (float)i / (float)divs;
                verts->at(k) = evalBezierPatch(controlPoints, u, v);
                glm::vec3  dU = dUBezier(controlPoints, u, v);
                glm::vec3  dV = dVBezier(controlPoints, u, v);
                normals->at(k) = glm::normalize(glm::cross(dU, dV));
                st->at(k).x = u;
                st->at(k).y = v;
            }
        }


        std::unique_ptr<MeshTriangle> meshTriangle =
                std::make_unique<MeshTriangle>(
                        verts,
                        vertIndices,
                        st,
                        faceIndices,
                        normals,
                        o2w,
                        divs * divs);

        meshTriangle->materialType = MaterialType::COW;


        objects.push_back(std::move(meshTriangle));
    }
}

void loadTinyOBJFromFile(std::string filename,std::vector<std::unique_ptr<Hittable>> &objects, glm::mat4 o2w){
    std::string inputfile = filename;
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./"; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();

    int shapeIndex = 0;
    auto myShape = shapes[shapeIndex];
    auto numOfFaces = myShape.mesh.num_face_vertices.size();

    int numberOfPoints = 0;
    for (size_t f = 0; f < numOfFaces; f++) {
        int fv = myShape.mesh.num_face_vertices[f];
        for (size_t v = 0; v < fv; v++) {
            tinyobj::index_t idx = myShape.mesh.indices[3 * f + v];
            if (idx.vertex_index > numberOfPoints) {
                numberOfPoints = idx.vertex_index;
            }
        }
    }

    numberOfPoints++;

    std::cout << "numberOfPoints: " << numberOfPoints << std::endl;


    std::unique_ptr<std::vector<uint32_t>> faceIndex(new std::vector<uint32_t>(numOfFaces));
    std::unique_ptr<std::vector<glm::ivec3>> vertIndices(new std::vector<glm::ivec3>(numOfFaces));
    std::unique_ptr<std::vector<glm::vec3>> vertices(new std::vector<glm::vec3>(numberOfPoints));
    std::unique_ptr<std::vector<glm::vec3>> normals(new std::vector<glm::vec3>(numberOfPoints));
    std::unique_ptr<std::vector<glm::vec2>> st(new std::vector<glm::vec2>(numberOfPoints));


    int vertsArraySize = 0;
    size_t index_offset = 0;

    // Loop over shapes
    //for (size_t s = 0; s < shapes.size(); s++) {
    for (size_t faceInd = 0; faceInd < myShape.mesh.num_face_vertices.size(); faceInd++) {
        size_t numOfFacePoints = size_t(myShape.mesh.num_face_vertices[faceInd]);

        //std::cout << "numOfFacePoints: " << numOfFacePoints << std::endl;
        faceIndex->at(faceInd) = numOfFacePoints;

        std::cout << "Face " << faceInd << " has " << numOfFacePoints << " vertices" << std::endl;

        glm::ivec3 oneVertIndices;

        // Loop over vertices in the face.
        for (size_t vertInd = 0; vertInd < numOfFacePoints; vertInd++) {
            tinyobj::index_t idx = myShape.mesh.indices[index_offset + vertInd];

            oneVertIndices[vertInd] = idx.vertex_index;

            std::cout << "  Point " << vertInd << " has vertex index " << idx.vertex_index << std::endl;

            // access to vertex index
            std::cout << "      vertex_index: " << idx.vertex_index << std::endl;
            std::cout << "      normal_index: " << idx.normal_index << std::endl;
            std::cout << "      texcoord_index: " << idx.texcoord_index << std::endl;



            // access to vertex
            tinyobj::real_t vx = attrib.vertices[3*size_t(idx.vertex_index)+0];
            tinyobj::real_t vy = attrib.vertices[3*size_t(idx.vertex_index)+1];
            tinyobj::real_t vz = attrib.vertices[3*size_t(idx.vertex_index)+2];

            glm::vec3 vert = glm::vec3(vx, vy, vz);

            // transform the vertex
            vert = glm::vec3(o2w * glm::vec4(vert, 1.0f));

            std::cout << "      vx: " << vx << " vy: " << vy << " vz: " << vz << std::endl;

            vertices->at(idx.vertex_index) = vert;


            // Check if `normal_index` is zero or positive. negative = no normal data
            if (idx.normal_index >= 0) {
                tinyobj::real_t nx = attrib.normals[3*size_t(idx.normal_index)+0];
                tinyobj::real_t ny = attrib.normals[3*size_t(idx.normal_index)+1];
                tinyobj::real_t nz = attrib.normals[3*size_t(idx.normal_index)+2];

                std::cout << "      nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
                normals->at(idx.vertex_index) = glm::vec3(nx, ny, nz);
            }

            // Check if `texcoord_index` is zero or positive. negative = no texcoord data
            if (idx.texcoord_index >= 0) {
                tinyobj::real_t tx = attrib.texcoords[2*size_t(idx.texcoord_index)+0];
                tinyobj::real_t ty = attrib.texcoords[2*size_t(idx.texcoord_index)+1];

                std::cout << "      tx: " << tx << " ty: " << ty << std::endl;

                st->at(idx.vertex_index) = glm::vec2(tx, ty);
            }

            else {
                std::cout << "no texcoord data" << std::endl;
                exit(-1);
            }

            // Optional: vertex colors
            // tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
            // tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
            // tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
        }
        vertIndices->at(faceInd) = oneVertIndices;

//        std::cout << oneVertIndices.x << " " << oneVertIndices.y << " " << oneVertIndices.z << std::endl;
        index_offset += numOfFacePoints;

        // per-face material
        // myShape.mesh.material_ids[faceInd];

    }

    // exit app

    //std::unique_ptr<MeshTriangle> meshTriangle =
    //        std::make_unique<MeshTriangle>(
     //               numFaces, faceIndex, vertsIndex, verts, normals, st);

    //meshTriangle->materialType = MaterialType::COW;

    //objects.push_back(std::move(meshTriangle));

    std::cout << "num of points " << vertices->size() << std::endl;
    std::cout << "num of normals " << normals->size() << std::endl;
    std::cout << "num of st " << st->size() << std::endl;

    std::unique_ptr<MeshTriangle> meshTriangle =
            std::make_unique<MeshTriangle>(
                        vertices, vertIndices, st, normals);

    meshTriangle->materialType = MaterialType::COW;

    objects.push_back(std::move(meshTriangle));
}

void loadPolyMeshFromFile(const glm::mat4& o2w, std::vector<std::unique_ptr<Hittable>> &objects)
{
    std::ifstream ifs;

    ifs.open("cow.geo");
    if (ifs.fail()) throw;
    std::stringstream ss;
    ss << ifs.rdbuf();
    uint32_t numFaces;
    ss >> numFaces;
    std::unique_ptr<std::vector<uint32_t>> faceIndex(new std::vector<uint32_t>(numFaces));
    uint32_t vertsIndexArraySize = 0;
    // reading face index array
    for (uint32_t i = 0; i < numFaces; ++i) {
        ss >> faceIndex->at(i);
        vertsIndexArraySize += faceIndex->at(i);
    }
//    std::cout << "numFaces: " << numFaces << std::endl;

    std::unique_ptr<std::vector<uint32_t>> vertsIndex(new std::vector<uint32_t>(vertsIndexArraySize));
    uint32_t vertsArraySize = 0;
    // reading vertex index array
    for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
        ss >> vertsIndex->at(i);
        if (vertsIndex->at(i) > vertsArraySize) vertsArraySize = vertsIndex->at(i);
    }
//    std::cout << "vertsIndexArraySize: " << vertsIndexArraySize << std::endl;

    vertsArraySize += 1;
    // reading vertices
    std::unique_ptr<std::vector<glm::vec3>> verts(new std::vector<glm::vec3>(vertsArraySize));
    for (uint32_t i = 0; i < vertsArraySize; ++i) {
        ss >> verts->at(i).x >> verts->at(i).y >> verts->at(i).z;
        // o2w is the object to world transform
        glm::vec3 res;
        multVecMatrix(verts->at(i), res, o2w);
        verts->at(i) = res;
    }
//    std::cout << "vertsArraySize: " << vertsArraySize << std::endl;

    // reading normals
    std::unique_ptr<std::vector<glm::vec3>> normals(new std::vector<glm::vec3>(vertsIndexArraySize));
    for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
        ss >> normals->at(i).x >> normals->at(i).y >> normals->at(i).z;
    }
//    std::cout << "normalsArraySize: " << vertsArraySize << std::endl;

    // reading st coordinates
    std::unique_ptr<std::vector<glm::vec2>> st(new std::vector<glm::vec2>(vertsIndexArraySize));
    for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
        ss >> st->at(i).x >> st->at(i).y;
    }
//    std::cout << "stArraySize: " << vertsArraySize << std::endl;

    ifs.close();


    std::unique_ptr<MeshTriangle> meshTriangle =
            std::make_unique<MeshTriangle>(
                    numFaces, faceIndex, vertsIndex, verts, normals, st);

    meshTriangle->materialType = MaterialType::COW;

    objects.push_back(std::move(meshTriangle));


}



HittableList loadSceneFromFile(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename + ".yaml");

    HittableList objects;

    std::map<std::string, MaterialType> materials;

    materials["DIFFUSE_AND_GLOSSY"] = DIFFUSE_AND_GLOSSY;
    materials["REFLECTION_AND_REFRACTION"] = REFLECTION_AND_REFRACTION;
    materials["REFLECTION"] = REFLECTION;
    materials["METAL"] = METAL;



    for (auto object : config["objects"]) {
        auto type = object["type"].as<std::string>();
        auto material_name = object["material"].as<std::string>();
        auto material = materials[material_name];

        if(type == "sphere") {
            auto center_v = object["center"].as<std::vector<float>>();
            auto radius = object["radius"].as<float>();
            auto color = object["color"].as<std::vector<float>>();
            Sphere s(glm::vec3(center_v[0], center_v[1], center_v[2]), radius);
            s.materialType = material;
            s.diffuseColor = glm::vec3(color[0], color[1], color[2]);
            objects.push_back(std::make_unique<Sphere>(s));
        } else if(type == "plane") {
            auto up_left_v = object["up_left"].as<std::vector<float>>();
            auto up_right_v = object["up_right"].as<std::vector<float>>();
            auto down_left_v = object["down_left"].as<std::vector<float>>();
            auto down_right_v = object["down_right"].as<std::vector<float>>();
//            material_name = object["material"].as<std::string>();
            //auto material = materials[material_name];

            std::vector<glm::vec3> vertices;
            vertices.emplace_back(up_left_v[0], up_left_v[1], up_left_v[2]);
            vertices.emplace_back(up_right_v[0], up_right_v[1], up_right_v[2]);
            vertices.emplace_back(glm::vec3(down_right_v[0], down_right_v[1], down_right_v[2]));
            vertices.emplace_back(glm::vec3(down_left_v[0], down_left_v[1], down_left_v[2]));

            std::vector<glm::ivec3> vertIndices;
            vertIndices.emplace_back(glm::vec3(0, 1, 3));
            vertIndices.emplace_back(glm::vec3(1, 2, 3));

            std::vector<glm::vec2> uvIndices;
            uvIndices.emplace_back(glm::vec2(0, 0));
            uvIndices.emplace_back(glm::vec2(1, 0));
            uvIndices.emplace_back(glm::vec2(1, 1));
            uvIndices.emplace_back(glm::vec2(0, 1));


            objects.emplace_back(std::make_unique<MeshTriangle>(
                    vertices, vertIndices, uvIndices));

        } else if(type == "bezier") {
            const auto numOfObjects = objects.size();

            glm::mat4 o2w = glm::translate(glm::mat4(1.0f), glm::vec3(3 ,-5, -7));
            //to rotate it, we need to rotate it by 90 degrees around the x-axis
            o2w = glm::rotate(o2w, glm::radians(-90.0f), glm::vec3(1, 0, 0));
            o2w = glm::rotate(o2w, glm::radians(-70.0f), glm::vec3(0, 0, 1));
            //to scale it, we need to scale it by 0.5
            o2w = glm::scale(o2w, glm::vec3(0.7, 0.7, 0.7));


            createPolyTeapot(o2w, objects);
            std::cout << "Created poly teapot, which has " << objects.size() - numOfObjects << " objects" << std::endl;


        } else if(type == "mesh") {

            const auto numOfObjects = objects.size();

            glm::mat4  o2w = glm::translate(glm::mat4(1.0f), glm::vec3(-3, -5, -8));
            o2w = glm::scale(o2w, glm::vec3(0.5, 0.5, 0.5));

            loadPolyMeshFromFile(o2w,objects);

            std::cout << "Created poly mesh, which has " << objects.size() - numOfObjects << " objects" << std::endl;

        } else {
            throw std::runtime_error("Unknown object type");
        }
    }

    glm::mat4 o2w = glm::translate(glm::mat4(1.0f), glm::vec3(0 ,-0, -5));
    //to rotate it, we need to rotate it by 90 degrees around the x-axis
    o2w = glm::rotate(o2w, glm::radians(0.0f), glm::vec3(1, 0, 0));
    o2w = glm::rotate(o2w, glm::radians(140.0f), glm::vec3(0, 1, 0));
    //to scale it, we need to scale it by 0.5
    o2w = glm::scale(o2w, glm::vec3(0.1, 0.1, 0.1));

    loadTinyOBJFromFile("dino.obj", objects, o2w);


    return objects;
}


// Returns true if the ray intersects a Hittable, false otherwise.
// \param orig is the ray origin
// \param dir is the ray direction
// \param objects is the list of objects the scene contains
// \param[out] tNear contains the distance to the closest intersected Hittable.
// \param[out] index stores the index of the intersect triangle if the interested Hittable is a mesh.
// \param[out] uv stores the u and v barycentric coordinates of the intersected point
// \param[out] *hitObject stores the pointer to the intersected Hittable (used to retrieve material information, etc.)
// \param isShadowRay is it a shadow ray. We can return from the function sooner as soon as we have found a hit.
bool trace(
    const glm::vec3 &orig, const glm::vec3 &dir,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    hit_record &rec)
{
    for (const auto &object : objects) {
        float tNearK = kInfinity;
        uint32_t indexK;
        glm::vec2 uvK;
        if (object->intersect(
                orig,
                dir,
                tNearK,
                indexK,
                uvK,
                rec) && tNearK < rec.t) {
            rec.t = tNearK;
            rec.triIndex = indexK;
            rec.uv = uvK;
            rec.p = orig + dir * rec.t;
            rec.t = tNearK;
            rec.uv = uvK;
            rec.object = object.get();
        }
    }

    return (rec.object != nullptr);
}


glm::vec3 castRay(
    const glm::vec3 &orig, const glm::vec3 &dir,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights,
    const Options &options,
    uint32_t depth);


glm::vec3 getReflectionColor (
    hit_record& rec, 
    const glm::vec3& dir, 
    const  Options& options,
    unsigned depth,
    const std::vector<std::unique_ptr<Hittable>>& objects,
    const std::vector<std::unique_ptr<Light>>& lights);

glm::vec3 randomInUnitSphere()
{
    glm::vec3 p;
    do {
        p = 2.0f * glm::vec3(drand48(), drand48(), drand48()) - glm::vec3(1, 1, 1);
    } while (glm::dot(p, p) >= 1.0);
    return p;
}

glm::vec3 getRefractionColor (
    hit_record& rec, 
    const glm::vec3& dir, 
    const Options& options,
    unsigned depth,
    const std::vector<std::unique_ptr<Hittable>>& objects,
    const std::vector<std::unique_ptr<Light>>& lights)
{
    glm::vec3 refractionDirection = normalize(refract(dir, rec.normal, rec.object->ior));
                
    glm::vec3 refractionRayOrig = (glm::dot(refractionDirection, rec.normal) < 0) ?
        rec.p - rec.normal * options.bias :
        rec.p + rec.normal * options.bias;

    return castRay(
        refractionRayOrig, 
        refractionDirection, 
        objects, 
        lights, 
        options, 
        depth + 1);
    
}

// Implementation of the Whitted-style light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected Hittable is either reflective or reflective and refractive,
// then we compute the reflection/refraction direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refraction depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is diffuse/glossy we use the Phong illumination model to compute the color
// at the intersection point.
glm::vec3 castRay(
    const glm::vec3 &orig, const glm::vec3 &dir,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights,
    const Options &options,
    unsigned depth)
{
    if (depth > options.maxDepth) {
        return options.backgroundColor;
    }

    glm::vec3 hitColor = options.backgroundColor;
    hit_record rec;
    
    if (trace(orig, dir, objects, rec)) {

        rec.object->getSurfaceProperties(dir, rec);
        
        glm::vec2 uv = rec.uv;
        glm::vec3 N = rec.normal;
        glm::vec2 st = rec.st;
        glm::vec3 hitPoint = rec.p;

        switch (rec.object->materialType) {
            case REFLECTION_AND_REFRACTION:
            {
                float kr = fresnel(dir, N, rec.object->ior);

                glm::vec3 reflectionColor = getReflectionColor(
                    rec, dir, options, depth, objects, lights);

                glm::vec3 refractionColor = getRefractionColor(
                    rec, dir, options, depth, objects, lights);


                hitColor = reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            case REFLECTION:
            {
                float kr = fresnel(dir, N, rec.object->ior);

                glm::vec3 reflectionColor = getReflectionColor(
                    rec, dir, options, depth, objects, lights);

                hitColor = reflectionColor * kr;
                break;
            }
            case METAL:
            {
                float roughness = 0.1f;
                float kr = fresnel(dir, N, rec.object->ior);

                glm::vec3 reflectionColor = getReflectionColor(
                    rec,
                    dir + roughness * randomInUnitSphere(),
                    options,
                    depth,
                    objects,
                    lights);

                hitColor = reflectionColor * kr;
                break;
            }
            case COW:
            {
                float NdotView = std::max(0.f, glm::dot(N, -dir));
                const int M = 10;
                float checker = (fmod(st.x * M, 1.0f) > 0.5f) ^ (fmod(st.y * M, 1.0f) < 0.5f);
                float c = 0.3 * (1 - checker) + 0.7 * checker;

                hitColor = { c * NdotView,c * NdotView,c * NdotView }; //Vec3f(uv.x, uv.y, 0);
                break;
            }
            default:
            {
//                hitColor = rec.object->evalDiffuseColor(st);
//                break;

                // We use the Phong illumination model int the default case. The phong model
                // is composed of a diffuse and a specular reflection component.

                glm::vec3 lightAmt = glm::vec3(0,0,0), specularColor =  glm::vec3(0,0,0);
                glm::vec3 shadowPointOrig = (glm::dot(dir, N) < 0) ?
                    hitPoint + N * options.bias :
                    hitPoint - N * options.bias;

                // Loop over all lights in the scene and sum their contribution up
                // We also apply the lambert cosine law here though we haven't explained yet what this means.
                for (const auto &light : lights) {
                    glm::vec3 lightDir = light->position - hitPoint;
                    lightDir = normalize(lightDir);
                    float cosine = std::max(0.f, glm::dot(lightDir, N));
                    // cosine is the cosine of the angle between the light direction and the normal
                    // at the intersection point. If cosine is 0, then the light is perpendicular to the
                    // surface and thus does not contribute to the illumination of the surface.
                    // If cosine is 1, then the light is parallel to the surface and thus contributes
                    // the maximum amount of illumination to the surface.
                    
                    //Hittable *shadowHitObject = nullptr;
                    hit_record shadowRec;
                    // is the point in shadow, and is the nearest occluding Hittable closer to the Hittable than the light itself?
                    bool inShadow = trace(shadowPointOrig, lightDir, objects, shadowRec) &&
                        shadowRec.t * shadowRec.t < glm::dot(lightDir, lightDir);

                    if(!inShadow) {
                       lightAmt += light->intensity * cosine;
                    }

                    glm::vec3 reflectionDirection = reflect(-lightDir, N);
                    specularColor += powf(std::max(0.f, -glm::dot(reflectionDirection, dir)), rec.object->specularExponent) * light->intensity;
                }
                hitColor = lightAmt * rec.object->evalDiffuseColor(st) * rec.object->Kd + specularColor * rec.object->Ks;
                break;
            }
        }
    }

    return hitColor;
}

glm::vec3 getReflectionColor(hit_record &rec, const glm::vec3 &dir, const Options &options, unsigned int depth,
                             const std::vector<std::unique_ptr<Hittable>> &objects,
                             const std::vector<std::unique_ptr<Light>> &lights) {
    glm::vec3 reflectionDirection = normalize(reflect(dir, rec.normal));
    glm::vec3 reflectionRayOrig = (glm::dot(reflectionDirection, rec.normal) < 0) ?
                                  rec.p - rec.normal * options.bias :
                                  rec.p + rec.normal * options.bias;

    return castRay(
            reflectionRayOrig,
            reflectionDirection,
            objects,
            lights,
            options,
            depth + 1);
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void render(
        Options &options,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights)
{
    glm::vec3 orig(0);
    multVecMatrix(glm::vec3(0), orig,options.cameraToWorld);
    auto framebuffer = std::vector<glm::vec3>(options.width * options.height);
    auto scale = (float)tan(glm::radians(options.fov * 0.5));
    float imageAspectRatio = (float)options.width / (float)options.height;

    int index = 0;
    for (unsigned j = 0; j < options.height; ++j) {
        for (unsigned i = 0; i < options.width; ++i) {
            // generate primary ray direction
            float x = (2.0f * ((float)i + 0.5f) / (float)options.width - 1.0f) * imageAspectRatio * scale;
            float y = (1.0f - 2.0f * ((float)j + 0.5f) / (float)options.height) * scale;
            // glm::vec3 dir = normalize(glm::vec3(x, y, -1));
            glm::vec3 dir;
            multDirMatrix(glm::vec3(x, y, -1), dir,options.cameraToWorld);
            dir = glm::normalize(dir);
            framebuffer[index++] = castRay(orig, dir, objects, lights, options, 0);
        }
        std::cout << "\rRendering (" << 1 << " spp) " << 100.f * (float)j / (float)(options.height - 1) << "%" << std::flush;
    }

    // save framebuffer to file
    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (uint32_t i = 0; i < options.height * options.width; ++i) {
        char r = (char)(255 * clamp_(0, 1, framebuffer[i].x));
        char g = (char)(255 * clamp_(0, 1, framebuffer[i].y));
        char b = (char)(255 * clamp_(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }

    ofs.close();
}

// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image width and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().


//void createCurveGeometry(std::vector<std::unique_ptr<Hittable>> &objects)
//{
//    unsigned ndivs = 16;
//    unsigned ncurves = 1 + (curveNumPts - 4) / 3;
//    std::vector<glm::vec3> pts(4);
//
//    auto P = std::make_unique<std::vector<glm::vec3>>((ndivs + 1) * ndivs * ncurves + 1);
//    auto N = std::make_unique<std::vector<glm::vec3>>((ndivs + 1) * ndivs * ncurves + 1);
//    auto st = std::make_unique<std::vector<glm::uvec2>>((ndivs + 1) * ndivs * ncurves + 1);
//
//
//
//    for (uint32_t i = 0; i < ncurves; ++i) {
//        for (uint32_t j = 0; j < ndivs; ++j) {
//            pts[0] = curveData[i * 3];
//            pts[1] = curveData[i * 3 + 1];
//            pts[2] = curveData[i * 3 + 2];
//            pts[3] = curveData[i * 3 + 3];
//            float s = j / (float)ndivs;
//            glm::vec3 pt = evalBezierCurve(pts, s);
//            glm::vec3 tangent = glm::normalize(derivBezier(pts, s));
//            bool swap = false;
//
//            uint8_t maxAxis;
//            if (std::abs(tangent.x) > std::abs(tangent.y))
//                if (std::abs(tangent.x) > std::abs(tangent.z))
//                    maxAxis = 0;
//                else
//                    maxAxis = 2;
//            else if (std::abs(tangent.y) > std::abs(tangent.z))
//                maxAxis = 1;
//            else
//                maxAxis = 2;
//
//            glm::vec3 up, forward, right;
//
//            switch (maxAxis) {
//                case 0:
//                case 1:
//                    up = tangent;
//                    forward = glm::vec3(0, 0, 1);
//                    right = glm::cross(up,forward);
//                    forward = glm::cross(right,up);
//                    break;
//                case 2:
//                    up = tangent;
//                    right = glm::vec3(0, 0, 1);
//                    forward = glm::cross(right,up);
//                    right = glm::cross(up,forward);
//                    break;
//                default:
//                    break;
//            };
//
//            float sNormalized = (i * ndivs + j) / float(ndivs * ncurves);
//            float rad = 0.1 * (1 - sNormalized);
//            for (uint32_t k = 0; k <= ndivs; ++k) {
//                float t = k / (float)ndivs;
//                float theta = t * 2 * M_PI;
//                glm::vec3 pc(cos(theta) * rad, 0, sin(theta) * rad);
//                float x = pc.x * right.x + pc.y * up.x + pc.z * forward.x;
//                float y = pc.x * right.y + pc.y * up.y + pc.z * forward.y;
//                float z = pc.x * right.z + pc.y * up.z + pc.z * forward.z;
//                P[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = glm::vec3(pt.x + x, pt.y + y, pt.z + z);
//                N[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = glm::normalize(glm::vec3(x, y, z));
//                st[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = glm::vec2(sNormalized, t);
//            }
//        }
//    }
//    P[(ndivs + 1) * ndivs * ncurves] = curveData[curveNumPts - 1];
//    N[(ndivs + 1) * ndivs * ncurves] = (curveData[curveNumPts - 2] - glm::normalize(curveData[curveNumPts - 1]));
//    st[(ndivs + 1) * ndivs * ncurves] = glm::vec2(1, 0.5);
//    uint32_t numFaces = ndivs * ndivs * ncurves;
//    std::unique_ptr<uint32_t []> verts(new uint32_t[numFaces]);
//    for (uint32_t i = 0; i < numFaces; ++i)
//        verts[i] = (i < (numFaces - ndivs)) ? 4 : 3;
//    std::unique_ptr<uint32_t []> vertIndices(new uint32_t[ndivs * ndivs * ncurves * 4 + ndivs * 3]);
//    uint32_t nf = 0, ix = 0;
//    for (uint32_t k = 0; k < ncurves; ++k) {
//        for (uint32_t j = 0; j < ndivs; ++j) {
//            if (k == (ncurves - 1) && j == (ndivs - 1)) { break; }
//            for (uint32_t i = 0; i < ndivs; ++i) {
//                vertIndices[ix] = nf;
//                vertIndices[ix + 1] = nf + (ndivs + 1);
//                vertIndices[ix + 2] = nf + (ndivs + 1) + 1;
//                vertIndices[ix + 3] = nf + 1;
//                ix += 4;
//                ++nf;
//            }
//            nf++;
//        }
//    }
//
//    for (uint32_t i = 0; i < ndivs; ++i) {
//        vertIndices[ix] = nf;
//        vertIndices[ix + 1] = (ndivs + 1) * ndivs * ncurves;
//        vertIndices[ix + 2] = nf + 1;
//        ix += 3;
//        nf++;
//    }
//
//    objects.push_back(std::unique_ptr<MeshTriangle>(new MeshTriangle(
//            st,
//            verts,
//            N,
//            vertIndices,
//            P,
//            glm::mat4(1.0f),
//            )));
//}
//





int main()
{
    YAML::Node config = YAML::LoadFile("config.yaml");
    
    // creating the scene (adding objects and lights)
    std::vector<std::unique_ptr<Hittable>> objects;
    std::vector<std::unique_ptr<Light>> lights;

    auto fileName = config["scene"].as<std::string>();
    YAML::Node input = YAML::LoadFile(fileName + ".yaml");
    objects = loadSceneFromFile(fileName);

    unsigned numOfObjects = objects.size();
    //std::cout << "Loaded " << numOfObjects << " objects" << std::endl;



//    createCurveGeometry(objects);



    lights.push_back(std::make_unique<Light>(glm::vec3(0,4,-11), glm::vec3(1,1,1)));
    lights.push_back(std::make_unique<Light>(glm::vec3(0,0,5), glm::vec3(1,1,1)));

    // setting up options
    Options options{};
    options.width = 300;
    options.height =300;
    options.fov = 80;
    options.backgroundColor = glm::vec3(0.235294, 0.67451, 0.843137);
    options.maxDepth = 2;
    options.bias = 0.00001;
    //cameraToWorld is the inverse of the camera matrix
    //it is the matrix that transforms from world space to camera space
    //we create it by rotating around the y-axis by 180 degrees and then translating by (0,0,5)
    options.cameraToWorld = glm::mat4(1.0f);
//    options.cameraToWorld = glm::translate(options.cameraToWorld, glm::vec3(17,7,15));
//    options.cameraToWorld = glm::rotate(options.cameraToWorld, glm::radians(50.0f), glm::vec3(0,1,0));

    // finally, render
    render(options, objects, lights);

    

    return 0;
}