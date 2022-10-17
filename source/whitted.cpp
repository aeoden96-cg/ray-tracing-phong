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



struct Options
{
    uint32_t width;
    uint32_t height;
    float fov;
    float imageAspectRatio;
    uint8_t maxDepth;
    glm::vec3 backgroundColor;
    float bias;
};

class Light
{
public:
    Light(const glm::vec3 &p, const glm::vec3 &i) : position(p), intensity(i) {}
    glm::vec3 position;
    glm::vec3 intensity;
};





// [comment]
// Compute reflection direction
// [/comment]
glm::vec3 reflect(const glm::vec3 &I, const glm::vec3 &N)
{
    return I - 2 * glm::dot(I, N) * N;
}

// [comment]
// Compute refraction direction using Snell's law
//
// We need to handle with care the two possible situations:
//
//    - When the ray is inside the Hittable
//
//    - When the ray is outside.
//
// If the ray is outside, you need to make cosi positive cosi = -N.I
//
// If the ray is inside, you need to invert the refractive indices and negate the normal N
// [/comment]
glm::vec3 refract(const glm::vec3 &I, const glm::vec3 &N, const float &ior)
{
    float cosi = clamp_(-1, 1, glm::dot(I, N));
    float etai = 1, etat = ior;
    glm::vec3 n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? glm::vec3(0,0,0) : eta * I + (eta * cosi - sqrtf(k)) * n;
}

// [comment]
// Compute Fresnel equation
//
// \param I is the incident view direction
//
// \param N is the normal at the intersection point
//
// \param ior is the mateural refractive index
//
// \param[out] kr is the amount of light reflected
// [/comment]
void fresnel(const glm::vec3 &I, const glm::vec3 &N, const float &ior, float &kr)
{
    float cosi = clamp_(-1, 1, glm::dot(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

// [comment]
// Returns true if the ray intersects an Hittable, false otherwise.
//
// \param orig is the ray origin
//
// \param dir is the ray direction
//
// \param objects is the list of objects the scene contains
//
// \param[out] tNear contains the distance to the cloesest intersected Hittable.
//
// \param[out] index stores the index of the intersect triangle if the interesected Hittable is a mesh.
//
// \param[out] uv stores the u and v barycentric coordinates of the intersected point
//
// \param[out] *hitObject stores the pointer to the intersected Hittable (used to retrieve material information, etc.)
//
// \param isShadowRay is it a shadow ray. We can return from the function sooner as soon as we have found a hit.
// [/comment]
bool trace(
    const glm::vec3 &orig, const glm::vec3 &dir,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    float &tNear, uint32_t &index, glm::vec2 &uv, Hittable **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        glm::vec2 uvK;
        if (objects[k]->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear) {
            *hitObject = objects[k].get();
            tNear = tNearK;
            index = indexK;
            uv = uvK;
        }
    }

    return (*hitObject != nullptr);
}

// [comment]
// Implementation of the Whitted-syle light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected Hittable is either reflective or reflective and refractive,
// then we compute the reflection/refracton direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refractin depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is duffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
// [/comment]
glm::vec3 castRay(
    const glm::vec3 &orig, const glm::vec3 &dir,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights,
    const Options &options,
    uint32_t depth,
    bool test = false)
{
    if (depth > options.maxDepth) {
        return options.backgroundColor;
    }
    
    glm::vec3 hitColor = options.backgroundColor;
    float tnear = kInfinity;
    glm::vec2 uv;
    uint32_t index = 0;
    Hittable *hitObject = nullptr;
    if (trace(orig, dir, objects, tnear, index, uv, &hitObject)) {
        glm::vec3 hitPoint = orig + dir * tnear;
        glm::vec3 N; // normal
        glm::vec2 st; // st coordinates
        hitObject->getSurfaceProperties(hitPoint, dir, index, uv, N, st);
        glm::vec3 tmp = hitPoint;
        switch (hitObject->materialType) {
            case REFLECTION_AND_REFRACTION:
            {
                glm::vec3 reflectionDirection = normalize(reflect(dir, N));
                glm::vec3 refractionDirection = normalize(refract(dir, N, hitObject->ior));
                glm::vec3 reflectionRayOrig = (glm::dot(reflectionDirection, N) < 0) ?
                    hitPoint - N * options.bias :
                    hitPoint + N * options.bias;
                glm::vec3 refractionRayOrig = (glm::dot(refractionDirection, N) < 0) ?
                    hitPoint - N * options.bias :
                    hitPoint + N * options.bias;
                glm::vec3 reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1, 1);
                glm::vec3 refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1, 1);
                float kr;
                fresnel(dir, N, hitObject->ior, kr);
                hitColor = reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            case REFLECTION:
            {
                float kr;
                fresnel(dir, N, hitObject->ior, kr);
                glm::vec3 reflectionDirection = reflect(dir, N);
                glm::vec3 reflectionRayOrig = (glm::dot(reflectionDirection, N) < 0) ?
                    hitPoint + N * options.bias :
                    hitPoint - N * options.bias;
                hitColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1) * kr;
                break;
            }
            default:
            {
                // [comment]
                // We use the Phong illumation model int the default case. The phong model
                // is composed of a diffuse and a specular reflection component.
                // [/comment]
                glm::vec3 lightAmt = glm::vec3(0,0,0), specularColor =  glm::vec3(0,0,0);
                glm::vec3 shadowPointOrig = (glm::dot(dir, N) < 0) ?
                    hitPoint + N * options.bias :
                    hitPoint - N * options.bias;
                // [comment]
                // Loop over all lights in the scene and sum their contribution up
                // We also apply the lambert cosine law here though we haven't explained yet what this means.
                // [/comment]
                for (uint32_t i = 0; i < lights.size(); ++i) {
                    glm::vec3 lightDir = lights[i]->position - hitPoint;
                    // square of the distance between hitPoint and the light
                    float lightDistance2 = glm::dot(lightDir, lightDir);
                    lightDir = normalize(lightDir);
                    float LdotN = std::max(0.f, glm::dot(lightDir, N));
                    Hittable *shadowHitObject = nullptr;
                    float tNearShadow = kInfinity;
                    // is the point in shadow, and is the nearest occluding Hittable closer to the Hittable than the light itself?
                    bool inShadow = trace(shadowPointOrig, lightDir, objects, tNearShadow, index, uv, &shadowHitObject) &&
                        tNearShadow * tNearShadow < lightDistance2;

                    if(!inShadow) {
                       lightAmt += lights[i]->intensity * LdotN;
                    }


                    glm::vec3 reflectionDirection = reflect(-lightDir, N);
                    specularColor += powf(std::max(0.f, -glm::dot(reflectionDirection, dir)), hitObject->specularExponent) * lights[i]->intensity;
                }
                hitColor = lightAmt * hitObject->evalDiffuseColor(st) * hitObject->Kd + specularColor * hitObject->Ks;
                break;
            }
        }
    }

    return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void render(
    const Options &options,
    const std::vector<std::unique_ptr<Hittable>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights)
{
    glm::vec3 *framebuffer = new glm::vec3[options.width * options.height];
    glm::vec3 *pix = framebuffer;
    float scale = tan(glm::radians(options.fov * 0.5));
    float imageAspectRatio = options.width / (float)options.height;
    glm::vec3 orig(0);
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
            glm::vec3 dir = normalize(glm::vec3(x, y, -1));
            *(pix++) = castRay(orig, dir, objects, lights, options, 0);
        }
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

    delete [] framebuffer;
}

// [comment]
// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image widht and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().
// [/comment]
int main(int argc, char **argv)
{
    // creating the scene (adding objects and lights)
    std::vector<std::unique_ptr<Hittable>> objects;
    std::vector<std::unique_ptr<Light>> lights;
    
    Sphere *sph1 = new Sphere(glm::vec3(-1, 0, -12), 2);
    sph1->materialType = DIFFUSE_AND_GLOSSY;
    sph1->diffuseColor = glm::vec3(0.6, 0.7, 0.8);
    Sphere *sph2 = new Sphere(glm::vec3(0.5, -0.5, -8), 1.5);
    sph2->ior = 1.5;
    sph2->materialType = REFLECTION_AND_REFRACTION;
    
    objects.push_back(std::unique_ptr<Sphere>(sph1));
    objects.push_back(std::unique_ptr<Sphere>(sph2));

    glm::vec3 verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}};
    uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3};
    glm::vec2 st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    MeshTriangle *mesh = new MeshTriangle(verts, vertIndex, 2, st);
    mesh->materialType = DIFFUSE_AND_GLOSSY;
    
    objects.push_back(std::unique_ptr<MeshTriangle>(mesh));

    lights.push_back(std::unique_ptr<Light>(new Light(glm::vec3(-20, 70, 20), glm::vec3(0.5f, 0.5f, 0.5f))));
    lights.push_back(std::unique_ptr<Light>(new Light(glm::vec3(30, 50, -12), glm::vec3(1, 1, 1))));
    
    // setting up options
    Options options;
    options.width = 640;
    options.height = 480;
    options.fov = 90;
    options.backgroundColor = glm::vec3(0.235294, 0.67451, 0.843137);
    options.maxDepth = 5;
    options.bias = 0.00001;
    
    // finally, render
    render(options, objects, lights);

    return 0;
}