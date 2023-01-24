#include "texture.h"
#include "render.h"
#include "randoms.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include "light.h"

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 hitPoint = ray.origin + ray.direction * ray.t;
    glm::vec3 I = lightColor;

    glm::vec3 Kd = hitInfo.material.kd;
    glm::vec3 Ks = hitInfo.material.ks;

    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(lightPosition - hitPoint);
    glm::vec3 R = normalize(2.0f * (N * dot(N, L) - L) + L);
    glm::vec3 V = normalize(ray.origin - hitPoint);

    float s = hitInfo.material.shininess;
    
    if (hitInfo.material.kdTexture && features.enableTextureMapping) {
        Kd = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
    }

    glm::vec3 D = I * Kd * dot(N, L);

    glm::vec3 S = I * Ks * pow(dot(R, V), s);

    glm::vec3 result = glm::vec3(0.0);

    if (dot(N, ray.direction) >= 0) {
        return result;
    }

    if (dot(N, L) >= 0.0) {
        result += D;
    }

    if (dot(R, V) >= 0.0) {
        result += S;
    }

    return result;
}


const Ray computeReflectionRay(Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    if (hitInfo.material.ks == glm::vec3(0.0)) {
        return Ray { glm::vec3(0.0), glm::vec3(0.0), -1 };
    }

    Ray reflectionRay;

    
    reflectionRay.debugColor = ray.debugColor + 1;
    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(-1.0f * ray.direction);
    reflectionRay.direction = normalize(2.0f * (N * dot(N, L) - L) + L);

    reflectionRay.origin = ray.origin + ray.direction * ray.t + reflectionRay.direction * 0.0001f;

    reflectionRay.t = std::numeric_limits<float>::max();

    return reflectionRay;
}

//calculates a pertubed ray based on the perfect reflection.
const Ray computeGlossyReflectionRay(Ray reflection)
{
    //draw the unpertubed reflection in red.
    drawRay(reflection, glm::vec3(1, 0, 0));
    
    //create a basis for the sqaure perpendicular to the reflection ray.
    glm::vec3 t = { 1, 0, 0 };
    glm::vec3 u = glm::normalize(glm::cross(t, reflection.direction));
    glm::vec3 v = glm::cross(reflection.direction, u);

    //define the size of the sqaure, and thus the magnatude of the reflection glossyness.  
    float a = 0.5f;

    //take two random floats between 0 and 1;
    float r1 = (float)rand() / RAND_MAX;
    float r2 = (float)rand() / RAND_MAX;

    //calculate the scalars for pertubing the reflection in each direction. 
    float p = -0.5f * a + r1 * a;
    float q = -0.5f * a + r2 * a;

    //calculate the new reflection
    glm::vec3 r = glm::normalize(reflection.direction + p * u + q * v);
    Ray newReflection = Ray(reflection.origin, r, reflection.t);

    return newReflection;
}