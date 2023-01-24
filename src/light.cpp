#include "light.h"
#include "config.h"
#include "common.h"
#include "texture.h"
#include "render.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    //take random float between 0 and 1.
    float random = (float) rand() / RAND_MAX;

    //calculate random point on the segment using the random float.
    glm::vec3 startEnd = segmentLight.endpoint1 - segmentLight.endpoint0;
    position = segmentLight.endpoint0 + random * startEnd;

    //interpolate the color using the same random float.
    color = (1 - random) * segmentLight.color0 + random * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    //take two random floats betweem 0 and 1.
    float r1 = (float) rand() / RAND_MAX;
    float r2 = (float) rand() / RAND_MAX;

    //calculate a random point on the parallelogram using the random floats.
    position = parallelogramLight.v0 + r1 * parallelogramLight.edge01 + r2 * parallelogramLight.edge02;

    //use bilinear interpolation to calculate the color at the point.
    glm::vec3 colorA = (1 - r1) * parallelogramLight.color0 + r1 * parallelogramLight.color1;
    glm::vec3 colorB = (1 - r1) * parallelogramLight.color2 + r1 * parallelogramLight.color3;
    color = (1 - r2) * colorA + r2 * colorB;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    //return clause for disabled Hard Shadows.
    if (!features.enableHardShadow) {
        return 1.0f;
    }

    glm::vec3 point = ray.origin + ray.t * ray.direction;
    glm::vec3 normal = hitInfo.normal / glm::length(hitInfo.normal);

    //calculate on which side of the hit plane the camera (view) and the light are.
    bool viewside = glm::dot(normal, (ray.origin - point) / glm::length(ray.origin - point)) > 0.0f;
    bool lightside = glm::dot(normal, (samplePos - point) / glm::length(samplePos - point)) > 0.0f;

    //point is in shadow if the light shines on a different side of the plane than the viewray hits
    if (viewside != lightside) {
        return 0.0f;
    }

    // calculate the distance from point to the light and the t-value of the ray from the point to the light.
    glm::vec3 lightdir = samplePos - point;
    float dLightPoint = glm::length(lightdir);
    lightdir = lightdir / dLightPoint;
    // move the point slightly to avoid intersecting with itself.
    Ray lightray = Ray(point + 0.000001f * lightdir, lightdir, dLightPoint);
    bvh.intersect(lightray, hitInfo, features);

    //draw the lightray in its correct color if it reaches the light.
    glm::vec3 color = debugColor;
    float result = 1.0f;
    
    //Point is in shadow when the lightray is blocked (t < dLightPoint).
    //lightray.t < dLightPoint compensated for floating point inaccuracy.
    if (lightray.t - dLightPoint < -0.00001) {
        color = { 1, 0, 0 };
        result = 0.0f;
    }

    //restore the t-value for when it is set to infinity by the intersect function.
    if (lightray.t > dLightPoint) {
        lightray.t = dLightPoint;
    }

    drawRay(lightray, color);
    return result;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {

        glm::vec3 totalLight = glm::vec3(0, 0, 0);

        //initialize random-seed used for sampling.
        std::srand(std::chrono::system_clock::now().time_since_epoch().count());

        //loop over all the lights.
        for (const auto& light : scene.lights) {

            if (std::holds_alternative<PointLight>(light)) {
                const PointLight p = std::get<PointLight>(light);

                //calculate the color the point refects after being hit with light from the lightsource.
                glm::vec3 color = computeShading(p.position, p.color, features, ray, hitInfo) * testVisibilityLightSample(p.position, p.color, bvh, features, ray, hitInfo);
                totalLight += color;
            } 
            
            else if (std::holds_alternative<SegmentLight>(light) && features.enableSoftShadow) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);

                glm::vec3 totalLightSampled = glm::vec3(0.0);
                float totalVisible = 0.0;

                //take 20 random samples from the segment.
                for (int i = 0; i <20 ; i++) {
                    glm::vec3 position = glm::vec3(0.0);
                    glm::vec3 color = glm::vec3(0.0);
                    sampleSegmentLight(segmentLight, position, color);

                    //if a sample is not visible from the point, its light should not count towards the average.
                    float visible = testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);
                    glm::vec3 c = computeShading(position, color, features, ray, hitInfo) * visible;
                    totalLightSampled += c;
                    totalVisible += visible;
                }
                //take the average of all visible samples.
                if (totalVisible > 0.0) {
                    totalLight += totalLightSampled / totalVisible;
                }
               
            } 
            
            else if (std::holds_alternative<ParallelogramLight>(light) && features.enableSoftShadow) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);

                glm::vec3 totalLightSampled = glm::vec3(0.0);
                float totalVisible = 0.0;
                
                // take 100 random samples from the Parallelogram.
                for (int i = 0; i < 100; i++) {
                    glm::vec3 position = glm::vec3(0.0);
                    glm::vec3 color = glm::vec3(0.0);
                    sampleParallelogramLight(parallelogramLight, position, color);

                    // if a sample is not visible from the point, its light should not count towards the average.
                    float visible = testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);
                    glm::vec3 c = computeShading(position, color, features, ray, hitInfo) * visible;
                    totalLightSampled += c;
                    totalVisible += visible;
                }
                // take the average of all visible samples.
                if (totalVisible > 0.0) {
                    totalLight += totalLightSampled / totalVisible;
                }
            }
        }

        return totalLight;
    } 
    
    else {
        // Instead of material.kd return acquireTexel when textures are enabled
        if (features.enableTextureMapping && hitInfo.material.kdTexture) {
            return acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        }

        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
