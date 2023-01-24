#include "render.h"
#include "intersect.h"
#include "light.h"
#include "randoms.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

bool AntialisingAlgorithm = 0;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    hitInfo.normal = glm::vec3 { 0 };
    hitInfo.barycentricCoord = glm::vec3 { 0 };
    hitInfo.texCoord = glm::vec2 { 0 };

    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        //calculate reflections.
        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);

            //calculate gloss only for relfective surfaces.
            if (features.extra.enableGlossyReflection && reflection.t != -1) {

                //set the random seed for sampling.
                std::srand(std::chrono::system_clock::now().time_since_epoch().count());

                //calculate 1000 perturbed reflections.
                glm::vec3 totalColor = { 0, 0, 0 };
                for (int i = 0; i < 1000; i++) {
                    Ray r = computeGlossyReflectionRay(reflection);
                    bool hit = bvh.intersect(r, hitInfo, features);
                    if (r.t >= 0) {
                        drawRay(r, { 0, 1, 0 });
                        if (hit) {
                            totalColor += computeLightContribution(scene, bvh, features, r, hitInfo);
                            
                        }
                    }
                }
                //average the color of all rays. 
                Lo += totalColor / 1000.0f;

            //perfect reflection
            } else {
                if (reflection.direction != glm::vec3 { 0.0f }) {
                    bool hit = bvh.intersect(reflection, hitInfo, features);

                    if (reflection.t >= 0) {

                        // There is a reflection (Ks is not black)
                        drawRay(reflection, { 0, 1, 0 });

                        if (hit) {
                            Lo += computeLightContribution(scene, bvh, features, reflection, hitInfo);
                        }
                    }
                }
                
            }
        }

        // When shading is enabled color the ray the color
        // of the point it hits
        if (features.enableShading) {
            drawRay(ray, Lo);
        }

        // When shading is disabled color the ray the color
        else {
            drawRay(ray);
        }

        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {

            if (!features.extra.enableDepthOfField && !features.extra.enableMultipleRaysPerPixel) {
                // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
                const glm::vec2 normalizedPixelPos {
                    float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y) / float(windowResolution.y) * 2.0f - 1.0f
                };
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
                continue;
            }

            glm::vec3 color(0.0f);
            // option algorithm 1 random sampling but inside a grid of samples x samples which can be set form UI
            if (camera.getAlgorithm() == 1) {

                for (int i = 0; i < camera.getSamples(); i++)
                    for (int j = 0; j < camera.getSamples(); j++) {
                        const glm::vec2 normalizedPixelPos {
                            float(x + (i + random_float()) / camera.getSamples()) / float(windowResolution.x) * 2.0f - 1.0f,
                            float(y + (j + random_float()) / camera.getSamples()) / float(windowResolution.y) * 2.0f - 1.0f
                        };
                        const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                        color += getFinalColor(scene, bvh, cameraRay, features);
                    }

            }
            // option algorithm 0 truly random sampling 
            else {

                for (int i = 0; i < camera.getSamples(); i++) {

                    const glm::vec2 normalizedPixelPos {
                        float(x + random_float()) / float(windowResolution.x) * 2.0f - 1.0f,
                        float(y + random_float()) / float(windowResolution.y) * 2.0f - 1.0f
                    };
                    const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                    color += getFinalColor(scene, bvh, cameraRay, features);
                }

            }
           

            // we get the color that needs to be averaged to the number of samples
            screen.setPixel(x, y, color, camera.getSamples(), camera.getAlgorithm());
        }
    }
}
