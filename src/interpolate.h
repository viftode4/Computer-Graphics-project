#pragma once

#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()

// Given the three triangle vertices and a point, compute the corresponding barycentric coordinates of the point
// returns a vec3 with the barycentric coordinate (alpha, beta, gamma)
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p);

// Interpolate three normals using barycentric coordinates
// returns the interpolated normal
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord);

// Interpolate three texture coordinates using barycentric coordinates
// returns the interpolated texture coordinates
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord);

// First do linear interpolation in the x-direction for both y values. Then interpolate over the y-direction with the results
// of the first two interpolations
glm::vec3 bilinearInterpolation(const glm::vec3& x1y1, const glm::vec3& x2y1, const glm::vec3& x1y2, const glm::vec3& x2y2,
    float& x1, float& x2, float& x, float& y1, float& y2, float& y);
