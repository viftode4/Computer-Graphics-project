#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 A = v1 - v0;
    glm::vec3 B = v2 - v0;
    glm::vec3 C = p - v0;

    float dAA = dot(A, A);
    float dAB = dot(A, B);
    float dAC = dot(A, C);
    float dBB = dot(B, B);
    float dBC = dot(B, C);

    float D = dAA * dBB - dAB * dAB;

    glm::vec3 bary;

    bary.y = (dBB * dAC - dAB * dBC) / D;     // alpha
    bary.z = (dAA * dBC - dAB * dAC) / D;     // beta
    bary.x = 1.0f - bary.y - bary.z;          // 1 - alpha - beta

    return bary;
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return n0 * barycentricCoord.x + n1 * barycentricCoord.y + n2 * barycentricCoord.z;
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return t0 * barycentricCoord.x + t1 * barycentricCoord.y + t2 * barycentricCoord.z;
}

glm::vec3 bilinearInterpolation(const glm::vec3& x1y1, const glm::vec3& x2y1, const glm::vec3& x1y2, const glm::vec3& x2y2,
    float& x1, float& x2, float& x, float& y1, float& y2, float& y)
{
    if (x1 == x2 && y1 == y2) {
        return x1y1;
    }

    if (x1 == x2 && y1 != y2) {
        return (y2 - y) / (y2 - y1) * x1y1 + (y - y1) / (y2 - y1) * x1y2;
    }

    if (x1 != x2 && y1 == y2) {
        return (x2 - x) / (x2 - x1) * x1y1 + (x - x1) / (x2 - x1) * x2y1;
    }

    glm::vec3 xDirection1 = (x2 - x) / (x2 - x1) * x1y1 + (x - x1) / (x2 - x1) * x2y1;
    glm::vec3 xDirection2 = (x2 - x) / (x2 - x1) * x1y2 + (x - x1) / (x2 - x1) * x2y2;

    return (y2 - y) / (y2 - y1) * xDirection1 + (y - y1) / (y2 - y1) * xDirection2;
}