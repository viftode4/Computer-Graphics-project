#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    float a, b, c;

    a = glm::dot(glm::cross(v2 - v1, p - v1), n) / glm::dot(glm::cross(v1 - v0, v2 - v0), n);
    b = glm::dot(glm::cross(v0 - v2, p - v2), n) / glm::dot(glm::cross(v1 - v0, v2 - v0), n);
    c = glm::dot(glm::cross(v1 - v0, p - v0), n) / glm::dot(glm::cross(v1 - v0, v2 - v0), n);

    if (a < 0 || b < 0 || a + b > 1)
        return false;
    return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float a, b;

    b = glm::dot(plane.normal, ray.direction);
    if (!b)
        return false;

    a = plane.D - glm::dot(plane.normal, ray.origin);
    if (a / b < 0)
        return false;

    ray.t = a / b;
    return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 n = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    plane.normal = n;
    plane.D = glm::dot(v0, n);

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    float t = ray.t;

    if (!intersectRayWithPlane(trianglePlane(v0, v1, v2), ray)) {
        ray.t = t;
        return false;
    }

    if (!pointInTriangle(v0, v1, v2, trianglePlane(v0, v1, v2).normal, ray.origin + ray.direction * ray.t)) {
        ray.t = t;
        return false;
    }

    if (ray.t > t) {
        ray.t = t;
        return false;
    }

    return true;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float radius = sphere.radius;

    float a = glm::dot(ray.direction, ray.direction);
    float b = 2.0f * glm::dot(ray.direction, ray.origin - sphere.center);
    float c = glm ::dot(ray.origin - sphere.center, ray.origin - sphere.center) - radius * radius;

    float d = b * b - 4.0f * a * c;

    if (d < 0)
        return false;

    float t1 = (-b + std::sqrt(d)) / (2.0f * a);
    float t2 = (-b - std::sqrt(d)) / (2.0f * a);

    if (std::abs(t2) < std::abs(t1))
        t1 = t2;

    if (ray.t > t1)
        ray.t = t1;

    return true;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float tmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tmax = (box.upper.x - ray.origin.x) / ray.direction.x;

    if (tmin > tmax)
        std::swap(tmin, tmax);

    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;

    if (tymin > tymax)
        std::swap(tymin, tymax);

    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    if (tzmin > tzmax)
        std::swap(tzmin, tzmax);

    tmin = std::max({ tmin, tymin, tzmin });
    tmax = std::min({ tmax, tymax, tzmax });

    if (tmin > tmax || tmax < 0)
        return false;

    //if (tmin < 0 && tmax > 0)
    //    ray.t = 0.0f;

    if (ray.t > tmin)
        ray.t = tmin;

    return true;
}
