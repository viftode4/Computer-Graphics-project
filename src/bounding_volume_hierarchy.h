#pragma once
#include "common.h"
#include "colors.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct Node {
    AxisAlignedBox box;
    bool leaf;
    std::vector<std::pair<int, int>> indx;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    void setDebugTraversal(bool x);

    void draw(const int index, const int& level, const int depth);

    void build(std::vector<Node>& bvhVec, std::vector<std::pair<int, int>> triangles, int root, int depth);

    glm::vec3 getAxis(int i);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    

private:
    // x axis, y axis, z axis
    glm::vec3 vecAxis[3] = { { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f } };
    std::vector<int> leaves;
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> bvh;
};
