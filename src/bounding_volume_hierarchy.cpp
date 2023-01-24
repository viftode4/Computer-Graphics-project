#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.cpp"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include "randoms.h"
#include <glm/glm.hpp>
#include <queue>
#include <stack>
#include <iostream>

bool debugBVH = 0;
// AABB helpers
AxisAlignedBox combine(const AxisAlignedBox& a, const AxisAlignedBox& b)
{
    AxisAlignedBox box;

    box.lower.x = std::min(a.lower.x, b.lower.x);
    box.lower.y = std::min(a.lower.y, b.lower.y);
    box.lower.z = std::min(a.lower.z, b.lower.z);

    box.upper.x = std::max(a.upper.x, b.upper.x);
    box.upper.y = std::max(a.upper.y, b.upper.y);
    box.upper.z = std::max(a.upper.z, b.upper.z);

    return box;
}

AxisAlignedBox getTriangleBox(const Vertex& a, const Vertex& b, const Vertex& c)
{
    AxisAlignedBox box;
    box.lower.x = std::min({ a.position.x, b.position.x, c.position.x });
    box.lower.y = std::min({ a.position.y, b.position.y, c.position.y });
    box.lower.z = std::min({ a.position.z, b.position.z, c.position.z });

    box.upper.x = std::max({ a.position.x, b.position.x, c.position.x });
    box.upper.y = std::max({ a.position.y, b.position.y, c.position.y });
    box.upper.z = std::max({ a.position.z, b.position.z, c.position.z });

    return box;
}

AxisAlignedBox getMeshesAABB(const std::vector<std::pair<int, int>>& triangles, const std::vector<Mesh>& meshes)
{
    AxisAlignedBox allBox;
    bool firstBox = true;

    for (int i = 0; i < triangles.size(); i++) {

        if (firstBox) {
            allBox = getTriangleBox(meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].x],
                meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].y],
                meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].z]);
            firstBox = false;
        }
        auto aabb = getTriangleBox(meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].x],
            meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].y],
            meshes[triangles[i].first].vertices[meshes[triangles[i].first].triangles[triangles[i].second].z]);
        allBox = combine(allBox, aabb);
    }

    return allBox;
}

glm::vec3 BoundingVolumeHierarchy::getAxis(int i) {
    return vecAxis[i];
}
// bvh builder
void BoundingVolumeHierarchy::build(std::vector<Node>& bvhVec, std::vector<std::pair<int,int>> triangles, int root, int depth)
{

    bvhVec[root].indx.resize(0);
    bvhVec[root].indx.shrink_to_fit();
    // first we get the AABB bounding box for the triangles in the curent node
    bvhVec[root].box = getMeshesAABB(triangles, m_pScene->meshes);
    // we first assume this node is not a leaf
    bvhVec[root].leaf = false;

    // formula for depth of BST from index
    // we stop when the depth exceeds 20 (computation waste)
    // the node that reaches this is automatically a leaf
    // if depth is 20 (big) because we don't lose that much speed compared to cost
    if (depth >= 20 || triangles.size() <= 1) {
        bvhVec[root].leaf = true;
        bvhVec[root].indx = triangles;
        m_numLeaves++;
        leaves.push_back(root);

        return;
    }
    bvhVec[root].indx.push_back({ 2 * root + 1, 2 * root + 2 });
    // for spliting by median we need to sort by a custom
    // comparator by each axis
    // this way we can easily pick our axis;

    auto axis = (depth % 3 == 0) ? 0 : (depth % 3 == 1) ? 1 : 2;

    std::vector<std::pair<float, int>> centroids;
    std::vector<std::pair<int, int>> leftTriangles, rightTriangles;

    /*centroids.resize(0);
    centroids.shrink_to_fit();

    leftTriangles.resize(0);
    leftTriangles.shrink_to_fit();

    rightTriangles.resize(0);
    rightTriangles.shrink_to_fit();*/

    for (int i = 0; i < triangles.size(); i++) 
        centroids.push_back({ glm::dot(getAxis(axis), this->m_pScene->meshes[triangles[i].first].getCentroid(triangles[i].second)), i});
    
    std::sort(centroids.begin(), centroids.end());

    for (int i = 0; i < centroids.size(); i++)
        if (i < centroids.size() / 2)
            leftTriangles.push_back(triangles[centroids[i].second]);
        else rightTriangles.push_back(triangles[centroids[i].second]);

    triangles.resize(0);
    triangles.shrink_to_fit();

    centroids.resize(0);
    centroids.shrink_to_fit();

    build(bvhVec, leftTriangles, 2 * root + 1, depth + 1);
    build(bvhVec, rightTriangles, 2 * root + 2,  depth + 1);
}
// returns a vector of (index mesh, index triangle in mesh) pairs
void getTriangles(std::vector<std::pair<int, int>>& sol, const std::vector<Mesh>& meshes)
{

    for (int i = 0; i < meshes.size(); i++)
        for (int j = 0; j < meshes[i].triangles.size(); j++)
            sol.push_back({ i, j });
}

//above are helper functions
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene): m_pScene(pScene)
{   
    // getting all the triangles in the format:
    // pair<int,int> x  
    // x.first -> the index of the mesh
    // x.second -> index of the triangle

    this->bvh.resize(0);
    this->bvh.shrink_to_fit();

    leaves.resize(0);
    leaves.shrink_to_fit();

    std::vector<std::pair<int, int>> triangles;

    getTriangles(triangles, pScene->meshes);

    m_numLevels = std::min(20.0, ceil(log2(triangles.size()))) + 1;
    m_numLeaves = 0;

    this->bvh.resize((1 << (m_numLevels+1))-1);

    build(this->bvh, triangles, 0, 0);
}


// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{

   /* level = std::min(level, m_numLevels);

    int startingIndex = 0;
    for (int i = 0; i < level; i++) {
        startingIndex += pow(2, i);
    }

    int cap = startingIndex + pow(2, level);

    AxisAlignedBox aabb;

    for (startingIndex; startingIndex < cap; startingIndex++) {
        aabb = { bvh[startingIndex].box.lower, bvh[startingIndex].box.upper };

        if (bvh[startingIndex].indx.size())
            drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.8f);
    }*/



    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    draw(0, level, 0);
    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
   // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}
void BoundingVolumeHierarchy::draw(const int index, const int& level, const int depth)
{
    if (depth >= level) {
        if (this->bvh[index].indx.size())
        drawAABB(this->bvh[index].box, DrawMode::Wireframe);
        return;
    }

    if (this->bvh[index].leaf)
        return;

    draw(this->bvh[index].indx[0].first, level, depth + 1);
    draw(this->bvh[index].indx[0].second, level, depth + 1);
}
 
// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    if (leafIdx <= 0)
        leafIdx++;
    if (leafIdx >= leaves.size())
        leafIdx = leaves.size();
    drawAABB(bvh[leaves[leafIdx-1]].box, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f));
    for (auto it : bvh[leaves[leafIdx-1]].indx)
        drawTriangle(Vertex(m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second].x].position),
            Vertex(m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second].y].position),
            Vertex(m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second].z].position));
    
    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


void drawNodes(std::stack<AxisAlignedBox>& nodes, const std::vector<glm::vec3>& triangle, int color ) {
    while (!nodes.empty()) {
        drawAABB(nodes.top(), DrawMode::Wireframe, getColor(color));
        nodes.pop();
    }
    if (triangle.size() == 3)
        drawTriangle(Vertex(triangle[0]), Vertex(triangle[1]), Vertex(triangle[2]));
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                    hitInfo.normal = cross((v1.position - v0.position), (v2.position - v0.position));

                    if (features.enableNormalInterp) {
                        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);

                        // Visual debug
                        Ray normalRay = Ray { ray.origin + ray.direction * ray.t, hitInfo.normal, 10 };
                        Ray Ray0 = Ray { v0.position, v0.normal, 5 };
                        Ray Ray1 = Ray { v1.position, v1.normal, 5 };
                        Ray Ray2 = Ray { v2.position, v2.normal, 5 };
                        drawRay(normalRay, glm::vec3(1.0, 0, 0));
                        drawRay(Ray0, glm::vec3(0, 1.0, 0));
                        drawRay(Ray1, glm::vec3(0, 1.0, 0));
                        drawRay(Ray2, glm::vec3(0, 1.0, 0));
                    }

                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        if (!intersectRayWithShape(this->bvh[0].box, ray))
            return false;

        std::priority_queue<std::pair<float, int>> queue;
        std::stack<AxisAlignedBox> intersected;
        std::vector<glm::vec3> hitTriangle;

        // for bvh debug
        float oldT = FLT_MAX;
        ray.t = FLT_MAX;
        bool hit = 0;
        //insert root doesnt really matter
        queue.push({ FLT_MAX, 0 });
        
        while (!queue.empty()) {
            float t = -queue.top().first;
            int root = queue.top().second;

            queue.pop();

            if (t > ray.t)
                break;

            if (debugBVH)
                intersected.push(bvh[root].box);

            if (!bvh[root].leaf) {
                oldT = ray.t;
                ray.t = FLT_MAX;

                if (intersectRayWithShape(this->bvh[root*2+1].box, ray)) {
                    queue.push({ -ray.t, root * 2 + 1});
                }

                ray.t = FLT_MAX;
                if (intersectRayWithShape(this->bvh[root*2+2].box, ray)) {
                    queue.push({ -ray.t, root*2 + 2});
                }
                ray.t = oldT;
                continue;
            }

            oldT = ray.t;

            for (int i = 0; i < bvh[root].indx.size(); i++) {
                const auto& triangle = bvh[root].indx[i];
                const auto& mesh = m_pScene->meshes[triangle.first];

                ray.t = FLT_MAX;

                const auto v0 = mesh.vertices[mesh.triangles[triangle.second].x];
                const auto v1 = mesh.vertices[mesh.triangles[triangle.second].y];
                const auto v2 = mesh.vertices[mesh.triangles[triangle.second].z];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo))
                    if (oldT >= ray.t) {
                        hit = 1;

                        hitInfo.material = mesh.material;
                        hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);
                        hitInfo.texCoord = interpolateTexCoord(v0.position, v1.position, v2.position, hitInfo.barycentricCoord);
                        hitInfo.normal = cross((v1.position - v0.position), (v2.position - v0.position));
                        
                        if (features.enableNormalInterp) {
                            hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);

                            // Visual debug
                            Ray normalRay = Ray { ray.origin + ray.direction * ray.t, hitInfo.normal, 10 };
                            Ray Ray0 = Ray { v0.position, v0.normal, 5 };
                            Ray Ray1 = Ray { v1.position, v1.normal, 5 };
                            Ray Ray2 = Ray { v2.position, v2.normal, 5 };
                            drawRay(normalRay, glm::vec3(1.0, 0, 0));
                            drawRay(Ray0, glm::vec3(0, 1.0, 0));
                            drawRay(Ray1, glm::vec3(0, 1.0, 0));
                            drawRay(Ray2, glm::vec3(0, 1.0, 0));
                    }

                        if (debugBVH) 
                            hitTriangle = { mesh.vertices[mesh.triangles[triangle.second].x].position,
                                mesh.vertices[mesh.triangles[triangle.second].y].position,
                                mesh.vertices[mesh.triangles[triangle.second].z].position };
                        
                        oldT = ray.t;
                    }
            }
            ray.t = oldT;

        }

        if (debugBVH)
            drawNodes(intersected, hitTriangle, ray.debugColor);

        return hit;
    }
}