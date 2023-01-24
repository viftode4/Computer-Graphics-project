#pragma once
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/quaternion.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include "ray.h"

class Window;

class Trackball {
public:
	// NOTE(Mathijs): field of view in radians! (use glm::radians(...) to convert from degrees to radians).
	Trackball(Window* pWindow, float fovy, float distanceFromLookAt = 4.0f, float rotationX = 0.0f, float rotationY = 0.0f);
	Trackball(Window* pWindow, float fovy, const glm::vec3& lookAt, float distanceFromLookAt = 4.0f, float rotationX = 0.0f, float rotationY = 0.0f);
	~Trackball() = default;

	static void printHelp();

	void disableTranslation();

	[[nodiscard]] glm::vec3 left() const;
	[[nodiscard]] glm::vec3 up() const;
	[[nodiscard]] glm::vec3 forward() const;

	[[nodiscard]] glm::vec3 position() const; // Position of the camera.
    [[nodiscard]] glm::vec3 positionWithOffset(const glm::vec3& offset) const; // Position of the camera with offset.
	[[nodiscard]] glm::vec3 lookAt() const; // Point that the camera is looking at / rotating around.
	[[nodiscard]] glm::mat4 viewMatrix() const;
	[[nodiscard]] glm::mat4 projectionMatrix() const;
        [[nodiscard]] glm::vec3 rotationEulerAngles() const;
        [[nodiscard]] float distanceFromLookAt() const;
        [[nodiscard]] int getSamples() const;
        [[nodiscard]] float getFocusDistance() const;
        [[nodiscard]] float getLensRadius() const;
        [[nodiscard]] bool getAlgorithm() const;

	void setCamera(const glm::vec3 lookAt, const glm::vec3 rotations, const float dist); // Set the position and orientation of the camera.
	
	//for DOF to be able to change focus range
	void setDOFradius(const float r);

	// for DOF to be able to change focus distance
    void setDOFfocusDistance(const float r);

	//for Antialising to set the number of samples we want per pixel
	void setSamples(const int no);
	
	void setAlgorithm(const bool x);

	// Generate ray given pixel in NDC space (ranging from -1 to +1. (-1,-1) at bottom left, (+1, +1) at top right).
	[[nodiscard]] Ray generateRay(const glm::vec2& pixel) const;

private:

	void mouseButtonCallback(int button, int action, int mods);
	void mouseMoveCallback(const glm::vec2& pos);
	void mouseScrollCallback(const glm::vec2& offset);

private:
	// for DOF
    // focus distance is distance from look at let's say.
	// the blue ball in raster view
	float lens_radius = 0;
    float focus_distance = 1;
	//for Antialising
    int samples = 1;
    bool random = 0;

	const Window* m_pWindow;
	float m_fovy;
	float m_halfScreenSpaceHeight;
	float m_halfScreenSpaceWidth;
	bool m_canTranslate { true };

	glm::vec3 m_lookAt{ 0.0f }; // Point that the camera is looking at / rotating around.
	float m_distanceFromLookAt;
	glm::vec3 m_rotationEulerAngles{ 0 }; // Rotation as euler angles (in radians).

	glm::vec2 m_prevCursorPos; // Cursor position.
};
