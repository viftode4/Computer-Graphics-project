#include "randoms.h"
#include <glm/geometric.hpp>

float random_float()
{
    static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    static std::mt19937 generator;
    return distribution(generator);
}

//return a random vector in a unit circle lens (-1,1), (-1,1)
glm::vec3 random_unit()
{
    return glm::normalize(glm::vec3(random_float() * 2.0f - 1.0, random_float() * 2.0f - 1.0, 0)) * random_float();
}
