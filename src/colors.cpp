#include "colors.h"

glm::vec3 colors[10];

void generateColors()
{
    for (int i = 0; i < 10; i++) {
        (void)random_float();
        colors[i] = { random_float(),
            random_float(),
            random_float() };
    }
}
glm::vec3 getColor(int i)
{
    if (i < 0)
        i = 0;
    if (i >= 10)
        i = 0;
    return colors[i];
}