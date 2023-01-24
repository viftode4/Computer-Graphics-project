#include "texture.h"
#include <framework/image.h>
#include <interpolate.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    float x = texCoord.x * image.width;
    float y = texCoord.y * image.height;

    int index = floor(x) * image.width + floor(y);

    if (!features.extra.enableBilinearTextureFiltering) {
        return image.pixels[index];
    }
    
    // Establish where the texCoord lies in comparison
    // to the closest pixel
    glm::vec2 center = { floor(x), floor(y) };

    // Translate the point down by (0.5, 0.5) to 
    // get the center of all texels alligned on the grid
    x = std::max(0.0, x - 0.5);
    y = std::max(0.0, y - 0.5);

    bool right = x > center.x;
    bool above = y > center.y;

    // Default value of all x's and y's is the 
    // same as that of the closest pixel. 
    float x1 = center.x;
    float x2 = center.x;
    float y1 = center.y;
    float y2 = center.y;

    glm::vec3 x1y1, x2y1, x1y2, x2y2;

    if (right) {
        x2 = std::min((float)image.width, x2 + 1);
    }

    // Left
    else {
        x1 = std::max(0.0, x1 - 1.0);
    }

    if (above) {
        y2 = std::min((float)image.width, y2 + 1);
    }

    // Beneath
    else {
        y1 = std::max(0.0, y1 - 1.0);
    }

    x1y1 = image.pixels[x1 * image.width + y1];
    x2y1 = image.pixels[x2 * image.width + y1];
    x1y2 = image.pixels[x1 * image.width + y2];
    x2y2 = image.pixels[x2 * image.width + y2];

    return bilinearInterpolation(x1y1, x2y1, x1y2, x2y2, x1, x2, x, y1, y2, y);
}