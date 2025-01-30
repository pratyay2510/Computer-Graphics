#include "parse.h"
#include "texture.h"
#include "dump_png.h"
#include "misc.h"
#include <cmath>
#include <algorithm>

Texture::Texture(const Parse* parse, std::istream& in)
{
    std::string filename;
    in >> name >> filename >> use_bilinear_interpolation;
    Read_png(data, width, height, filename.c_str());
}

Texture::~Texture()
{
    delete[] data;
}

// Helper function to wrap floating-point values into the range [0, 1)
inline double Wrap_Float(double value, double max)
{
    double wrapped = std::fmod(value, max);
    if (wrapped < 0) wrapped += max;
    return wrapped;
}

vec3 Texture::Get_Color(const vec2& uv) const
{
    // 1. Wrap texture coordinates to ensure they are in the range [0, 1)
    double u = Wrap_Float(uv[0], 1.0);  // Handle wrapping for u
    double v = Wrap_Float(uv[1], 1.0);  // Handle wrapping for v

    // Debugging: Log texture coordinates
    // Pixel_Print("texture (u,v): (", u, " ", v, ")");

    // 2. Compute the pixel indices
    int i = static_cast<int>(std::floor(u * width)) % width;   // Ensure 0 <= i < width
    int j = static_cast<int>(std::floor(v * height)) % height; // Ensure 0 <= j < height

    // 3. Handle edge cases for negative wrapping
    if (i < 0) i += width;   // Wrap negative indices
    if (j < 0) j += height;  // Wrap negative indices

    // Debugging: Log computed pixel indices
    // Pixel_Print("(i,j): ", i, " ", j);

    // 4. Access the pixel at the calculated index
    const Pixel& pixel = data[j * width + i];

    // 5. Convert the Pixel (assuming it is in ARGB format) to vec3
    // Extract RGB components and normalize to [0, 1] range
    double r = (pixel >> 24) & 0xFF;  // Extract red component
    double g = (pixel >> 16) & 0xFF;  // Extract green component
    double b = (pixel >> 8) & 0xFF;   // Extract blue component

    // Debugging: Log extracted color values
    // Pixel_Print("color: (", r / 255.0, " ", g / 255.0, " ", b / 255.0, ")");

    return vec3(r, g, b) / 255.0;  // Scale RGB to [0, 1] range
}




