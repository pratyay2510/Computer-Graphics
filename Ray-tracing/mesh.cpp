#include "mesh.h"
#include <fstream>
#include <limits>
#include <string>
#include <algorithm>
#include <cassert>

static const double weight_tolerance = 1e-4;

Mesh::Mesh(const Parse* parse, std::istream& in)
{
    std::string file;
    in >> name >> file;
    Read_Obj(file.c_str());
}

// Read in a mesh from an obj file. Populates the bounding box and registers
// one part per triangle (by setting number_parts).
void Mesh::Read_Obj(const char* file)
{
    std::ifstream fin(file);
    if (!fin)
    {
        exit(EXIT_FAILURE);
    }

    std::string line;
    ivec3 e, t;
    vec3 v;
    vec2 u;

    while (getline(fin, line))
    {
        if (sscanf(line.c_str(), "v %lg %lg %lg", &v[0], &v[1], &v[2]) == 3)
        {
            vertices.push_back(v);
        }

        if (sscanf(line.c_str(), "f %d %d %d", &e[0], &e[1], &e[2]) == 3)
        {
            for (int i = 0; i < 3; i++) e[i]--;
            triangles.push_back(e);
        }

        if (sscanf(line.c_str(), "vt %lg %lg", &u[0], &u[1]) == 2)
        {
            uvs.push_back(u);
        }

        if (sscanf(line.c_str(), "f %d/%d %d/%d %d/%d", &e[0], &t[0], &e[1], &t[1], &e[2], &t[2]) == 6)
        {
            for (int i = 0; i < 3; i++) e[i]--;
            triangles.push_back(e);
            for (int i = 0; i < 3; i++) t[i]--;
            triangle_texture_index.push_back(t);
        }
    }
    num_parts = triangles.size();
}

// Check for an intersection against the ray.
Hit Mesh::Intersection(const Ray& ray, int part) const
{
    Hit closest_hit;
    closest_hit.dist = std::numeric_limits<double>::infinity();

    if (part >= 0)
    {
        closest_hit = Intersect_Triangle(ray, part);
    }
    else
    {
        // Check all triangles
        for (int i = 0; i < (int)triangles.size(); i++)
        {
            Hit hit = Intersect_Triangle(ray, i);
            if (hit.dist >= small_t && hit.dist < closest_hit.dist)
            {
                closest_hit = hit;
            }
        }
    }

    if (closest_hit.dist == std::numeric_limits<double>::infinity())
    {
        closest_hit.dist = -1; // No intersection found
    }

    return closest_hit;
}

// Compute the normal direction for the triangle with index part.
vec3 Mesh::Normal(const Ray& ray, const Hit& hit) const
{
    assert(hit.triangle >= 0);

    ivec3 e = triangles[hit.triangle];
    vec3 A = vertices[e[0]];
    vec3 B = vertices[e[1]];
    vec3 C = vertices[e[2]];

    vec3 normal = cross(B-A, C-A).normalized();
    return normal;
}

Hit Mesh::Intersect_Triangle(const Ray& ray, int tri) const
{
    Hit hit;
    hit.triangle = -1;
    hit.dist = -1;

    // Retrieve the triangle vertices
    ivec3 e = triangles[tri];
    vec3 A = vertices[e[0]];
    vec3 B = vertices[e[1]];
    vec3 C = vertices[e[2]];

    // Compute the normal for the triangle
    vec3 normal = cross(B - A, C - A);
    double denominator = dot(normal, ray.direction);

    // Check if the ray is parallel to the triangle
    if (std::abs(denominator) < small_t) return hit;

    // Compute barycentric coordinates
    vec3 v = B - A;
    vec3 w = C - A;
    vec3 y = ray.endpoint - A;
    vec3 u = ray.direction;

    double beta = dot(cross(u, w), y) / dot(cross(u, w), v);
    double gamma = dot(cross(u, v), y) / dot(cross(u, v), w);
    double alpha = 1.0 - beta - gamma;
    double t = -dot(cross(v, w), y) / dot(cross(v, w), u);

    // Pixel_Print("mesh M triangle ", tri, " intersected; weights: (", alpha, " ", beta, " ", gamma, "); dist ", t);

    if (alpha >= -weight_tolerance && beta >= -weight_tolerance && gamma >= -weight_tolerance)
    {
        hit.dist = t;
        hit.triangle = tri;

        // Compute interpolated texture coordinates
        if (!triangle_texture_index.empty() && tri < (int)triangle_texture_index.size())
        {
            ivec3 tex_idx = triangle_texture_index[tri];
            vec2 uvA = uvs[tex_idx[0]];
            vec2 uvB = uvs[tex_idx[1]];
            vec2 uvC = uvs[tex_idx[2]];

            hit.uv = alpha * uvA + beta * uvB + gamma * uvC;
            // Pixel_Print("Computed UV: (", hit.uv[0], " ", hit.uv[1], ")");
        }
    }

    return hit;
}

std::pair<Box, bool> Mesh::Bounding_Box(int part) const
{
    if (part < 0)
    {
        Box box;
        box.Make_Empty();
        for (const auto& v : vertices)
            box.Include_Point(v);
        return {box, false};
    }

    ivec3 e = triangles[part];
    vec3 A = vertices[e[0]];
    Box b = {A, A};
    b.Include_Point(vertices[e[1]]);
    b.Include_Point(vertices[e[2]]);
    return {b, false};
}
