#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

extern bool enable_acceleration;

Render_World::~Render_World()
{
    for (auto a : all_objects) delete a;
    for (auto a : all_shaders) delete a;
    for (auto a : all_colors) delete a;
    for (auto a : lights) delete a;
}

// Find and return the Hit structure for the closest intersection. Ensure that hit.dist >= small_t.
std::pair<Shaded_Object, Hit> Render_World::Closest_Intersection(const Ray& ray) const
{
    Hit closest_hit;
    closest_hit.dist = std::numeric_limits<double>::infinity(); // Initialize to infinity
    Shaded_Object closest_object;

    // Iterate through all objects to find the closest intersection
    for (const auto& obj : objects)
    {
        Hit hit = obj.object->Intersection(ray, -1); // Check intersection
        if (hit.Valid() && hit.dist < closest_hit.dist && hit.dist >= small_t)
        {
            closest_hit = hit;
            closest_object = obj;
        }
    }

    return {closest_object, closest_hit};
}

// Set up the initial view ray and call Cast_Ray
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    Ray ray;
    ray.endpoint = camera.position; // Camera position as the ray origin
    ray.direction = (camera.World_Position(pixel_index) - camera.position).normalized(); // Direction toward the pixel

    vec3 color = Cast_Ray(ray, 1); // Cast ray with recursion depth = 1
    camera.Set_Pixel(pixel_index, Pixel_Color(color)); // Set the pixel color
}

void Render_World::Render()
{
    for (int j = 0; j < camera.number_pixels[1]; j++)
    {
        for (int i = 0; i < camera.number_pixels[0]; i++)
        {
            Render_Pixel(ivec2(i, j)); // Render each pixel
        }
    }
}

// Cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray, int recursion_depth) const
{
    if (recursion_depth > recursion_depth_limit)
        return vec3(0, 0, 0); // Return black if recursion depth exceeds the limit

    auto [closest_object, closest_hit] = Closest_Intersection(ray);

    if (closest_object.object)
    {
        // Calculate the intersection point and normal
        vec3 intersection_point = ray.Point(closest_hit.dist);
        vec3 normal = closest_object.object->Normal(ray, closest_hit);

        // Shade the surface using the object's shader
        return closest_object.shader->Shade_Surface(*this, ray, closest_hit, intersection_point, normal, recursion_depth);
    }
    else if (background_shader)
    {
        // Use the background shader if no object is intersected
        return background_shader->Shade_Surface(*this, ray, {}, {}, {}, recursion_depth);
    }

    return vec3(0, 0, 0); // Return black if no object and no background shader
}
