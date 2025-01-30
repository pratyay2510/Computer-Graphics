#include "reflective_shader.h"
#include "parse.h"
#include "ray.h"
#include "render_world.h"

Reflective_Shader::Reflective_Shader(const Parse* parse, std::istream& in)
{
    in >> name;
    shader = parse->Get_Shader(in);
    in >> reflectivity;
    reflectivity = std::max(0.0, std::min(1.0, reflectivity));
}

vec3 Reflective_Shader::
Shade_Surface(const Render_World& render_world, const Ray& ray, const Hit& hit,
    const vec3& intersection_point, const vec3& normal, int recursion_depth) const
{
    // Base color from the underlying shader
    vec3 color = shader->Shade_Surface(render_world, ray, hit, intersection_point, normal, recursion_depth);

    // Define a small epsilon offset to avoid self-intersection
    const double epsilon = small_t;

    // Calculate reflection direction
    vec3 v_ray = ray.direction.normalized();
    vec3 r_dir = 2.0 * dot(-v_ray, normal) * normal + v_ray;
    Ray reflected_ray(intersection_point + normal*epsilon, r_dir);

    // Handle reflection contribution
    if (recursion_depth < render_world.recursion_depth_limit)
    {
        vec3 reflected_color = render_world.Cast_Ray(reflected_ray, recursion_depth + 1);
        color = (1 - reflectivity) * color + reflectivity * reflected_color;
    }
    else 
        color = (1-reflectivity)*color;

    return color;
}
