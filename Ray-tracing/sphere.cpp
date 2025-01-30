#include "sphere.h"
#include "ray.h"
#include <cmath>
#include <limits>

Sphere::Sphere(const Parse* parse, std::istream& in)
{
    in >> name >> center >> radius;

    // Ensure a valid radius
    if (radius <= 0)
    {
        throw std::runtime_error("Sphere radius must be greater than zero.");
    }
}

Hit Sphere::Intersection(const Ray& ray, int part) const
{
    vec3 oc = ray.endpoint - center;
    double a = dot(ray.direction, ray.direction);
    double b = 2 * dot(ray.direction, oc);
    double c = dot(oc, oc) - radius * radius;

    double discriminant = b * b - 4 * a * c;

    Hit hit;
    hit.dist = -1; // Initialize to invalid value
    hit.triangle = part; // For compatibility with mesh-based systems

    if (discriminant >= 0) // Valid intersection exists
    {
        double sqrt_discriminant = sqrt(discriminant);
        double t1 = (-b - sqrt_discriminant) / (2 * a);
        double t2 = (-b + sqrt_discriminant) / (2 * a);

        // Select the smallest positive t that is >= small_t
        if (t1 >= small_t && t2 >= small_t)
        {
            hit.dist = std::min(t1, t2);
        }
        else if (t1 >= small_t)
        {
            hit.dist = t1;
        }
        else if (t2 >= small_t)
        {
            hit.dist = t2;
        }
    }

    return hit; // Return a valid or invalid hit
}

vec3 Sphere::Normal(const Ray& ray, const Hit& hit) const
{
    // Ensure hit.dist is valid
    if (hit.dist < small_t)
    {
        throw std::runtime_error("Invalid hit distance in Sphere::Normal");
    }

    vec3 intersection_point = ray.Point(hit.dist);
    return (intersection_point - center).normalized();
}

std::pair<Box,bool> Sphere::Bounding_Box(int part) const
{
    return {{center-radius,center+radius},false};
}
