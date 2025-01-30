#include "acceleration.h"
#include "object.h"
#include "hit.h"

extern int acceleration_grid_size;

Acceleration::Acceleration()
{
    domain.Make_Empty();
    num_cells.fill(acceleration_grid_size);
}

void Acceleration::Add_Object(const Object* obj, int id)
{
    TODO;
}

ivec3 Acceleration::Cell_Index(const vec3& pt) const
{
    TODO;
    return {};
}

void Acceleration::Initialize()
{
    TODO;
}

std::pair<int,Hit> Acceleration::Closest_Intersection(const Ray& ray) const
{
    TODO;
    return {};
}
