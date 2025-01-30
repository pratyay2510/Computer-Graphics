#ifndef __MISC_H__
#define __MISC_H__

#include "vec.h"
#include "ray.h"
#include <iostream>
#include <iomanip>
#include <sstream>

// Prints out a TODO message at most once.
#define TODO {static std::ostream& todo=std::cout<<"TODO: "<<__FUNCTION__<<" in "<<__FILE__<<std::endl;(void)todo;}

typedef unsigned int Pixel;

inline Pixel Pixel_Color(const vec3& color)
{
    unsigned int r=std::min(color[0],1.0)*255;
    unsigned int g=std::min(color[1],1.0)*255;
    unsigned int b=std::min(color[2],1.0)*255;
    return (r<<24)|(g<<16)|(b<<8)|0xff;
}

inline vec3 From_Pixel(Pixel color)
{
    return vec3(color>>24,(color>>16)&0xff,(color>>8)&0xff)/255.;
}

// Useful for creating indentation in pixel traces
struct Debug_Scope
{
    static bool enable;
    static int level;

    Debug_Scope(){level++;}
    ~Debug_Scope(){level--;}
};

// This routine is useful for generating pixel traces.  It only prints when the
// desired pixel is being traced.
// This routine is useful for generating pixel traces. It only prints when the
// desired pixel is being traced.
template<class... Args>
static void Pixel_Print(Args&&... args)
{
    if (!Debug_Scope::enable) return;
    for (int i = 0; i < Debug_Scope::level; i++) std::cout << "  ";
    (std::cout << ... << std::forward<Args>(args)) << std::endl;
}

// Macro for debugging at function entry
#define DEBUG_ENTER_FUNCTION(func_name)                       \
    Debug_Scope scope;                                        \
    if (Debug_Scope::enable)                                  \
        Pixel_Print("Entering function: ", func_name);

// Macro for printing variable values
#define DEBUG_VARIABLE(var_name, value)                      \
    if (Debug_Scope::enable)                                 \
        Pixel_Print(#var_name, " = ", value);

// Helper for printing vectors with formatting
inline std::string Vec_To_String(const vec3& v)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << "("
        << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return oss.str();
}

// Debugging example usage
inline void Debug_Ray(const std::string& label, const Ray& ray)
{
    Pixel_Print(label, " origin: ", Vec_To_String(ray.endpoint), 
                      ", direction: ", Vec_To_String(ray.direction));
}




inline int wrap(int i,int n)
{
    int k=i%n;
    if(k<0) k+=n;
    return k;
}

#endif

