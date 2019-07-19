#include "radprop/coord.h" 

#include <GeographicLib/PolarStereographic.hpp> 

void radprop::SurfaceCoord::convert(Mode mode) 
{

  if (m == mode) return; 

  if (m == WGS84)
  {
    double new_x, new_y; 
    GeographicLib::PolarStereographic::UPS().Forward(mode == NP_Stereographic, y,x, new_x, new_y); 
    m = mode; 
    x = new_x; 
    y = new_y; 
  }
  else
  {
    double new_x, new_y; 
    GeographicLib::PolarStereographic::UPS().Reverse(mode == NP_Stereographic, x,y, new_y, new_x); 
    m = mode; 
    x = new_x; 
    y = new_y; 
  }
}

