#ifndef _radprop_coord_h
#define _radprop_coord_h

namespace radprop
{

  struct SurfaceCoord
  {

    public: 

      enum Mode
      {
        WGS84, 
        SP_Stereographic, 
        NP_Stereographic
      };

      // x = lon or easting
      // y = lat or northing 
      SurfaceCoord(double x = 0, double y = 0, Mode mode = WGS84)  
      :  x(x), y(y), m(mode)
      {}

      static SurfaceCoord fromWGS84( double lat , double lon) { return SurfaceCoord(lon,lat,WGS84); } 
      void getWGS84( double & lat , double & lon) const { SurfaceCoord tmp = as(WGS84); lat = tmp.y; lon = tmp.x ; } 

      //convert to new mode
      void to(Mode mode) { if (mode == m) return; convert(mode); }

      // prefer to use above since it inlines the check that we're the same :) 
      void convert(Mode mode); 

      //return as new mode 
      SurfaceCoord as(Mode mode) const {  SurfaceCoord tmp(x,y,m); tmp.to(mode); return tmp; } 


      double x;
      double y; 
      Mode m; 
  }; 



}

#endif
