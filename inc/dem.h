#ifndef _radprop_dem_h 
#define _radprop_dem_h 

/** Digital Elevation Model class 
 * implemented as a TH2F. May be saved. 
 *
 * heights are in m above 
 *
 * See utilities for conversion. 
 ***/ 

#include "TH2.h" 

namespace GeographicLib
{
  class GeoCoords; 
}

namespace radprop
{



  class DEM 
  {
    public: 

      DEM() : stereo(0) {;}
      virtual ~DEM() {;} 

      /** Create a DEM from an input raster file.
       *  Requires GDAL for geotif files. 
       *
       *  If it's a directory, it is assumed to contain s3tm (or s3tm-hd) files and those may be read without gdal. 
       *  stereo should be zero if the file is in lat/lon. Otherwise it is assumed to be north (+1) or south (-1) pole steroegraphic. 
       *  If bounds is non-null, it is assumed to be a rectangle [xmin ymin xmax ymax] in whatever units the file has. 
       */
      DEM(const char * file, int stereo=0, 
          const double * bounds = NULL, bool verbose = false); 

      double getHeight(const GeographicLib::GeoCoords & where,  bool MSL=true ) const { return 0;}
      double * getHeightsBetween(int howmany,  const GeographicLib::GeoCoords & start, const GeographicLib::GeoCoords & stop, double * fill = 0,  bool MSL = true) const { return 0; }


      const TH2 * getHist() const { return &the_hist; } 
      TH2 * getHist() { return &the_hist; } 

    private: 
      TH2F the_hist; 
      int stereo; //0 for lon lat alt, -1 for south pole stereo, 1 for north pole stero


      ClassDef(DEM,1); 
  }; 



}




#endif

