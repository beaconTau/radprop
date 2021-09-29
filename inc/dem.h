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
#include "radprop/coord.h" 


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

      /** Returns the height at the given position (interpolating the DEM 
       * Returns -999 if out of range. 
       *
       * */ 
      double getHeight(const SurfaceCoord & where,  bool MSL=false ) const;
      

      class Path
      {
        public: 

          Path(const SurfaceCoord & start, const SurfaceCoord & end);
          Path(const SurfaceCoord & start, double az, double distance); 
          double getDistance() const { return distance; } 
          double getAzimuth() const { return az; } 
          double getReverseAzimuth() const { return az_back; } 

          //guaranteed to be in WGS84 
          const SurfaceCoord & getStart() const { return start; } 
          const SurfaceCoord & getEnd() const { return end; } 

        private: 
          SurfaceCoord start; 
          SurfaceCoord end; 
          double az; 
          double az_back; 
          double distance; 
      }; 

      /** Retrieves howmany heights between two points. If fill is passed, it
       * will be filled (and returned) otherwise a new array is allocated. if X
       * is passed, the distances are from the start are filled there (useful
       * for plotting). If pts is passed, they are filled with the points.        */ 

      enum HeightMode
      {
        Height_WGS84, 
        Height_MSL, 
        Height_Relative, 
        Height_RelWithCurvature
      };

      double * getHeightsBetween(int howmany,  const Path & path, double * dx = 0, 
                                  double * fill = 0,  HeightMode  mode = Height_WGS84, double * X =0,   SurfaceCoord  * pts = 0) const; 


      TH2 * makeElevationAngleMap(const SurfaceCoord & where, int naz, double az_min, double az_max, double nr, double rmin, double rmax, 
          double height_off_ground = 0, double noval = -999, 
          bool tolatlon = false, double dlat = 0.01, double dlon = 0.01) const; 

      TH2 * makeInverseElevationAngleMap(const SurfaceCoord & where, int naz, double az_min, double az_max, int nel, double elmin, double elmax,
          double height_off_ground = 0, double rmax = 300e3, double noval = -999, 
           TH2** lat_map = 0, TH2** lon_map = 0) const; 


      const TH2 * getHist() const { return &the_hist; } 
      TH2 * getHist() { return &the_hist; } 

    private: 
      TH2F the_hist; 
      int stereo; //0 for lon lat alt, -1 for south pole stereo, 1 for north pole stero


      ClassDef(DEM,1); 
  }; 



}




#endif

