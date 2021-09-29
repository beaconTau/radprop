#ifndef _radprop_prop_h
#define _radprop_prop_h

#include "TH2.h" 
#include "TGraph2D.h" 
#include "TGraph.h" 
#include <string> 



/** The actual propagation utilities ! */ 

namespace radprop 
{

  class SurfaceCoord; 
  class DEM; 


  struct PropagationOptions
  {
    enum 
    {
      METHOD_ITWOMv3, 
      METHOD_ITM
    } method = METHOD_ITWOMv3; 

    enum 
    {
      POL_H, 
      POL_V,
      POL_CIRCULAR
    } polarization = POL_H; 

    enum 
    {
      CLIMATE_NONE, 
      CLIMATE_EQUATORIAL,
      CLIMATE_CONTINENTAL_SUBTROPICAL,
      CLIMATE_MARITIME_TROPICAL,
      CLIMATE_DESERT,
      CLIMATE_CONTINENTAL_TEMPERATE,
      CLIMATE_MARITIME_TEMPERATE_OVER_LAND,
      CLIMATE_MARITIME_TEMPERATE_OVER_SEA
    } radio_climate = CLIMATE_CONTINENTAL_TEMPERATE;  

    double surface_dielectric = 15;  // normal ground
    double surface_conductivity = 0.005; //
    double sea_level_refractivity = 301; 
    double frequency = 100; //MHz
    double fraction_of_situations = 0.5; 
    double fraction_of_time = 0.5; 

    double clutter_height = 22.5;  //only valid for itwomv3 
    double clutter_density = 1;    //only valid for itwomv3 
    double clutter_refractivity = 1e3; //only valid for itwomv3 
    int mode_var = 1; // not sure what this is, but only used by itwomv3. It says FCC =1, splat=12. Who knows. 
    double delta_h_diff = 0; //something above past-horizon propagation? only used in itwomv3

    //global instantiation of the default
    static const PropagationOptions & defaultPropagation(); 
  };


  
  struct PropagationResult 
  {
    double dBloss; 
    enum  PropagationMode
    {
      PROPAGATION_MODE_LOS, 
      SINGLE_HORIZON,
      SINGLE_HORIZON_DIFFRACTION_DOMINANT,
      SINGLE_HORIZON_TROPOSPHERIC_DOMINANT,
      SINGLE_HORIZON_PEAK, 
      DOUBLE_HORIZON, 
      DOUBLE_HORIZON_DIFFRACTION_DOMINANT,
      DOUBLE_HORIZON_TROPOSPHERIC_DOMINANT,
      DOUBLE_HORIZON_PEAK 
    } imode;
    std::string mode; 
    int err; 
  }; 

  /* Most basic path loss */ 

  
  int propagate(int npts, double dx, const double * elevation_profile,
                PropagationResult & result, 
                double tx_height =1, double rx_height=1, 
                const PropagationOptions & opt = PropagationOptions::defaultPropagation()); 


  /* Slightly simpler API */ 
  int propagate(int npts, const SurfaceCoord & tx, const SurfaceCoord & rx, const DEM & dem, 
                PropagationResult & result, 
                double tx_height =1, double rx_height=1, 
                const PropagationOptions & opt = PropagationOptions::defaultPropagation()); 



  // For doing a vertical slice in a given direction. 
  struct VerticalSliceResult : public TObject
  {
    TH2F pathloss; 
    TGraph terrain_profile; 
    TH2C mode; 
    TH2C err; 
    bool variable_to_fixed; 
    void Draw() 
    {
      pathloss.Draw("col2z"); 
      terrain_profile.Draw("b1 same"); 
    }

    virtual ~VerticalSliceResult() {;} 
    ClassDef(VerticalSliceResult,1); 
  }; 

  int  propagateVerticalSlice( VerticalSliceResult & result, 
                               const SurfaceCoord & fixed_pos, 
                               const SurfaceCoord & max_variable_pos, 
                               const DEM & dem,
                               int nxbins = 1000, 
                               int nybins = 100, 
                               double fixed_height =1, 
                               double max_variable_above_fixed = 1000,
                               bool variable_to_fixed = true, 
                               const PropagationOptions & opt = PropagationOptions::defaultPropagation(), 
                               SurfaceCoord * pts = 0 
                               );

  int propagateVerticalSlice( VerticalSliceResult & result, 
                              const SurfaceCoord & fixed_pos, 
                              const DEM & dem,
                              double distance = 100e3, 
                              double bearing = 0, 
                              int nxbins = 5001, 
                              int nybins = 201, 
                              double fixed_height =1, 
                              double max_variable_above_fixed = 1000, 
                              bool variable_to_fixed = true, 
                              const PropagationOptions & opt = PropagationOptions::defaultPropagation(),
                              SurfaceCoord * pts = 0 
                              ); 

  /*
   * Horizontal Wedge. 
   * Not really horizontal (since Earth curves). 
   *
   * note that the tx height can be relative to RX or relative to ground!! 
   *
   */

   struct HorizontalWedgeResult : public TObject
   {
     TGraph2D pathloss; 
     TGraph2D terrain; 
     TGraph2D mode; 
     TGraph2D err; 
     TGraph wedge_bounds; 
     double xmin,xmax,ymin,ymax; 
     bool tx_relative_to_rx; 
     ClassDef(HorizontalWedgeResult, 1); 
   };

  int propagateHorizontalWedge (HorizontalWedgeResult & result, 
                                const SurfaceCoord & fixed_pos, 
                                const DEM & dem, 
                                double max_distance = 100e3, 
                                int n_distance_points = 5001, 
                                double min_bearing = 0, 
                                double max_bearing = 360, 
                                double dbearing = 0.5 , 
                                double fixed_height = 1, 
                                double variable_height_= 0,
                                bool variable_relative_to_fixed = true, //otherwise, relative to ground! For SPLAT! like behavior, this should be false!
                                bool variable_to_fixed = false, 
                                const PropagationOptions & opt = PropagationOptions::defaultPropagation() 
                               );



  /** Makes a horizontal box around the RX instead 
   * of a wedge. Unlike the wedge, many more terrain profiles
   * must be generated. But the result is easier to use! 
   *
   * Of course, you can interpolate a 2D hist from the wedge instead... 
   */

  struct HorizontalBoxResult
  {
    TH2F pathloss; 
    TH2F terrain; 
  };

  /*
  int propagateHorizontalBox( HorizontalBoxResult & result, 
                             const SurfaceCoord & rx_pos, 
                             double north_distance = 100e3, 
                             double east_distance = 100e3, 
                             double south_distance = 100e3, 
                             double west_distance = 100e3, 
                             double dx = 1000, 
                             double rx_height = 1, 
                             double tx_height_= 0,
                             bool tx_relative_to_rx = true, //otherwise, relative to ground! For SPLAT! like behavior, this should be false!
                             const PropagationOptions & opt = PropagationOptions::defaultPropagation() ){return 1;} 
                             */


} 


#endif
