#include "radprop/prop.h" 
#include "radprop/dem.h" 
#include "radprop/coord.h" 
#include "GeographicLib/Geodesic.hpp" 

#include "prop_private.h" 

const radprop::PropagationOptions & radprop::PropagationOptions::defaultPropagation()
{
  static PropagationOptions opt; 
  return opt; 
}



int radprop::propagate(int npts, double dx, const double *elev_in, 
                       PropagationResult & r, 
                       double tx_height, double rx_height, 
                       const PropagationOptions & opt)
{

  std::vector<double> elev(npts+2); 
  elev[0]=npts-1; 
  elev[1] = dx; 
  //this memcpy is annoying! TODO: change itwomv3.0.cpp to take a pointer instead 
  memcpy(&elev[0] +2, elev_in, npts*sizeof(double)); //blah 
 
  char strmode[128];  
  if (opt.method == PropagationOptions::METHOD_ITWOMv3) 
  {
    point_to_point(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                   opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization,
                   opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err, 
                   opt.clutter_refractivity, opt.clutter_height, opt.clutter_density, opt.mode_var, opt.delta_h_diff);
  }
  else 
  {
    point_to_point_ITM(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                       opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization,
                       opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err);
  }

  r.mode = strmode; 
  return r.err; 
}

int radprop::propagate(int npts, const SurfaceCoord &  start, const SurfaceCoord & stop, const DEM & dem, 
                       PropagationResult & r, 
                       double tx_height, double rx_height, 
                       const PropagationOptions & opt)
{
  //avoid a copy here by directly filling the array! 

  std::vector<double> elev(npts+2); 
  elev[0]=npts-1; 
  dem.getHeightsBetween(npts, start, stop, &elev[1], &elev[2]);
  
 
  char strmode[128];  
  if (opt.method == PropagationOptions::METHOD_ITWOMv3) 
  {
    point_to_point(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                   opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization, 
                   opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err,
                   opt.clutter_refractivity, opt.clutter_height, opt.clutter_density, opt.mode_var, opt.delta_h_diff);
  }
  else 
  {
    point_to_point_ITM(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                       opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization,
                       opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err);
  }

  r.mode = strmode; 
  return r.err; 
}



int radprop::propagateVerticalSlice (
   VerticalSliceResult & result, 
   const SurfaceCoord & fixed_pos, 
   const SurfaceCoord & max_variable_pos, 
   const DEM & dem, 
   int nxbins, 
   int nybins, 
   double fixed_height,
   double max_variable_above_fixed,
   bool variable_to_fixed, 
   const PropagationOptions & opt)
{

  result.terrain_profile.Set(nxbins); 

  result.terrain_profile.SetTitle("Terrain Profile; meters from fixed; elevation"); 
  double dx; 
  dem.getHeightsBetween(nxbins, fixed_pos, max_variable_pos, &dx, result.terrain_profile.GetY(), false, result.terrain_profile.GetX()); 
  
  double distance = (nxbins-1)*dx; 

  // if we are going from variable to fixed, we need to reverse the profile 
  std::vector<double> maybe_profile; 
  const double  * profile = result.terrain_profile.GetY();
  double fixed_alt = profile[0]; 
  double max_alt = fixed_alt + max_variable_above_fixed; 
  double min_alt = fixed_alt;
  if (variable_to_fixed) //we need to reverse it for that
  {
    maybe_profile.resize(nxbins); 
    for (int i = 0; i < nxbins; i++)
    {
      maybe_profile[i] = profile[nxbins-i-1]; 
      if (maybe_profile[i] < min_alt) min_alt = maybe_profile[i]; 
    }
    profile = &maybe_profile[0]; 
  }
  else //we still need to calculate the minimum
  {
    for (int i = 0; i < nxbins; i++) 
    {
      if (profile[i] < min_alt) min_alt = profile[i]; 
    }
  }

  result.pathloss.SetBins(nxbins,-dx/2,distance+dx/2, nybins, min_alt,max_alt); 

  float * arr = result.pathloss.GetArray(); 
  result.pathloss.SetEntries(nxbins*nybins); 

#ifdef ENABLE_OPENMP
#pragma omp parallel for 
#endif
  for (int y = 1; y <= nybins; y++) 
  {
    PropagationResult res; 
    double alt = result.pathloss.GetYaxis()->GetBinCenter(y); 
    for (int x =1 ; x <= nxbins; x++) 
    {
      int bin = variable_to_fixed ? y * (nxbins+2) + (nxbins+1-x) :  y * (nxbins+2) + x; 

      if (alt <profile[x-1]) 
      {
        arr[bin] = -1; 
      }
      else
      {
        propagate( variable_to_fixed ? nxbins-x+1 : x, dx, 
                   variable_to_fixed ? profile+x-1 : profile, 
                   res, 
                   variable_to_fixed ? alt - profile[x-1] : fixed_height, 
                   variable_to_fixed ? fixed_height : alt - profile[x-1], 
                   opt); 
        if (isnan(res.dBloss) || res.dBloss < 0) arr[bin] = -1; 
        else
        {
//          printf("%d %d %d: %g\n", x,y,bin,res.dBloss); 
          arr[bin] = res.dBloss; 
        }

      }

      

    }
  }

  result.pathloss.SetStats(false); 
  result.pathloss.SetTitle("Path loss; surface distance (m); altitude (m); path loss (dB)"); 
  result.terrain_profile.SetFillStyle(1); 
  result.terrain_profile.SetFillColor(1); 
  result.terrain_profile.SetLineWidth(2); 
  
  //todo, error handling
  return 0; 
}


int radprop::propagateVerticalSlice(VerticalSliceResult & result, 
                                    const SurfaceCoord & fixed_pos, 
                                    const DEM & dem, 
                                    double distance,
                                    double bearing, 
                                    int nxbins, 
                                    int nybins, 
                                    double fixed_height, double max_variable_above_fixed,
                                    bool variable_to_fixed, const PropagationOptions & opt)

{

  SurfaceCoord x0 = fixed_pos.as(SurfaceCoord::WGS84); 
  SurfaceCoord x1(0,0,SurfaceCoord::WGS84); 
  GeographicLib::Geodesic::WGS84().Direct( x0.y, x0.x, bearing, distance, x1.y, x1.x); 
  return propagateVerticalSlice(result, x0,x1, dem, nxbins, nybins, fixed_height, max_variable_above_fixed, variable_to_fixed, opt); 


}
