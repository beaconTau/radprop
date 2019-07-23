#include "radprop/prop.h" 
#include "radprop/dem.h" 
#include "radprop/coord.h" 
#include "GeographicLib/Geodesic.hpp" 
#include "GeographicLib/GeodesicLine.hpp" 

#include "prop_private.h" 

ClassImp(radprop::VerticalSliceResult); 
ClassImp(radprop::HorizontalWedgeResult); 

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
  int imode;
  if (opt.method == PropagationOptions::METHOD_ITWOMv3) 
  {
    point_to_point(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                   opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization,
                   opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err,imode, 
                   opt.clutter_refractivity, opt.clutter_height, opt.clutter_density, opt.mode_var, opt.delta_h_diff);
  }
  else 
  {
    point_to_point_ITM(&elev[0], tx_height, rx_height, opt.surface_dielectric, opt.surface_conductivity, 
                       opt.sea_level_refractivity, opt.frequency, opt.radio_climate, opt.polarization,
                       opt.fraction_of_situations, opt.fraction_of_time, r.dBloss, strmode, r.err, imode);
  }

  r.mode = strmode; 
  r.imode = (radprop::PropagationResult::PropagationMode) imode; 
  return r.err; 
}

int radprop::propagate(int npts, const SurfaceCoord &  start, const SurfaceCoord & stop, const DEM & dem, 
                       PropagationResult & r, 
                       double tx_height, double rx_height, 
                       const PropagationOptions & opt)
{
  double dx; 
  std::vector<double> elev(npts); 
  dem.getHeightsBetween(npts, start, stop, &dx, &elev[0]);
  //this makes an extra copy, but otherwise we have to duplicate a lot of code... 
  return propagate(npts, dx, &elev[0], r, tx_height, rx_height, opt); 
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
   const PropagationOptions & opt, 
   SurfaceCoord * pts 
   )
{

  result.terrain_profile.Set(nxbins); 

  result.terrain_profile.SetTitle("Terrain Profile; meters from fixed; elevation"); 
  double dx; 
  dem.getHeightsBetween(nxbins, fixed_pos, max_variable_pos, &dx, result.terrain_profile.GetY(), false, result.terrain_profile.GetX(), pts); 
  
  double distance = (nxbins-1)*dx; 

  // if we are going from variable to fixed, we need to reverse the profile 
  std::vector<double> maybe_profile; 
  const double  * profile = result.terrain_profile.GetY();
  double fixed_alt = profile[0]; 
  double max_alt = fixed_alt + max_variable_above_fixed + fixed_height; 
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
  result.mode.SetBins(nxbins,-dx/2,distance+dx/2, nybins, min_alt,max_alt); 
  result.err.SetBins(nxbins,-dx/2,distance+dx/2, nybins, min_alt,max_alt); 

  float * arr = result.pathloss.GetArray(); 
  char * err_arr = result.err.GetArray(); 
  char * mode_arr = result.mode.GetArray(); 
  result.pathloss.SetEntries(nxbins*nybins); 
  result.mode.SetEntries(nxbins*nybins); 
  result.err.SetEntries(nxbins*nybins); 

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
        err_arr[bin] = -1; 
        mode_arr[bin] = -1; 
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
        mode_arr[bin] = (char) res.imode;
        err_arr[bin] = res.err;

      }

      

    }
  }

  result.pathloss.SetStats(false); 
  result.mode.SetStats(false); 
  result.err.SetStats(false); 
  result.pathloss.SetTitle("Path loss; surface distance (m); altitude (m); path loss (dB)"); 
  result.err.SetTitle("Path loss; surface distance (m); altitude (m); error code"); 
  result.mode.SetTitle("Path loss; surface distance (m); altitude (m); mode"); 
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
                                    bool variable_to_fixed, const PropagationOptions & opt, SurfaceCoord * pts)

{

  SurfaceCoord x0 = fixed_pos.as(SurfaceCoord::WGS84); 
  SurfaceCoord x1(0,0,SurfaceCoord::WGS84); 
  GeographicLib::Geodesic::WGS84().Direct( x0.y, x0.x, bearing, distance, x1.y, x1.x); 
  return propagateVerticalSlice(result, x0,x1, dem, nxbins, nybins, fixed_height, max_variable_above_fixed, variable_to_fixed, opt,pts); 


}


int radprop::propagateHorizontalWedge(HorizontalWedgeResult & result, const SurfaceCoord & fixed_pos, const DEM & dem, 
                                      double max_distance, int n_distance_points,  double min_bearing , double max_bearing , 
                                      double dbearing, 
                                      double fixed_height, double variable_height, 
                                      bool variable_relative_to_fixed, 
                                      bool variable_to_fixed, 
                                      const PropagationOptions & opt) 
{

  //first, let's get all the profiles 

  int nbearings = (max_bearing-min_bearing)/dbearing;
  std::vector<double> bearings(nbearings); 
  std::vector<std::vector<double> > elevation(nbearings); 
  std::vector<std::vector<std::pair<double,double>> > lonlat(nbearings); 

  SurfaceCoord x0 = fixed_pos.as(SurfaceCoord::WGS84); 
  double dx = max_distance / (n_distance_points-1); 

  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 

#ifdef ENABLE_OPENMP
#pragma omp parallel for 
#endif
  for (int i = 0; i < nbearings; i++) 
  {
    double bearing = min_bearing + dbearing*i; 
    bearings[i] = bearing; 
    elevation[i].resize(n_distance_points); 
    lonlat[i].resize(n_distance_points); 
    GeographicLib::GeodesicLine l(geod, x0.y,x0.x, bearing); 
    SurfaceCoord x(0,0,SurfaceCoord::WGS84); 
    for (int j = 0; j < n_distance_points; j++) 
    {
      l.Position(j*dx, x.y, x.x); 
      lonlat[i][j] = std::pair<double,double>(x.x,x.y); 
      elevation[i][j] = dem.getHeight(x);
    }
  }

  double fixed_elevation = elevation[0][0]; 

  ///now let's create the requisite TGraph2Ds
  //
  result.pathloss.SetTitle("Path loss (dB); longitude; latitude");
  result.terrain.SetTitle("Elevation (m); longitude; latitude");
  result.mode.SetTitle("Propagation Mode; longitude; latitude");
  result.err.SetTitle("Error; longitude; latitude");
  int tgraph2d_points = 1+nbearings * (n_distance_points-1); //avoid multiply counting the center point! 

  result.pathloss.Set(tgraph2d_points); 
  result.terrain.Set(tgraph2d_points);
  result.mode.Set(tgraph2d_points); 
  result.err.Set(tgraph2d_points); 

  result.pathloss.SetNpx(2*n_distance_points);
  result.pathloss.SetNpy(2*n_distance_points);
  result.terrain.SetNpx(2*n_distance_points);
  result.terrain.SetNpy(2*n_distance_points);
  result.mode.SetNpx(2*n_distance_points);
  result.mode.SetNpy(2*n_distance_points);
  result.err.SetNpx(2*n_distance_points);
  result.err.SetNpy(2*n_distance_points);


  
  int ii = 0;
  for (int i = 0; i < nbearings; i++) 
  {
    for (int j = 0; j < n_distance_points; j++) 
    {
      if ( j == 0 &&  i > 0) continue; //avoiid multiply counting the first point
      double lon = lonlat[i][j].first;
      double lat = lonlat[i][j].second;
      result.pathloss.SetPoint(ii, lon,lat,-1);
      result.terrain.SetPoint(ii, lon,lat,elevation[i][j]);
      result.mode.SetPoint(ii, lon,lat,-1);
      result.err.SetPoint(ii, lon,lat,-1);
      ii++; 
    }
    //now reverse the profile if we're doing variable to fixed
    if (variable_to_fixed) std::reverse(elevation[i].begin(), elevation[i].end());
    
  }


  double * loss = result.pathloss.GetZ();
  double * err = result.err.GetZ();
  double * mode = result.mode.GetZ();
  //alright, now actually start the propagation! 

#ifdef ENABLE_OPENMP
#pragma omp parallel for
#endif 
  for (int i = 0; i < nbearings; i++) 
  {
    PropagationResult res; 
    for (int j = 0; j < n_distance_points; j++) 
    {
      //skip either the first or last point for i > 0 to avoid multiply counting the origin
      if (j == (variable_to_fixed ? n_distance_points -1 : 0)  &&  i > 0) continue;

      //slightly confusing because of the extra point in the first iteration 
      int index =  i * (n_distance_points-1) + (i > 0)  + ( variable_to_fixed ?  n_distance_points - j - 1  : j);

      double var_height = variable_relative_to_fixed ?  fixed_height+variable_height + fixed_elevation - elevation[i][j] : variable_height; 
      if (var_height  < 0) continue; 

      double tx_height = variable_to_fixed ? var_height : fixed_height;
      double rx_height = variable_to_fixed ? fixed_height: var_height; 
      int npts = variable_to_fixed ? n_distance_points - j : j; 
      double * prf = variable_to_fixed ? &elevation[i][j] : &elevation[i][0]; 
      propagate(npts, dx, prf, res, tx_height, rx_height, opt); 

      if (isnan(res.dBloss) || res.dBloss < 0) loss[index] = -1; 
      else loss[index] = res.dBloss; 
      mode[index] = (double) res.imode;
      err[index] = res.err;
    }
  }

  
  return 0; 
}

