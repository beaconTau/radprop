#include "radprop/dem.h" 
#include "TSystem.h" 
#include "TSystemFile.h" 
#include "TSystemDirectory.h" 
#include <vector> 
#include "TMath.h" 
#include <zlib.h> 


#ifdef HAVE_GDAL
#include <gdal_priv.h> 
#endif 

#include <GeographicLib/Geoid.hpp> 
#include <GeographicLib/GeodesicLine.hpp> 
#include <GeographicLib/Geodesic.hpp> 
#include <GeographicLib/Ellipsoid.hpp> 

static int dem_counter = 0;


static GeographicLib::Geoid geoid("egm96-5","",true,true); 


#ifdef HAVE_GDAL
static int setupFromGDAL(const char * file, TH2F & h, const double * bounds, int stereo, bool verbose)
{

  GDALDataset * ds; 
  static bool registered = false; 

  if (!registered) 
  {
    GDALAllRegister(); 
    registered = true; 
  }


  ds = (GDALDataset*) GDALOpen(file, GA_ReadOnly); 
  if (!ds) return 1; 

  int nx = ds->GetRasterXSize(); 
  int ny = ds->GetRasterYSize(); 
  int nlevels = ds->GetRasterCount(); //for

  if (nlevels > 1) 
  {
    fprintf(stderr,"WARNING, for now only read first raster layer...\n"); 
  }

  const char * projection = ds->GetProjectionRef(); 

  bool probably_stereo = strcasestr(projection,"stereographic"); 
  if (!stereo && probably_stereo) 
  {
    fprintf(stderr,"WARNING, this looks like it's a stereographic projection!\n"); 
  }
  else if (stereo && !probably_stereo) 
  {
    fprintf(stderr,"WARNING, this looks like it might not be a stereographic projection!\n"); 
  }

  if (verbose && projection) printf("Projection: %s\n", projection); 

  double affineTransform[6]; 
  ds->GetGeoTransform(affineTransform); 

  double xmin = affineTransform[0]; 
  double ymax = affineTransform[3]; //top, not bottom 
  double dx = affineTransform[1];
  double dy = fabs(affineTransform[5]);
  //TODO check no other affine transforms 
  
  double xmax = dx *nx + xmin;
  double ymin = -dy*ny + ymax; 

//  printf("x:[%g %g] y:[%g %g]\n",xmin,xmax,ymin,ymax); 

  int imin = 0;
  int imax = nx; 
  int jmin = 0;
  int jmax = ny; 

  if (bounds) 
  {
    //the order is important here since we're being stingy about the variables 
    if (bounds[2] < xmax) 
    {
      imax = 1+(bounds[2] - xmin) /dx;
      xmax = bounds[2]; 
    }

    if (bounds[3] < ymax) 
    {
      jmax = 1+(bounds[3] - ymin) /dy;
      ymax = bounds[3]; 
    }

    if (bounds[0] > xmin) 
    {
      imin = (bounds[0]-xmin)/dx;
      xmin = bounds[0]; 
    }
    if (bounds[1] > ymin) 
    {
      jmin = (bounds[1]-ymin)/dy; 
      ymin = bounds[1]; 
    }
  }

  int hist_nx = (imax-imin);
  int hist_ny = (jmax-jmin);

  h.SetBins( hist_nx, xmin, xmax, hist_ny, ymin, ymax); 

  //now we loop over rows 
  GDALRasterBand * band = ds->GetRasterBand(1); 
  h.SetEntries(hist_nx*hist_ny); 
  for (int j = jmin; j < jmax; j++) 
  {
    //figure out our offset inside the histogram
    float * arr = h.GetArray()  +  (1+(j-jmin)) * (hist_nx+2) + 1 ; 
    band->RasterIO(GF_Read, 
                   imin, ny-j-1, hist_nx, 1, //read one line 
                   arr, hist_nx, 1, 
                   GDT_Float32,0,0); 
    
    //set any crazy values to 0 for now (usually this is the sea)
    for (int i = 0; i < hist_nx; i++) 
    {
      if (arr[i] < -1000) arr[i]=0; 
    }
  }

  delete ds; 

  return 0;
}
#endif 

static int setupFromSRTM(const char * dir, TH2 & h, const double * bounds, bool verbose)
{

  //let's figure out what files we have 
 
  TSystemDirectory d(dir,dir); 
  TList * l_files = d.GetListOfFiles(); 

  TIter next(l_files); 
  
  int n_hd = 0; 
  int n_normal = 0; 
  std::vector<TString> files; 
  std::vector<std::pair<double,double> > corners; 
  std::vector<bool> hds; 

  TSystemFile * f= 0;
  double min_lat = 999;
  double max_lat = -999; 
  double min_lon = 999;
  double max_lon = -999; 
  FileStat_t st;; 

  while ((f = (TSystemFile*) next()))
  {

    TString fname = f->GetName(); 

    if (!fname.EndsWith(".hgt") && !fname.EndsWith(".hgt.gz")) continue; 
    gSystem->GetPathInfo(Form("%s/%s", dir, fname.Data()), st); 

    bool gzipped = fname.EndsWith("gz"); 


    if (verbose) printf("Considering %s/%s\n", dir, fname.Data()); 
    
    char lon_sign; 
    char lat_sign;
    int ilon; 
    int ilat; 
    sscanf(fname.Data(),gzipped ? "%c%d%c%d.hgt.gz" : "%c%d%c%d.hgt", &lat_sign, &ilat, &lon_sign,&ilon); 

    double lon = ilon;
    double lat = ilat; 
    if (lon_sign == 'W') lon*=-1; 
    if (lat_sign == 'S') lat*=-1; 

    if (bounds && lat+1 < bounds[1]) continue; 
    if (bounds && lat >  bounds[3]) continue; 
    if (bounds && lon+1 < bounds[0]) continue; 
    if (bounds && lon >  bounds[2]) continue; 
 
    if (verbose) printf("  in bounds!\n"); 

    int size = 0; 

    if (gzipped) 
    {
      //read the last four bytes and hope! 
      FILE * f = fopen(Form("%s/%s",dir,fname.Data()),"r"); 
      fseek(f, st.fSize-4,SEEK_SET); 
      fread(&size,4,1,f); 
      fclose(f); 
    }
    else
    {
      size = st.fSize; 
    }


    if (size == 25934402 ) 
    {
      n_hd++; 
      if (verbose) printf(" is HD!\n"); 
      hds.push_back(true); 
    }
    else if (size == 2884802 )
    {
      hds.push_back(false); 
      if (verbose) printf(" is SD!\n"); 
      n_normal++; 
    }
    else
    {
      if (verbose) printf(" Wrong size (%d)\n", size); 
      continue; 
    }

   
    if (min_lat > lat) min_lat = lat;
    if (max_lat < lat+1) max_lat = lat+1; 
    if (min_lon > lon) min_lon = lon;
    if (max_lon < lon+1) max_lon = lon+1; 

    corners.push_back(std::pair<double,double>(lon,lat)); 
    files.push_back(fname); 
    if (verbose) 
      printf("Adding %s SRTM file %s (%g,%g) to list \n", hds[hds.size()-1] ? "HD" : "Normal" , fname.Data(), lat, lon); 

  }

  bool hd = n_hd; 
  if (n_hd && n_normal) 
  {
    fprintf(stderr,"WARNING: found a mix of HD and normal files... will use normal resolution\n"); 
    hd = false; 
  }

  if (!n_hd && !n_normal) 
  {
    fprintf(stderr,"WARNING: found no valid files!\n"); 
    return 1; 
  }


  double xmin = min_lon; 
  double xmax = max_lon; 
  double ymin = min_lat; 
  double ymax = max_lat; 
  
  if (bounds) 
  {
    if (bounds[0] > xmin) xmin = bounds[0]; 
    if (bounds[1] > ymin) ymin = bounds[1]; 
    if (bounds[2] < xmax) xmax = bounds[2]; 
    if (bounds[3] < ymax) ymax = bounds[3]; 
  }

  double resolution = hd? (1./3600) : (3./3600); 

  int nx = (xmax-xmin)/resolution+1; 
  int ny = (ymax-ymin)/resolution+1; 
  
  h.SetBins( nx, xmin, xmax, ny, ymin, ymax); 

  //alright, now let's do the loop over files 

  std::vector<uint16_t> data(n_hd ? 3601*3601 : 1201*1201); 
  for (unsigned i = 0; i < files.size(); i++) 
  {

    //check if in bounds
    if (corners[i].first+1 < xmin) continue; 
    if (corners[i].second+1 < ymin) continue; 
    if (corners[i].first > xmax) continue; 
    if (corners[i].second > ymax) continue; 


    if (verbose) 
       printf("USING  %s SRTM file %s\n", hds[hds.size()-1] ? "HD" : "Normal" , files[i].Data()); 

    //otherwise, let's load it into memory

    TString this_file = Form("%s/%s",dir,files[i].Data()); 
    bool gzipped = this_file.EndsWith(".gz"); 

    int N = hds[i]? 3601 : 1201; 
    int rd = 0;
    if (gzipped) 
    {
      gzFile f= gzopen(this_file.Data(),"r"); 
      rd = gzread(f,&data[0],2*N*N); 
      gzclose(f); 
    }
    else
    {
      FILE * f = fopen(this_file.Data(),"rb"); 
      rd = fread(&data[0],1 , 2*N*N, f); 
      fclose(f); 
    }
    if (rd != 2*N*N) 
    {
      fprintf(stderr,"Hmm, couldn't read all of %s\n", this_file.Data()); 
      return 1; 
    }
    
    int min_x_bin = 1 + (corners[i].first - xmin)/resolution; 
    int min_y_bin = 1 + (corners[i].second - ymin)/resolution; 


    for (int ii = 0; ii < N;  ii++) 
    {
      int ibin = min_x_bin + ii; 
      if (ibin < 1) continue; 
      if (ibin > nx) continue; 
      for ( int jj = 0; jj < N; jj++) 
      {
        int jbin = min_y_bin + jj; 
        if (jbin < 1) continue; 
        if (jbin > ny) continue; 

        int data_bin = !hd && hds[i] ? 1+3*ii+(1+3*(N-jj-1)*N) :  ii+(N-jj-1)*N; 
        double val = short ((data[data_bin] & 0xff)*256 + (data[data_bin] >> 8)); 
        if (val > -1000) 
          h.SetBinContent(ibin,jbin,val); 

      }
    }
  }
  return 0; 
}


radprop::DEM::DEM(const char * file, 
                  int stereo, const double * bounds, bool verbose) 
  :stereo(stereo) 
{

  FileStat_t st; 
  gSystem->GetPathInfo(file,st); 


  the_hist.SetName(Form("dem_%03d", dem_counter++)); 

  int ok = 0; 
  //we have a whole bunch of srtm files! 
  if (R_ISDIR(st.fMode)) 
  {
    ok = setupFromSRTM(file, the_hist, bounds, verbose); 
  }

  else
  {
#ifdef HAVE_GDAL
   ok = setupFromGDAL(file,the_hist,bounds, stereo,verbose); 
#else
   fprintf(stderr,"Not compiled with GDAL Support!!!\n\n"); 
   ok = 1
#endif
  }


  if (ok) 
  {
    the_hist.SetTitle("FAILED"); 
    return; 
  }

  the_hist.SetTitle(file); 
  the_hist.SetMinimum(-100); 

  if (stereo) 
  {
    the_hist.GetXaxis()->SetTitle("Easting"); 
    the_hist.GetYaxis()->SetTitle("Northing"); 
  }
  else
  {
    the_hist.GetXaxis()->SetTitle("Longitude");
    the_hist.GetYaxis()->SetTitle("Latitude");
  }
  the_hist.SetStats(false); 

}



double radprop::DEM::getHeight(const SurfaceCoord & where, bool msl) const
{
  double x,y; 



  if (!stereo)
  {
    SurfaceCoord wgs84 = where.as(SurfaceCoord::WGS84); 
    x = wgs84.x;//lon 
    y = wgs84.y;//lat 
  }
  else
  {
    SurfaceCoord s = where.as(stereo < 0 ? SurfaceCoord::SP_Stereographic : SurfaceCoord::NP_Stereographic); 
    x = s.x;
    y = s.y;
  }


  //why the hell is interpolate not const? 
  TH2 * h = (TH2*) (&the_hist); 

  //if out of bounds, assume 0? 
  if (x > h->GetXaxis()->GetXmax() || x < h->GetXaxis()->GetXmin() || y > h->GetYaxis()->GetXmax() || y < h->GetYaxis()->GetXmin() )
  {
    if (msl) return 0; 


    SurfaceCoord wgs84 = where.as(SurfaceCoord::WGS84); 
    double alt = geoid.ConvertHeight(wgs84.y, wgs84.x, 0, GeographicLib::Geoid::ELLIPSOIDTOGEOID); 

    return alt; 

  }

  

  double alt = h->Interpolate(x,y); 

  if (msl) 
  {
    SurfaceCoord wgs84 = where.as(SurfaceCoord::WGS84); 
    //hopefully this is the right way! 
    alt = geoid.ConvertHeight(wgs84.y, wgs84.x, alt, GeographicLib::Geoid::ELLIPSOIDTOGEOID); 
  }

  return alt; 
}


double * radprop::DEM::getHeightsBetween(int howmany,  const Path & path, double * fill_dx, 
                                  double * fill,  HeightMode mode , double * X , SurfaceCoord * pts) const 
{

  double * H = fill ?: new double[howmany]; 


  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 

  //now make the geodesic line
  GeographicLib::GeodesicLine  l(geod, path.getStart().y, path.getStart().x, path.getAzimuth()); 


  double Rc = 0; 
  if (mode == Height_RelWithCurvature) 
  {
    Rc = GeographicLib::Ellipsoid::WGS84().NormalCurvatureRadius(path.getStart().y, path.getAzimuth()); 
  }

  double dx = path.getDistance()/(howmany-1); 
  if (fill_dx) *fill_dx = dx; 


  double H0 = 0; 
  for (int i = 0; i < howmany; i++) 
  {
    double lat, lon; 
    double x = i*dx; 
    l.Position(x, lat,lon); 

    H[i] = getHeight( SurfaceCoord(lon,lat, SurfaceCoord::WGS84), mode == Height_MSL); 
    if (i == 0) H0 = H[0]; 

    if (X) 
    {
      X[i] = x; 
    }

    if (mode == Height_RelWithCurvature) 
    {
      
      double dh = Rc/cos(x/Rc)-Rc; 
      H[i] -= dh; 
      H[i]-=H0; // make it all relative to the start
      if (X) //X in this case is horizontal distance
      {
        X[i] = sqrt(2*Rc*dh+dh*dh); 
      }
    }
    else if (mode == Height_Relative) 
    {
      H[i]-=H0; // make it all relative to the start
    }


    if (pts) 
    {
      pts[i] = SurfaceCoord(lon,lat, SurfaceCoord::WGS84); 
    }
  }

  return H; 
}




radprop::DEM::Path::Path(const SurfaceCoord & start, const SurfaceCoord & end) 
  : start(start.as(SurfaceCoord::WGS84)), end(end.as(SurfaceCoord::WGS84)) 
{

  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 
  geod.Inverse(start.y,start.x, end.y,end.x, distance,az, az_back); 
}


radprop::DEM::Path::Path(const SurfaceCoord & start, double azi, double distance) 
  : start(start.as(SurfaceCoord::WGS84)), az(azi), distance(distance) 
{
  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 
  geod.Direct(start.y,start.x,azi, distance, end.y, end.x, az_back ); 
}


TH2* radprop::DEM::makeElevationAngleMap(const SurfaceCoord & where, int naz, double az_min, double az_max, double nr, double rmin, double rmax, double height, double noval, bool latlon, double dlat, double dlon) const
{

  SurfaceCoord x0 = where.as(SurfaceCoord::WGS84); 
  TH2 * M = new TH2F("elevmap",Form("Elevation Angle Map around %g,%g\n ; Azimuth; Distance [m]; Elevation Angle ", x0.y,x0.x), naz, az_min, az_max, nr, rmin, rmax); 
  M->SetDirectory(0); 

  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 
  double daz = (az_max-az_min)/naz; 
  double dr = (rmax-rmin)/nr; 

  int nrsegs = rmax/dr; 

  std::vector<double> hs(nrsegs+1); 
  std::vector<double> ls(nrsegs+1); 



  double lat_min = x0.y;
  double lat_max = x0.y;
  double lon_min = x0.x; 
  double lon_max = x0.x; 

  for (int i = 0; i < naz; i++) 
  {
    double az = az_min + (i+0.5) * daz; 
    Path path(where,az,rmax); 

    if (path.getEnd().x < lon_min) lon_min = path.getEnd().x; 
    if (path.getEnd().y < lat_min) lat_min = path.getEnd().y; 
    if (path.getEnd().x > lon_max) lon_max = path.getEnd().x; 
    if (path.getEnd().y > lat_max) lat_max = path.getEnd().y; 

    getHeightsBetween(nrsegs+1,path, 0,&hs[0], Height_RelWithCurvature,&ls[0]); 

    for (int j = 0; j < nr; j++) 
    {
      int jj = j-(rmin/dr); 
      double h = 0.5*(hs[jj]+hs[jj+1]); 
      double l = 0.5*(ls[jj]+ls[jj+1]); 
      double el = atan2(h-height, l); 

      //ok now see if this is occluded; 

      double y = height; 
      double x = 0; 
      bool ok = true; 
      double dx = 10; // 10 meter step? 
      double dy = dx * sin(el); 
      int seg = 0; 
      double l0 = ls[0]; 
      double l1 = ls[1]; 
      double dl = l1-l0; 
      while (x < l)
      {
        if ( x > l1 ) 
        {
          l0 = l1;
          seg++; 
          l1 = ls[1+seg]; 
          dl = l1-l0; 
        }

        double frac = (x-l0) / dl; 
        double yterrain = (1-frac)*hs[seg]  + frac*hs[seg+1]; 
        if ( y < yterrain) 
        {
          ok = false; 
          break; 
        }

        x+=dx; 
        y+=dy; 
      }

      if (!ok) 
      {
        M->SetBinContent(i+1,j+1, noval); 
      }
      else
      {
        M->SetBinContent(i+1,j+1, TMath::RadToDeg() * el); 
      }

    }
  }

  if (!latlon) return M; 

  double nlat = (lat_max-lat_min)/dlat; 
  double nlon = (lon_max-lon_min)/dlon; 

  TH2 * M2 = new TH2F("elevmaplatlon", M->GetTitle(), nlon, lon_min,lon_max, nlat, lat_min, lat_max); 
  M2->SetDirectory(0); 
  M2->GetXaxis()->SetTitle("Longitude"); 
  M2->GetYaxis()->SetTitle("Latitude"); 


  for (int i = 1; i <=nlon; i++) 
  {
    double lon = M2->GetXaxis()->GetBinCenter(i); 
    for (int j = 1; j <=nlat; j++) 
    {
      double lat = M2->GetYaxis()->GetBinCenter(j); 
      double r, az; 
      geod.Inverse(x0.y, x0.x, lat, lon, r, az); 
      if (r < rmax && r > rmin && az < az_max && az > az_min) 
      {
        M2->SetBinContent(i,j,M->Interpolate(az,r)); 
      }
      else
      {
        M2->SetBinContent(i,j,noval); 
      }
    }
  }

  delete M; 
  return M2; 
}


TH2* radprop::DEM::makeInverseElevationAngleMap(const SurfaceCoord & where, int naz, double az_min, double az_max, 
                                                int nel, double elmin, double elmax, double height, 
                                                double max_r, double noval, TH2 ** lat_map, TH2 ** lon_map) const
{

  SurfaceCoord x0 = where.as(SurfaceCoord::WGS84); 
  TH2 * M = new TH2F("invelevmap",Form("Distance Map around %g,%g\n ; Azimuth [deg]; Elevation Angle [deg]; Surface Distance [m] ", x0.y,x0.x), naz, az_min, az_max, nel, elmin, elmax); 
  M->SetDirectory(0); 

  if (lat_map) 
  {
    *lat_map = new TH2F("invelevmaplat",Form("Latitude Map around %g,%g\n ; Azimuth [deg]; Elevation Angle [deg]; Latitude", x0.y,x0.x), naz, az_min, az_max, nel, elmin, elmax); 
    (*lat_map)->SetDirectory(0); 
  }

  if (lon_map) 
  {

    *lon_map = new TH2F("invelevmaplon",Form("Longitude Map around %g,%g\n ; Azimuth [deg]; Elevation Angle [deg]; Longitude ", x0.y,x0.x), naz, az_min, az_max, nel, elmin, elmax); 
    (*lon_map)->SetDirectory(0); 
  }


  const GeographicLib::Geodesic & geod = GeographicLib::Geodesic::WGS84(); 
  double daz = (az_max-az_min)/naz; 
  double del = (elmax-elmin)/nel; 

  int nrsegs = max_r * 0.01;  // 10 m segments? 

  std::vector<double> hs(nrsegs+1); 
  std::vector<double> ls(nrsegs+1); 


  double lonlat_noval = noval < -180 && noval > 180 ? noval : -999; 

  for (int i = 0; i < naz; i++) 
  {
    double az = az_min + (i+0.5) * daz; 
    Path path(where,az,max_r); 

    double dr; 
    getHeightsBetween(nrsegs+1,path, &dr,&hs[0], Height_RelWithCurvature,&ls[0]); 

    for (int j = 0; j < nel;  j++) 
    {
      double el = TMath::DegToRad() * (elmin + (j+0.5) * del); 

      //ok now see if this is occluded; 

      double y = height; 
      double x = 0; 
      bool ok = false; 
      double dx = 10; // 10 meter step? 
      double dy = dx * sin(el); 
      int seg = 0; 
      double l0 = ls[0]; 
      double l1 = ls[1]; 
      double dl = l1-l0; 
      double r = noval; 
      while (x < ls[nrsegs])
      {
        if ( x > l1 ) 
        {
          l0 = l1;
          seg++; 
          l1 = ls[1+seg]; 
          dl = l1-l0; 
        }

        double frac = (x-l0) / dl; 
        double yterrain = (1-frac)*hs[seg]  + frac*hs[seg+1]; 
        if ( y < yterrain) 
        {
          r = dr * (seg+frac); 
          ok = true; 

          break; 
        }

        x+=dx; 
        y+=dy; 
      }

      M->SetBinContent(i+1,j+1, r); 

      if ( lat_map || lon_map) 
      {
        double lat = lonlat_noval, lon = lonlat_noval; 
        if (ok) 
        {
          geod.Direct(x0.y,x0.x, az, r, lat,lon); 
        }

        if (lat_map) 
        {
          (*lat_map)->SetBinContent(i+1, j+1, lat);  
        }

        if (lon_map) 
        {
          (*lon_map)->SetBinContent(i+1, j+1, lon);  
        }
      }

    }
  }



  return M;
}
