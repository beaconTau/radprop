#include "radprop/dem.h" 
#include "TSystem.h" 
#include "TSystemFile.h" 
#include "TSystemDirectory.h" 
#include <vector> 


#ifdef HAVE_GDAL
#include <gdal_priv.h> 
#endif 




static int dem_counter = 0;



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

  if (stereo && probably_stereo) 
  {
    if (stereo > 0 && strcasestr(projection,"northing"))
    {
      fprintf(stderr,"WARNING, this looks like it might be a South-Pole stereographic projection !\n"); 
    }
    else if (stereo < 0 && strcasestr(projection,"southing"))
    {
      fprintf(stderr,"WARNING, this looks like it might be a North-Pole stereographic projection !\n"); 
    }
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

    if (!fname.EndsWith(".hgt")) continue; 
    gSystem->GetPathInfo(Form("%s/%s", dir, fname.Data()), st); 


    if (st.fSize == 25934402 ) 
    {
      n_hd++; 
      hds.push_back(true); 
    }
    else if (st.fSize == 2884802 )
    {
      hds.push_back(false); 
      n_normal++; 
    }
    else
    {
      continue; 
    }

    char lon_sign; 
    char lat_sign;
    int ilon; 
    int ilat; 
    sscanf(fname.Data(),"%c%d%c%d.hgt", &lat_sign, &ilat, &lon_sign,&ilon); 

    double lon = ilon;
    double lat = ilat; 
    if (lon_sign == 'W') lon*=-1; 
    if (lat_sign == 'S') lat*=-1; 
    
    if (min_lat > lat) min_lat = lat;
    if (max_lat < lat+1) max_lat = lat+1; 
    if (min_lon > lon) min_lon = lon;
    if (max_lon < lon+1) max_lon = lon+1; 

    corners.push_back(std::pair<double,double>(lon,lat)); 
    files.push_back(fname); 
    if (verbose) 
      printf("Loaded %s SRTM file %s (%g,%g)\n", hds[hds.size()-1] ? "HD" : "Normal" , fname.Data(), lat, lon); 

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

    //otherwise, let's load it into memory

    TString this_file = Form("%s/%s",dir,files[i].Data()); 
    FILE * f = fopen(this_file.Data(),"rb"); 

    int N = hds[i]? 3601 : 1201; 
    int rd = fread(&data[0], 2, N*N, f); 
    if (rd != N*N) 
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
    the_hist.GetYaxis()->SetTitle(stereo < 0 ? "Northing" : "Southing"); 
  }
  else
  {
    the_hist.GetXaxis()->SetTitle("Longitude");
    the_hist.GetYaxis()->SetTitle("Latitude");
  }
  the_hist.SetStats(false); 

}




