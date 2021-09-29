double lon =-118.2379382; 
double lat = 37.5898841;
double bounds[4] = {-120,36,-115,40};

TH2 * testElMap(double height = 3, double dr  = 100, double rmax = 300e3)
{

  radprop::SurfaceCoord x0(lon,lat); 
  radprop::DEM dem("srtm-hd/",0,bounds);

  return dem.makeElevationAngleMap(x0, 360,-180,180, rmax/dr, 0,rmax, height); 
}
