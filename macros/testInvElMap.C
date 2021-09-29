double lon =-118.2379382; 
double lat = 37.5898841;
double bounds[4] = {-121,35,-114,41};

void testInvElMap(double height = 3, double rmax = 300e3)
{

  radprop::SurfaceCoord x0(lon,lat); 
  radprop::DEM dem("srtm-hd/",0,bounds);

  TH2 * lon, *lat; 
  TH2 * h =  dem.makeInverseElevationAngleMap(x0, 360,-180,180,200,-10,10,height, rmax, -999,  &lon, &lat); 
  h->SetMinimum(0); 

  new TCanvas; 
  h->Draw("colz"); 

  new TCanvas; 
  lon->SetMinimum(bounds[0]);
  lon->Draw("colz"); 

  new TCanvas; 
  lat->SetMinimum(bounds[1]); 
  lat->Draw("colz"); 

}
