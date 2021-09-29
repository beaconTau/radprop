double lon =-118.23761867; 
double lat = 37.589339;
double bounds[4] = {-121,35,-114,41};

void testInvElMap(double height = 4, double rmax = 300e3)
{

  radprop::SurfaceCoord x0(lon,lat); 
  radprop::DEM dem("srtm-hd/",0,bounds);

  TFile f("beacon-site.root","RECREATE"); 
  TH2 * lon, *lat; 
  TH2 * h =  dem.makeInverseElevationAngleMap(x0, 3600,-180,180,400,-30,10,height, rmax, -999,  &lat, &lon); 
  h->SetMinimum(0); 
  h->SetStats(0); 

  new TCanvas; 
  h->Draw("colz"); 

  new TCanvas; 
  lon->SetMinimum(bounds[0]);
  lon->SetStats(0); 
  lon->Draw("colz"); 

  new TCanvas; 
  lat->SetStats(0); 
  lat->SetMinimum(bounds[1]); 
  lat->Draw("colz"); 

  h->Write(); 
  lon->Write(); 
  lat->Write(); 

}
