
void pointToPoint(double lon1, double lat1, double h1, 
                  double lon2, double lat2, double h2,
                  double freq = 50, bool itwomv3 = true,  
                  int Npoints = 10000, bool draw = true) 
{

  //figure out the bounds we need 

  double minlat = floor(TMath::Min(lat1,lat2)); 
  double maxlat = ceil(TMath::Max(lat1,lat2)); 
  double minlon = floor(TMath::Min(lon1,lon2)); 
  double maxlon = ceil(TMath::Max(lon1,lon2)); 


  double bounds[4] = {minlon,minlat, maxlon, maxlat}; 

  radprop::DEM * dem = new radprop::DEM("srtm-hd/",0,bounds);

  radprop::SurfaceCoord p1(lon1,lat1);
  radprop::SurfaceCoord p2(lon2,lat2);

  TGraph * profile = new TGraph(Npoints); 
  double dx; 
  dem->getHeightsBetween(Npoints, p1,p2, &dx, profile->GetY(), false, profile->GetX()); 

  if (draw) 
  {
    TCanvas * c = new TCanvas; 
    c->Divide(2,1); 
    c->cd(1); 
    dem->getHist()->Draw("col2z"); 
    TGraph * gpts = new TGraph(2); 
    gpts->SetPoint(0, lon1,lat1); 
    gpts->SetPoint(1, lon2,lat2); 
    gpts->SetLineColor(2); 
    gpts->Draw("lsame"); 
    c->cd(2); 
    profile->SetFillColor(1); 
    profile->SetFillStyle(1); 
    profile->SetLineWidth(2); 
    profile->SetTitle("Profile; surface distance (m); altitude (m)"); 
    profile->Draw("AB1"); 
  }



  radprop::PropagationResult r; 
  radprop::PropagationOptions opt; 
  opt.frequency = 50; 
  opt.method = itwomv3 ? radprop::PropagationOptions::METHOD_ITWOMv3 :radprop::PropagationOptions::METHOD_ITM; 
 // opt.clutter_height = 1; 
//gg  opt.sea_level_refractivity = 350; 
  radprop::propagate(Npoints, dx, profile->GetY(), r, h1, h2,opt); 

  printf(" Propagation from  (%g,%g), TX height = %g to (%g,%g), RX height = %g at %g MHz\n",lon1,lat1,h1, lon2,lat2,h2, freq); 
  printf("    Path loss: %g\n", r.dBloss);; 
  printf("    Mode: %s\n", r.mode.c_str());; 
  printf("    Err: %d\n", r.err); 
}


