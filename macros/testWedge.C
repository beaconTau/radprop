/* Example macro for beacon. 
 *
 * assumes relevant srtm-hd files are loaded in srtm-hd dir. 
 *
 * */ 



radprop::HorizontalWedgeResult * testWedge(bool vary_tx = false, double lon =-118.2379382, double lat = 37.5898841, double fixed_height = 3, double var_height = 1)
{

  gStyle->SetPalette(kRainBow); 
  double bounds[4] = {-120,36,-115,40};
  radprop::DEM dem("srtm-hd/",0,bounds);
  radprop::SurfaceCoord test(lon,lat);//approx beacon location 
  radprop::PropagationOptions opt; 
  opt.frequency = 50; //50 MHz 
  //opt.method = radprop::PropagationOptions::METHOD_ITM; 

  //make a 5 canvas thing 
  TCanvas * c = new TCanvas("c1","c1", 1800,1000); 
  c->Divide(2,1); 
  TFile f(Form("wedge_%s.root",vary_tx ? "vary_tx" : "vary_rx"), "RECREATE"); 


  radprop::HorizontalWedgeResult * r = new radprop::HorizontalWedgeResult;
  radprop::propagateHorizontalWedge(*r,test,dem, 100e3, 

                                     5001, 0,360,0.5, fixed_height, var_height,false,vary_tx,opt); 

  c->cd(1); 
  r->terrain.Draw("colz"); 
  c->cd(2); 
  gStyle->SetPalette(kRainBow);
  r->pathloss.Draw("colz"); 
  return r;
}
