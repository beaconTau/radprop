/* Example macro for beacon. 
 *
 * assumes relevant srtm-hd files are loaded in srtm-hd dir. 
 *
 * */ 



void testSlice(bool vary_tx = false, double lat = 37.582232, double lon= -118.248259, double ant_height = 3)
{

  gStyle->SetPalette(kRainBow); 
  double bounds[4] = {-119,36,-116,39};
  radprop::DEM dem("srtm-hd/",0,bounds);
  radprop::SurfaceCoord test(lon,lat);//approx beacon location 
  radprop::PropagationOptions opt; 
  opt.frequency = 50; //50 MHz 
  opt.clutter_height = 0; 
  //opt.method = radprop::PropagationOptions::METHOD_ITM; 

  //make a 5 canvas thing 
  TCanvas * c = new TCanvas("c1","c1", 1800,1000); 
  c->Divide(1,5); 
  int ci = 1; 
  TFile f(Form("test_%s.root",vary_tx ? "vary_tx" : "vary_rx"), "RECREATE"); 

  for (int bearing = 70; bearing <=110; bearing+=10)
  {
    c->cd(ci++); 
    radprop::VerticalSliceResult * r = new radprop::VerticalSliceResult;
    radprop::propagateVerticalSlice(*r,test,dem, 100e3, bearing, 2001, 201, ant_height, 1000, vary_tx,opt);

    r->pathloss.SetTitle(Form("path loss bearing = %d",bearing)); 
    r->Draw(); 
//    c->cd(ci++); 
 //   r->mode.SetTitle(Form("propagation mode, bearing = %d",bearing)); 
 //   r->mode.Draw("col2z"); 
    r->Write(Form("bearing_%d", bearing)); 
  }

}
