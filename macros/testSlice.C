/* Example macro for beacon. 
 *
 * assumes relevant srtm-hd files are loaded in srtm-hd dir. 
 *
 * */ 


{

  gStyle->SetPalette(kRainBow); 
  double bounds[4] = {-119,36,-116,39};
  radprop::DEM dem("srtm-hd/",0,bounds);
  radprop::SurfaceCoord test(-118.2379382, 37.5898841);//approx beacon location 
  radprop::PropagationOptions opt; 
  opt.frequency = 50; //50 MHz 

  //make a 5 canvas thing 
  TCanvas * c = new TCanvas("c1","c1", 1800,1000); 
  c->Divide(1,5); 
  int ci = 1; 
  for (int bearing = 70; bearing <=110; bearing+=10)
  {
    c->cd(ci++); 
    radprop::VerticalSliceResult * r = new radprop::VerticalSliceResult;
    radprop::propagateVerticalSlice(*r,test,dem, 100e3, bearing, 2000, 200, 1, 1000, true,opt);

    r->pathloss.SetTitle(Form("bearing = %d",bearing)); 
    r->Draw(); 
  }
}
