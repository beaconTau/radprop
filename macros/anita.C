/* Example macro for anita. 
 *
 * assumes rema file is present
 *
 * */ 



radprop::VerticalSliceResult * anita(bool vary_tx = true)
{

  gStyle->SetPalette(kRainBow); 

  radprop::SurfaceCoord anita(-323694.84, -82432.419, radprop::SurfaceCoord::SP_Stereographic);//anita for event 721X
  double anita_height=38575.090; 
  double bearing=341.59-203; 
  double around = 1200e3; 
  double bounds[4] = {anita.x-around,anita.y-around,anita.x+around,anita.y+around};
  radprop::DEM dem("REMA_200m_dem_filled.tif",-1,bounds);
  double fixed_height = anita_height-dem.getHeight(anita); 
  printf("anita height above ground: %g\n", fixed_height); 
  radprop::PropagationOptions opt; 
  opt.frequency = 300; //50 MHz 
  opt.clutter_height = 1; 
  opt.surface_dielectric = 2.25; 
  opt.surface_conductivity = 1e-3; 

  //opt.method = radprop::PropagationOptions::METHOD_ITM; 
  //
  //

  int npts = 10001; 
  std::vector<radprop::SurfaceCoord> path(npts); 


  TCanvas * c = new TCanvas("c1","c1", 1800,1000); 
  radprop::VerticalSliceResult * r = new radprop::VerticalSliceResult;
  radprop::propagateVerticalSlice(*r,anita,dem, 1000e3, bearing, npts, 1501, fixed_height, 1000, vary_tx,opt,&path[0]);
  r->Draw(); 

  new TCanvas; 

  dem.getHist()->DrawCopy("col2z"); 

  TGraph * g = new TGraph(npts); 
  for (int i = 0; i < npts; i++)
  {
    path[i].to(radprop::SurfaceCoord::SP_Stereographic); 
    g->SetPoint(i, path[i].x, path[i].y); 
  }
  TMarker * m = new TMarker(anita.x,anita.y,23); 
  g->Draw("lsame"); 
  m->Draw(); 

  return r;

}
