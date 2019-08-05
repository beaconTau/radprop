#! /usr/bin/env python3

# python port of point_to_point macro. Directly translated from a ROOT macro, with all its warts
# note that it assumes ../srtm-hd exists and has the right files in it 

import ROOT
import math
import array 
import sys

ROOT.gSystem.Load("../build/libradprop.so"); 


def point_to_point(lon1,lat1,h1,lon2,lat2,h2, freq=50, itwomv3=True, Npoints=10000,draw=True):

  #figure out the bounds we need 
  minlat = math.floor(min(lat1,lat2)); 
  maxlat = math.ceil(max(lat1,lat2)); 
  minlon = math.floor(min(lon1,lon2)); 
  maxlon = math.ceil(max(lon1,lon2)); 

  

  bounds = array.array('d', [minlon,minlat, maxlon, maxlat]) 

  dem = ROOT.radprop.DEM("../srtm-hd/",0,bounds);

  p1 = ROOT.radprop.SurfaceCoord(lon1,lat1);
  p2 = ROOT.radprop.SurfaceCoord(lon2,lat2);

  profile = ROOT.TGraph(Npoints); 
  dx = array.array('d',[1]); 
  dem.getHeightsBetween(Npoints, p1,p2, dx, profile.GetY(), False, profile.GetX()); 

  if draw:
    c = ROOT.TCanvas(); 
    c.Divide(2,1); 
    c.cd(1); 
    dem.getHist().Draw("col2z"); 
    gpts =ROOT.TGraph(2); 
    gpts.SetPoint(0, lon1,lat1); 
    gpts.SetPoint(1, lon2,lat2); 
    gpts.SetLineColor(2); 
    gpts.Draw("lsame"); 
    c.cd(2); 
    profile.SetFillColor(1); 
    profile.SetFillStyle(1); 
    profile.SetLineWidth(2); 
    profile.SetTitle("Profile; surface distance (m); altitude (m)"); 
    profile.Draw("AB1"); 



  r = ROOT.radprop.PropagationResult() 
  opt = ROOT.radprop.PropagationOptions() 
  opt.frequency = freq; 
  opt.method = ROOT.radprop.PropagationOptions.METHOD_ITWOMv3 if itwomv3 else  ROOT.radprop.PropagationOptions.METHOD_ITM 

  ROOT.radprop.propagate(Npoints, dx[0], profile.GetY(), r, h1, h2,opt); 

  print(" Propagation from  (%g,%g), TX height = %g to (%g,%g), RX height = %g at %g MHz" %(lon1,lat1,h1, lon2,lat2,h2, freq)) 
  print("    Path loss: %g" %( r.dBloss))
  print("    Mode: %s" %( r.mode)) 
  print("    Err: %d"%( r.err)); 



if __name__=="__main__": 
  lon1 = float(sys.argv[1])
  lat1 = float(sys.argv[2])
  h1 = float(sys.argv[3])
  lon2 = float(sys.argv[4])
  lat2 = float(sys.argv[5])
  h2 = float(sys.argv[6])
  freq = float(sys.argv[7]) if len (sys.argv) > 7 else 50
  point_to_point(lon1,lat1,h1,lon2,lat2,h2, freq,True,10000,False) 
  

