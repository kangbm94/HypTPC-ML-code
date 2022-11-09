#include "Math.hh"
#ifndef TPCGlobalFunctions_h
#define TPCGlobalFunctions_h
bool IsInsideTPC(double x, double y, double z){
	if(sqrt(x*x+z*z)<250) return true;
	else return false;
}
TGeoVolume* TPCGeometry() {

    //TPC Drawing
  Double_t flength = 586;
  Double_t fheight = 550;
  Double_t edge = flength/(1+sqrt(2));
  Double_t tan = TMath::Tan(22.5*TMath::Pi()/180.);
  
  TGeoManager *geom = new TGeoManager("K1.8HS","K1.8HS");
  TGeoMaterial *mat_world = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium *med_world = new TGeoMedium("med_world",1,mat_world);
  TGeoVolume *world = geom->MakeBox("world",med_world,flength*2,flength*2,flength*2);
  geom->SetTopVolume(world);
  world->SetVisibility(kFALSE);
  
  TGeoMaterial *mat_tpc = new TGeoMaterial("Al",26.98,13,2.7);
  TGeoMedium *med_tpc = new TGeoMedium("med_tpc",1,mat_tpc);
  TGeoXtru *xtru_tpc_inside = new TGeoXtru(2);
  Double_t zin[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		     tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t xin[] = { -flength/2, -edge/2, edge/2, flength/2,
		   flength/2, edge/2, -edge/2, -flength/2};  
  Double_t yin[] = { -fheight/2, fheight/2};
  Double_t scale[] = {1., 1.};
  Double_t x0[] = {0., 0.};
  Double_t z0[] = {0., 0.};

  Int_t nxz = sizeof(xin)/sizeof(Double_t);
  xtru_tpc_inside->DefinePolygon(nxz,xin,zin);
  
  Int_t n;
  Int_t ny = sizeof(yin)/sizeof(Double_t);
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_inside->DefineSection(n,yin[n],x0[n],z0[n],scale[n]);  

  TGeoVolume *tpc_inside = new TGeoVolume("tpc_inside",xtru_tpc_inside,med_tpc);

  ////////////////////////////////////////////////////////////////////////////
  flength = 586 - 10;
  edge = flength/(1+sqrt(2));

  TGeoXtru *xtru_tpc_outside = new TGeoXtru(2);
  Double_t z[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		     tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t x[] = { -flength/2, -edge/2, edge/2, flength/2,
		   flength/2, edge/2, -edge/2, -flength/2};  
  Double_t y[] = { -fheight/2, fheight/2};
  
  xtru_tpc_outside->DefinePolygon(nxz,x,z);  
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_outside->DefineSection(n,y[n],x0[n],z0[n],scale[n]);  
  
  TGeoVolume *tpc_outside = new TGeoVolume("tpc_outside",xtru_tpc_outside,med_tpc);

  //  new TGeoManager("Target","Tgt inside TPC");
  TGeoMaterial *mat_tgt = new TGeoMaterial("C",24,12,1.3);
  TGeoMedium *med_tgt = new TGeoMedium("med_tgt",1,mat_tgt);

  TGeoBBox *box_tgt = new TGeoBBox("box_tgt",30/2,20/2,20/2);
  TGeoVolume *tgt = new TGeoVolume("tgt",box_tgt,med_tgt);
  //  TGeoVolume *tgt = geom->MakeBox("tgt",med_tgt,20/2,30/2,20/2);

  //  tpc->SetLineColor(kBlack);
  //  tpc->SetLineWidth(1);
  TGeoVolumeAssembly *tpc = new TGeoVolumeAssembly("tpc");  
  tpc->AddNode(tgt,1,new TGeoTranslation(-143.,0.,0.));
  tpc->AddNode(tpc_outside,1);
  tpc->AddNode(tpc_inside,1);
  
  world->AddNode(tpc,1);

  geom->CloseGeometry();
  geom->SetVisLevel(4);

  return world;
}
const int max_nh=2500;
double Target_pos=-143,Target_x=34,Target_z=24;
double frame_width = 12.5;
enum{
	Else=0,
	L2PPi=1,
	L2NPi=2,
	KBeam=3
};
const int nbin=250;const int depth=1;
const double tpc_size=250;
short ToPixel(double x){
	x+=250;
	short x_pix = int(x* (double)nbin/tpc_size/2);
	return x_pix;
}
int ToShort(double y){
	y+=350;
	y*=10;
	return short(y);
}
bool IsInside(double z,double x){
	if(sqrt(z*z+x*x)<250){
		return true;
	}
	else{
		return false;
	}
}
bool IsTarget(double z, double x){
	if(abs(z-Target_pos)<Target_z/2&&abs(x)<Target_x/2){
		return true;
	}
	else{
		return false;
	}
}
bool IsActiveArea(double z, double x){
	double zr = z-Target_pos;
	if(abs(abs(z)-abs(x))<12.){
		return false;
	}
	else{
		return !IsTarget(z,x) and IsInside(z,x);
	}
}
static const int nhtpcmax = 300;

double StraightTrack(double x0,double u0,double z){
	return x0+u0*z;
}







#endif
