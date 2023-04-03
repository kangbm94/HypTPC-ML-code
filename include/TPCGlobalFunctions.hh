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
Double_t HypTPCdEdx(Double_t Z, Double_t *x, Double_t *p){
  //x : poq
  //p[0] : converting constant p[1] : density effect correction p[2] : mass
  Double_t me  = 0.5109989461;
  Double_t rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672); //[g cm-3]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t ZoverA = 17.2/37.6; //[mol g-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
  Double_t I2 = 0.9*188.0 + 0.1*41.7; I2 = I2*I2; //Mean excitaion energy [eV]
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t mom = 1000.*x[0]*Z; //MeV
  Double_t beta2 = mom*mom/(mom*mom+p[2]*p[2]);
  Double_t gamma2 = 1./(1.-beta2);
  Double_t Wmax = 2*me*beta2*gamma2/((me+p[2])*(me+p[2])+2*me*p[2]*(TMath::Sqrt(gamma2)-1));
  Double_t dedx = p[0]*constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2)-beta2-p[1]);
  return dedx;
}

Double_t HypTPCBethe(Double_t *x, Double_t *p){ return HypTPCdEdx(1, x, p); }
Int_t HypTPCdEdxPID_temp(Double_t dedx, Double_t poq){
  Double_t bethe_par[2] = {7195.92, -10.5616};
  Double_t limit = 0.6; //GeV/c
  Double_t mpi = 139.57039;
  Double_t mk  = 493.677;
  Double_t mp  = 938.2720813;
  Double_t md  = 1875.612762;
  TF1 *f_pim = new TF1("f_pim", HypTPCBethe, -3., 0., 3);
  TF1 *f_km = new TF1("f_km", HypTPCBethe, -3., 0., 3);
  TF1 *f_pip = new TF1("f_pip", HypTPCBethe, 0., 3., 3);
  TF1 *f_kp = new TF1("f_kp", HypTPCBethe, 0., 3., 3);
  TF1 *f_p = new TF1("f_p", HypTPCBethe, 0., 3., 3);
  TF1 *f_d = new TF1("f_d", HypTPCBethe, 0., 3., 3);

  f_pim -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_km -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_pip -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_kp -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_p -> SetParameters(bethe_par[0], bethe_par[1], mp);
  f_d -> SetParameters(bethe_par[0], bethe_par[1], md);

  Int_t pid[3] = {0};
  if(poq >= limit){
    pid[0]=1; pid[1]=1; pid[2]=1;
  }
  else if(limit > poq && poq >= 0.){
    Double_t dedx_d = f_d -> Eval(poq); Double_t dedx_p = f_p -> Eval(poq);
    Double_t dedx_kp = f_kp -> Eval(poq); Double_t dedx_pip = f_pip -> Eval(poq);
    if(dedx_d > dedx && dedx >= dedx_kp) pid[2]=1;
    if(dedx_p > dedx){
      pid[0]=1; pid[1]=1;
    }
  }
  else if(0.> poq && poq >= -limit){
    pid[0]=1; pid[1]=1;
  }
  else{
    pid[0]=1; pid[1]=1;
  }

  delete f_pim;
  delete f_km;
  delete f_pip;
  delete f_kp;
  delete f_p;
  delete f_d;

  Int_t output = pid[0] + pid[1]*2 + pid[2]*4;
  return output;
}







#endif
