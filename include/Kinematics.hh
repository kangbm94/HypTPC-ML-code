#include "TPCPadHelper.hh"
#ifndef Kinematics_h
#define Kinematics_h
const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = 0.90607;
const double& HS_field_Hall=  0.873800000;

std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
TF1 fint("fint",s_tmp.c_str(),-4.,4.);

TVector3 GlobalToTarget(TVector3 pos){
	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();
	double x_=-x;
	double y_=z-ZTarget; double z_=y; return TVector3(x_,y_,z_);
}
TVector3 TargetToGlobal(TVector3 pos){
	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();
	double x_=-x;
	double y_=z;
	double z_=y+ZTarget;
	return TVector3(x_,y_,z_);
}
TVector3 TargetToGlobalMom(TVector3 mom){
	double x = mom.X();
	double y = mom.Y();
	double z = mom.Z();
	double x_=-x;
	double y_=z;
	double z_=y;
	return TVector3(x_,y_,z_);
}
TVector3 GlobalToTargetMom(TVector3 mom){
	return TargetToGlobalMom(mom);
}

TVector3 VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
		Double_t& dist, Double_t& t1, Double_t& t2)
{
	//helix function 1
	//x = [0] + [3]*cos(t);
	//y = [1] + [3]*sin(t);
	//z = [2] + [3]*[4]*t;

	//helix function 2
	//x = [5] + [8]*cos(t);
	//y = [6] + [8]*sin(t);
	//z = [7] + [8]*[9]*t;

	TF2 fvert_helix("fvert_helix",
			"pow(([0]+[3]*cos(x))-([5]+[8]*cos(y)),2)+pow(([1]+[3]*sin(x))-([6]+[8]*sin(y)),2)+pow(([2]+[3]*[4]*x)-([7]+[8]*[9]*y),2)",
			-5.,5.,-5.,5.);

	fvert_helix.SetParameter(0, par1[0]);
	fvert_helix.SetParameter(1, par1[1]);
	fvert_helix.SetParameter(2, par1[2]);
	fvert_helix.SetParameter(3, par1[3]);
	fvert_helix.SetParameter(4, par1[4]);
	fvert_helix.SetParameter(5, par2[0]);
	fvert_helix.SetParameter(6, par2[1]);
	fvert_helix.SetParameter(7, par2[2]);
	fvert_helix.SetParameter(8, par2[3]);
	fvert_helix.SetParameter(9, par2[4]);

	Double_t close_zin, close_zout;
	fvert_helix.GetMinimumXY(close_zin, close_zout);
	t1 = close_zin;
	t2 = close_zout;
	dist = TMath::Sqrt(fvert_helix.GetMinimum());

	Double_t xin = par1[0]+par1[3]*cos(close_zin);
	Double_t xout = par2[0]+par2[3]*cos(close_zout);
	Double_t yin =  par1[1]+par1[3]*sin(close_zin);
	Double_t yout = par2[1]+par2[3]*sin(close_zout);
	Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
	Double_t zout =  par2[2]+par2[3]*par2[4]*close_zout;

	// Double_t vx = (par1[0]+par1[3]*cos(close_zin) + par2[0]+par2[3]*cos(close_zout))/2.;
	// Double_t vy = (par1[1]+par1[3]*sin(close_zin) + par2[1]+par2[3]*sin(close_zout))/2.;
	// Double_t vz = (par1[2]+par1[3]*par1[4]*close_zin + par2[2]+par2[3]*par2[4]*close_zout)/2.;
	Double_t vx = (xin+xout)/2.;
	Double_t vy = (yin+yout)/2.;
	Double_t vz = (zin+zout)/2.;

	Double_t dist2 = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
	dist = dist2;

	return   TargetToGlobal(TVector3(vx,vy,vz));

}
TVector3 VertexPointHelixLinear(const Double_t par1[5], const Double_t par2[4],
		Double_t& dist, Double_t& t1, Double_t& t2)
{
	//helix function 1
	//x = [0] + [3]*cos(t);
	//y = [1] + [3]*sin(t);
	//z = [2] + [3]*[4]*t;

	//helix function 2
	//x = [5] + [7]*t;
	//y = t;
	//z = [6]+[8]*t;

	TF2 fvert_helix_lin("fvert_helix_lin",
			"pow(([0]+[3]*cos(x))-([5]+[7]*y),2)+pow(([1]+[3]*sin(x))-y,2)+pow(([2]+[3]*[4]*x)-([6]+[8]*y),2)",
			-5.,5.,-250.,250.);
	fvert_helix_lin.SetParameter(0, par1[0]);
	fvert_helix_lin.SetParameter(1, par1[1]);
	fvert_helix_lin.SetParameter(2, par1[2]);
	fvert_helix_lin.SetParameter(3, par1[3]);
	fvert_helix_lin.SetParameter(4, par1[4]);
	fvert_helix_lin.SetParameter(5, par2[0]);
	fvert_helix_lin.SetParameter(6, par2[1]);
	fvert_helix_lin.SetParameter(7, par2[2]);
	fvert_helix_lin.SetParameter(8, par2[3]);

	Double_t close_zin, close_zout;
	fvert_helix_lin.GetMinimumXY(close_zin, close_zout);
	t1 = close_zin;
	t2 = close_zout;
	dist = TMath::Sqrt(fvert_helix_lin.GetMinimum());

	Double_t xin = par1[0]+par1[3]*cos(close_zin);
	Double_t xout = par2[0]+par2[2]*close_zout;
	Double_t yin =  par1[1]+par1[3]*sin(close_zin);
	Double_t yout = close_zout;
	Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
	Double_t zout = par2[1]+par2[3]*close_zout;

	// Double_t vx = (par1[0]+par1[3]*cos(close_zin) + par2[0]+par2[3]*cos(close_zout))/2.;
	// Double_t vy = (par1[1]+par1[3]*sin(close_zin) + par2[1]+par2[3]*sin(close_zout))/2.;
	// Double_t vz = (par1[2]+par1[3]*par1[4]*close_zin + par2[2]+par2[3]*par2[4]*close_zout)/2.;

	Double_t vx = (xin+xout)/2.;
	Double_t vy = (yin+yout)/2.;
	Double_t vz = (zin+zout)/2.;

/*
	double vx = xin;
	double vy = yin;
	double vz = zin;
	*/
	Double_t dist2 = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
	dist = dist2;
	return TargetToGlobal(TVector3(vx,vy,vz));
//	return TVector3(vertx, verty, vertz);
}


TVector3 CalcHelixMom(double par[5], double y)
{

	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

	double t = (y-par[2])/(par[3]*par[4]);
	double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c

	//From here!!!!
	double tmp_px = pt*(-1.*sin(t));
	double tmp_py = pt*(cos(t));
	double tmp_pz = pt*(par[4]);
	double px = -tmp_px*0.001;
	double py = tmp_pz*0.001;
	double pz = tmp_py*0.001;
	return TVector3(px,py,pz);
}



double GetTcal(double* par,TVector3 pos){
	TVector3 pos_(-pos.X(),
			pos.Z()-ZTarget,
			pos.Y());
	double fpar[8];
	for(int ip=0; ip<5; ++ip){
		fpar[ip] = par[ip];
	}
	fpar[5] = pos_.X();
	fpar[6] = pos_.Y();
	fpar[7] = pos_.Z();
	fint.SetParameters(fpar);
	double min_t = fint.GetMinimumX();
	return min_t;
}

TVector3 HelixPos(double* par, double t){
	double cx = par[0],cy=par[1],z0 = par[2],r = par[3],dz = par[4];
	double x = cx+r*cos(t);
	double y = cy+r*sin(t);
	double z = z0+r*dz*t;
	return TVector3(-x,z,y+ZTarget);
}

double MinHelixDistance(double* par, TVector3 pos){
	auto t = GetTcal(par,pos);
	auto v = HelixPos(par,t);
	return (pos-v).Mag();
}

vector<TVector3>gHitPos;
vector<TVector3>gRes;
static inline void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;   
  int dof = 0;
  
  double fpar[8];      
  //std::cout<<"paramter in fcn"<<std::endl;   
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  } 
 	int nh = gHitPos.size(); 
  for(int i=0; i<nh; ++i){ 
    TVector3 pos(-gHitPos[i].X(),
     gHitPos[i].Z()-ZTarget,
     gHitPos[i].Y());
    fpar[5] = pos.X(); 
    fpar[6] = pos.Y(); 
    fpar[7] = pos.Z();    
    fint.SetParameters(fpar);  
    double min_t = fint.GetMinimumX(); 
    double  x = par[0] + par[3]*cos(min_t);  
    double  y = par[1] + par[3]*sin(min_t);     
    double  z = par[2] + (par[4]*par[3]*min_t);        
    TVector3 fittmp(x, y, z);      
    TVector3 fittmp_(-1.*fittmp.X(),           
         fittmp.Z(),
         fittmp.Y()+ZTarget);      
    double tmp_t = atan2(pos.Y()-par[1], pos.X()-par[0]);  
    double  tmpx = par[0] + par[3]*cos(tmp_t);
    double  tmpy = par[1] + par[3]*sin(tmp_t); 
    double  tmpz = par[2] + (par[4]*par[3]*tmp_t); 
    TVector3 fittmp2_(-tmpx, tmpz,tmpy+ZTarget);
    TVector3 d = gHitPos[i] - fittmp_;
    TVector3 Res = gRes[i]; 
    double dxz = sqrt(d.x()*d.x()+d.z()*d.z());
    double Resxz = sqrt(Res.x()*Res.x()+Res.z()*Res.z());
    chisqr += pow(d.x()/Res.x(), 2) + pow(d.y()/Res.y(), 2) + pow(d.z()/Res.z(), 2);
    //    chisqr += pow(dxz/Resxz, 2) + pow(d.y()/Res.y(), 2);
    dof++;
    dof++;
  } 
  f = chisqr/(double)(dof-5);
}


#endif
