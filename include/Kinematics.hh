#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
#ifndef Kinematics_h
#define Kinematics_h
double MassSquare(double p, double l, double t){
	double beta = l / t / LightSpeed ;
	return p*p*(1-beta*beta)/beta/beta;
}


const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = 0.90607;
const double& HS_field_Hall=  0.873800000;
//const double& HS_field_Hall=  0.90820000;


std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
TF1 fint("fint",s_tmp.c_str(),-4.,7.);

TVector3 GlobalToTarget(TVector3 pos){
	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();
	double x_=-x;
	double y_=z-ZTarget;
	double z_=y;
	return TVector3(x_,y_,z_);
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
TLorentzVector ToTarget(TLorentzVector LV){
	double t = LV.T();
	double x = LV.X();
	double y = LV.Y();
	double z = LV.Z();
	double x_=-x;
	double y_=z;
	double z_=y;
	return TLorentzVector(x_,y_,z_,t);
}
TLorentzVector ToGlobal(TLorentzVector LV){
	return ToTarget(LV);
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
  double min_t = fint.GetMinimumX(-2.5*acos(-1),2.5*acos(-1)); 
	return min_t;
}
void GetHelixParameter(TVector3 Pos, TVector3 Mom,double charge,double* par){
	
	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

	auto Pos_helix = GlobalToTarget(Pos); 
	auto Mom_helix = GlobalToTargetMom(Mom); 
	double px = Mom_helix.X()*1000.;
	double py = Mom_helix.Y()*1000.;
	double pz = Mom_helix.Z()*1000.;

	double x = Pos_helix.X();
	double y = Pos_helix.Y();
	double z = Pos_helix.Z();

	double pt = sqrt(px*px+py*py);
	double t;
	if(charge>0){
		t = atan2(py,px)-0.5*acos(-1);
	}
	if(charge<0){
		t = atan2(py,px)+0.5*acos(-1);
	}
	double r = pt/(Const*dMagneticField);	
	double cx = x-r*cos(t);
	double cy = y-r*sin(t);
	double dz = pz / pt;
	double z0 =	z - r*dz*t;
	par[0]=cx;
	par[1]=cy;
	par[2]=z0;
	par[3]=r;
	par[4]=dz;
}

TVector3 VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
		Double_t& dist, Double_t& t1, Double_t& t2,double min_t1 = -4,double max_t1 = 7,double min_t2 = -4,double max_t2 = 7)
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
			min_t1,max_t1,min_t2,max_t2);

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

	Double_t xin = par1[0]+par1[3]*cos(t1);
	Double_t xout = par2[0]+par2[3]*cos(t2);
	Double_t yin =  par1[1]+par1[3]*sin(t1);
	Double_t yout = par2[1]+par2[3]*sin(t2);
	Double_t zin = par1[2]+par1[3]*par1[4]*t1;
	Double_t zout =  par2[2]+par2[3]*par2[4]*t2;

	Double_t vx = (xin+xout)/2.;
	Double_t vy = (yin+yout)/2.;
	Double_t vz = (zin+zout)/2.;

	Double_t distance = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
	dist = distance;

	return   TargetToGlobal(TVector3(vx,vy,vz));

}
TVector3 VertexPointHelixLinear(const Double_t par1[5], const Double_t par2[4],
		Double_t& dist, Double_t& t1, Double_t& t2)
{
	//helix function
	//{cx,cy,z0,r,dz}
	//x = [0] + [3]*cos(t); 
	//y = [1] + [3]*sin(t);  
	//z = [2] + [3]*[4]*t;

	//Linear function
	//{x0,z0,dxdy,dzdy} -> {-x0,y0,-u,v} in Global coordinate
	//x = [5] + [7]*t;
	//y = t;
	//z = [6]+[8]*t;

	TF2 fvert_helix_lin("fvert_helix_lin",
			"pow(([0]+[3]*cos(x))-([5]+[7]*y),2)+pow(([1]+[3]*sin(x))-y,2)+pow(([2]+[3]*[4]*x)-([6]+[8]*y),2)",
			-4.,7.,-100.,150.);
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
TVector3 CalcCircleMom(double par[5],TVector3 Pos )
{

	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

	auto HPos = GlobalToTarget(Pos);
	double x=HPos.x(),y=HPos.y();
	double x_c = x-par[0],y_c=y-par[1];
	double t = atan2(y_c,x_c);

//	double t = (y-par[2])/(par[3]*par[4]);
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

TVector3 HelixPos(double* par, double t){
	double cx = par[0],cy=par[1],z0 = par[2],r = par[3],dz = par[4];
	double x = cx+r*cos(t);
	double y = cy+r*sin(t);
	double z = z0+r*dz*t;
	return TVector3(-x,z,y+ZTarget);
}
TVector3 HelixPos(double* par, TVector3 pos){
	auto t = GetTcal(par,pos);
	auto v = HelixPos(par,t);
	return v;
}
double MinHelixDistance(double* par, TVector3 pos){
	auto t = GetTcal(par,pos);
	auto v = HelixPos(par,t);
	return (pos-v).Mag();
}

TVector3 HelixPosInTarget(double* par, double t){
	return GlobalToTarget(HelixPos(par,t));
}

double GetTcalWdist(double* par,TVector3 pos){
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
	double min_t=1e9 ;
	double distYMin = 1e9;
	for(int i = 0; i< 15; ++i){
		double min_t_temp = fint.GetMinimumX((-14.+2*i)*acos(-1),(-14.+2*(i+1))*acos(-1)); 
		double distY = (HelixPosInTarget(fpar,min_t_temp)-pos).Mag();
		if(distY<distYMin){
			distYMin = distY;
			min_t = min_t_temp;
		}
	}
	// min_t = fint.GetMinimumX(0.*acos(-1),2.*acos(-1)); 
	return min_t;
}

double GetTcalBeam(double* par,TVector3 pos){
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
  double min_t = fint.GetMinimumX(0.5*acos(-1),1.5*acos(-1)); 
	return min_t;
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
TVector3
VertexPoint(const TVector3& Xin, const TVector3& Xout,
            const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;
  Double_t x = 0.5*(x1+x2);
  Double_t y = 0.5*(y1+y2);
  if(std::isnan(x) || std::isnan(y) || std::isnan(z))
    return TVector3(NAN, NAN, NAN);

  return TVector3(x, y, z);

}
Double_t
CloseDist(const TVector3& Xin, const TVector3& Xout,
          const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;

  return TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
double GetTrackHitAngle(double* par, TVector3 pos){
	auto mom = CalcCircleMom(par,pos);
	double px = mom.x(),pz=mom.z();
	double norm_p = sqrt(px*px+pz*pz);
	double x = pos.x(),z=pos.z()-ZTarget;
	double norm_x = sqrt(x*x+z*z);
	double prd= (px*x+pz*z)/norm_p/norm_x;	
	if(prd<0)prd*=-1;
	return acos(prd);
}

/*
double VertexFunction(double* xyz,double* res, double* param){
	TVector3 pos(xyz[0],xyz[1],xyz[2]);	
	TVector3 h_pos = GlobalToTarget(pos);
	double hz = h_pos.z();
	h_pos.z()=0;
	TVector3 cent = TVector3(param[0],param[1],0);
	h_pos -= cent;
	double dir = h_pos * (1./h_pos.Mag());
	auto rad = param[2]*dir;
	auto dw = rad - h_pos;
	double t = atan2(h_pos.x(),-h_pos.y());
	double dw = cos(t)*dw.x() + sin(t)*dw.y();
	double dl = cos(t)*dw.y() - sin(t)*dw.x();

}

*/
static double MomToRad(double p, double v){
	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
	double pt  = p * 1./sqrt(1+v*v); 
	double rad = p * 1000. / (Const*dMagneticField); 	
	return rad;
}
static double RadToMom(double rad, double v){
	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
	double pt  = rad  * (Const * dMagneticField);
	double p = pt * hypot(1,v);
	return p/1000.; 

}


#endif
