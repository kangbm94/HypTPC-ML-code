#include "TPCPadHelper.hh"
#ifndef ReconTools_h
#define ReconTools_h
const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = 0.90607;
const double& HS_field_Hall=  0.873800000;

double PI = acos(-1);
double mpi = 139.570/1000;//GeV
double mk = 493.677/1000;
double mp = 938.272/1000;
double mL = 1115.68/1000;
double mXi = 1321.71/1000;
double mXiStar = 1535/1000;


std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
TF1 fint("fint",s_tmp.c_str(),-4.,4.);


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
TVector3 VertexPointHelixLinear(const Double_t par1[5], const Double_t par2[5],
		Double_t& dist, Double_t& t1, Double_t& t2)
{
	//helix function 1
	//x = [0] + [3]*cos(t);
	//y = [1] + [3]*sin(t);
	//z = [2] + [3]*[4]*t;

	//helix function 2
	//x = [5] + [7]*t;
	//y = [6] + [8]*t;
	//z =  t;

	TF2 fvert_helix_lin("fvert_helix_lin",
			"pow(([0]+[3]*cos(x))-([5]+[7]*y),2)+pow(([1]+[3]*sin(x))-([6]+[8]*y),2)+pow(([2]+[3]*[4]*x)-y,2)",
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
	Double_t yout = par2[1]+par2[3]*close_zout;
	Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
	Double_t zout = close_zout;

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
	Double_t vertx = -1.*vx;
	Double_t verty = vz;
	Double_t vertz = vy + ZTarget;
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
class Track :public TLorentzVector{
	private:
		int Q_=0;
		int pid_=0;
		double cd_ = 1e9;
		double hpar[5]={0};
		int tid_ = -1;
		bool isHelix = true;
	public: 
		Track(int pid, int Q, double* par,int tid){
			Q_  = Q;
			pid_= pid;
			tid_= tid;
			for(int i=0;i<5;++i) hpar[i] = par[i];
		}

		bool IsP(){
			if(pid_%8/4) return true;
			else return false;
		}
		bool IsK(){
			if(pid_%4/2) return true;
			else return false;
		}
		bool IsPi(){
			if(pid_%2/1) return true;
			else return false;
		}
		bool GetQ(){
			return Q_;
		}
		int GetID(){
			return tid_;
		}
		void SetCD(double cd){
			cd_ = cd;
		}
		double GetCD(){
			return cd_;
		}
		double* GetPar(){
			return hpar;
		}
};
bool InTarget(TVector3 Vect){
	double x = Vect.X();
	double y = Vect.Y();
	double z = Vect.Z() - ZTarget;
	if(abs(x)<15 and abs(y) < 25 and abs(z)<10) return true;
	else return false;
}

class Recon{
	private:
		vector<TLorentzVector>Daughters;
		TLorentzVector LV;
		TVector3 Vert;
		bool exist = false;
		bool proper = false;
		int CombID=-1;
		double par[4];
	public:
		Recon(vector<TLorentzVector> D,TVector3 vertex,int id1,int id2){
			LV.SetXYZM(0,0,0,0);Daughters = D;Vert=vertex;
			for(auto lv:D) LV+=lv;
			CombID = pow(2,id1)+pow(2,id2);
			exist = true;
			auto V_t = GlobalToTarget(Vert);
			auto mom = LV.Vect();
			auto dir_t = GlobalToTargetMom(mom);
			dir_t = dir_t* (1/dir_t.Z());
			auto u = dir_t.X();
			auto v = dir_t.Y();
			par[0]=V_t.X()-V_t.Z()*u,par[1]=V_t.Y()-V_t.Z()*v,par[2]= u,par[3]=v;
		}
		Recon(){}
		void Clear(){
			exist = false;LV.SetXYZM(0,0,0,0);Daughters.clear();Vert.SetXYZ(0,0,0);
		}
		TVector3 Vertex(){
			return Vert;
		}
		double* GetPar(){
			return par;
		}
		bool Exist(){
			return exist;
		}
		TLorentzVector GetLV(){
			return LV;
		}
		double Mass(){
			return LV.Mag();
		}
		TVector3 Momentum(){
			return LV.Vect();
		}
		bool Counted(Track p){
			int trid = p.GetID();
			return (CombID%int(pow(2,trid+1)))/int(pow(2,trid));//true if Reconstructed with track.
		}
		int GetID(){
			return CombID;
		}
};


class Vertex{
	protected:
		vector<Track> Tracks;
		TVector3 vert;
		vector<TVector3> verts;
		double cdcut = 10;
		void SetVert(){
			int n = verts.size();
			vert=TVector3(0,0,0);
			for(int i=0;i<n;++i)vert+=verts[i] * (1./n);
		}
		int Vert_id=0;
		vector<Recon>LdCand;
		vector<TVector3>Pmom;
		vector<TVector3>Pimom;
		vector<TVector3>Ldmom;
		vector<Recon>XiCand;
	public:
		Vertex(Track p){
			Tracks.push_back(p);Vert_id=pow(2,p.GetID());
//			cout<<"Vertex"<<endl;
		}
		Vertex(){}
		bool Counted(Track p){
			int trid = p.GetID();
			return (Vert_id%int(pow(2,trid+1)))/int(pow(2,trid));//true if Reconstructed with track.
		}
		virtual int NTrack(){
			return Tracks.size();}
		virtual bool AddTrack(Track p);
		void SearchLdCombination();
		Recon GetLd(){
			double comp = 9999;
			Recon val ;
			int num = 0;
			for(auto ldc : LdCand){ if( abs(mL-ldc.Mass())<comp) {comp=abs(mL-ldc.Mass());val=ldc;}}
			return val;
		}
};
class VertexLH:public Vertex{
	private:
		vector<Recon>XiCand;
		vector<TVector3>Pimom;
		vector<TVector3>Ldmom;
		vector<Recon>Recons;
		vector<Track>Tracks;
	public:
		VertexLH(Recon p){
			Recons.push_back(p);
	}
		virtual int NTrack(){
			return Recons.size()+Tracks.size();
		}
		bool AddTrack(Track p);
		void SearchXiCombination();
		Recon GetXi(){
			double comp = 9999;
			Recon val ;
			int num = 0;
			for(auto xic : XiCand){ if( abs(mXi-xic.Mass())<comp) {comp=abs(mXi-xic.Mass());val=xic;}}
			return val;
		}
};
#endif
