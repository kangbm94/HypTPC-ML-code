#include "Kinematics.hh"
#include "KinFit.hh"
#ifndef ReconTools_h
#define ReconTools_h

double PI = acos(-1);
double mpi = 139.570/1000;//GeV
double mpi0 = 134.976/1000;//GeV
double mk = 493.677/1000;
double mp = 938.272/1000;
double mL = 1115.677/1000;
double mXi = 1321.71/1000;
double mXiStar = 1535./1000;
double mH = 2150./1000;

class Track :public TLorentzVector{
	private:
		int Q_=0;
		int pid_=0;
		double cd_ = 1e9;
		double hpar[5]={0};
		int tid_ = -1;
		bool isHelix = true;
	public: 
		Track(){
		}
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
		void SetP(){
			pid_ = 4;
		}
		void SetK(){
			pid_ = 2;
		}
		void SetPi(){
			pid_ = 1;
		}
		void SetQ(int charge){
			Q_=charge;
		}

		int GetQ(){
			return Q_;
		}
		bool IsNegative(){
			if(Q_ > -1) return false;
			else return true;
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
	protected:
		vector<TLorentzVector>Daughters;
		TLorentzVector LV;
		TVector3 Vert;
		bool exist = false;
		bool proper = false;
		int CombID=-1;
		double par[5];
		double close_dist=-1;
		int trid1=-1,trid2=-1;
		Track ReconTrack;
		int charge;
	public:
		Recon(vector<TLorentzVector> D,TVector3 vertex,double clos_dist,int id1,int id2,int charge_=0);
		Recon(Recon P,Recon Q,double m1,double m2);
		Recon(){}
		void Clear(){
			exist = false;LV.SetXYZM(0,0,0,0);Daughters.clear();Vert.SetXYZ(0,0,0);
		}
		TVector3 Vertex(){
			return Vert;
		}
		void SetVertex(TVector3 v){
			Vert = v;
		}
		double* GetPar(){
				return par;
		}
		void SetExistance(bool flag){
			exist = flag;
		}
		int GetCharge(){
			return charge;
		}
		bool Exist(){
			return exist;
		}
		TLorentzVector GetLV(){
			return LV;
		}
		void SetLV(TLorentzVector LV_){
			LV = LV_;
		}
		double Mass(){
			return LV.Mag();
		}
		double GetCD(){
			return close_dist;
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
		int GetID1(){
			return trid1;
		}
		int GetID2(){
			return trid2;
		}
		TLorentzVector GetDaughter(int i){
			return Daughters.at(i);
		}
		void SetTrack(Track T){
			ReconTrack = T;
		}
		Track GetTrack(){
			return ReconTrack;
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
		vector<Recon>XiCandLdCor;
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
		void SetCdCut(double cd){
			cdcut = cd;
		}
};
class VertexLH:public Vertex{
	private:
		vector<Recon>XiCand;
		vector<Recon>XiCorCand;
		vector<TVector3>Pimom;
		vector<TVector3>Ldmom;
		vector<Recon>Recons;
		vector<Track>Tracks;
		KinematicFitter Fitter;
		double laglambda = 0;
		vector<Track> KuramaTracks  ;
		bool KuramaFlag = false;
	public:
		VertexLH(Recon p, bool KinematicFit = true){
			if(KinematicFit and p.Exist()){
				
				Fitter = KinematicFitter(mL);
				auto PLV = p.GetDaughter(0);
				auto PP = PLV.Vect();
				auto PPres = PP * (0.03);
				auto PiLV = p.GetDaughter(1);
				auto PiP = PiLV.Vect();
				auto PiPres = PiP * (0.03);
				Fitter.AssignLorentzVector(PLV,PiLV);
				Fitter.SetResolution(PPres,PiPres);
				Fitter.DoKinematicFit();
				auto PLVCor = Fitter.GetFittedLV().at(0);
				auto PiLVCor = Fitter.GetFittedLV().at(1);
				auto LVCor = PLVCor + PiLVCor;
				p.SetLV(LVCor);
				laglambda = Fitter.GetLambda();
			}
			else if (p.Exist()){
				auto LdLV = p.GetLV();
				auto LdP = LdLV.Vect();
				double x = LdP.X(),y=LdP.Y(),z=LdP.Z();
				auto LdLVCor = TLorentzVector();
				LdLVCor.SetXYZM(x,y,z,mL);
//				p.SetLV(LdLVCor);
				p.SetLV(LdLV);
			}
			Recons.push_back(p);
		}
		virtual int NTrack(){
			return Recons.size()+Tracks.size();
		}
		bool AddTrack(Track p);
		void AddKuramaTrack(Track p);
		void SearchXiCombination();
		double GetLagMulti(){
			return laglambda;
		}
		Recon GetLd(){
			return Recons[0];
		}
		Recon GetXi(){
			double comp = 9999;
			Recon val ;
			int num = 0;
			for(auto xic : XiCand){ if( abs(mXi-xic.Mass())<comp) {comp=abs(mXi-xic.Mass());val=xic;}}
			return val;
		}
		Recon GetXiCor(){
			double comp = 9999;
			Recon val ;
			int num = 0;
			for(auto xic : XiCorCand){ if( abs(mXi-xic.Mass())<comp) {comp=abs(mXi-xic.Mass());val=xic;}}
			return val;
		}
};
class XiStarRecon: public Recon{
	private:
		vector<TVector3>KMmom;
		vector<TVector3>KMPos;	
		vector<TVector3>KPmom;
		vector<TVector3>KPPos;	
		TVector3 TargetPos = TVector3(0,0,-143);
		double mass = mXiStar;
		TVector3 IniMom;
		int kpid = 0;
		int kmid = 0;
		Track K18Track;
		Track KuramaTrack;
	public:
		XiStarRecon(){}
		void Construct(vector<TVector3>KMX,vector<TVector3>KMP,vector<TVector3>KPX,vector<TVector3>KPP);
		TVector3 GetIniMom(){
			return IniMom;
		}
		void SetK18Track(Track T){
			K18Track = T;
		}
		void SetKuramaTrack(Track T){
			KuramaTrack = T;
		}
		TVector3 GetKMMom(){
			return KMmom.at(kmid);	
		}
		TVector3 GetKPMom(){
			return KPmom.at(kpid);	
		}
};

/*
	 class VertexLL:public Vertex{
	 private:
	 vector<Recon>HCand;
	 vector<Recon>HCorCand;
	 vector<Recon>Recons1;
	 vector<Recon>Recons2;
	 public:
	 VertexLL(Recon p){
	 Recons.push_back(p);
	 }
	 virtual int NTrack(){
	 return Recons.size()+Tracks.size();
	 }
	 bool AddTrack(Track p);
	 void SearchXiCombination();
	 Recon GetH(){
	 double comp = 9999;
	 Recon val ;
	 int num = 0;
	 for(auto Hc : HCand){ if( abs(mH-Hc.Mass())<comp) {comp=abs(mH-Hc.Mass());val=Hc;}}
	 return val;
	 }
	 Recon GetHCor(){
	 double comp = 9999;
	 Recon val ;
	 int num = 0;
	 for(auto Hc : HCorCand){ if( abs(mH-Hc.Mass())<comp) {comp=abs(mH-Hc.Mass());val=Hc;}}
	 return val;
	 }
	 };


*/
class Pi0Recon: public Recon{
	public:
		Pi0Recon(){};
		Pi0Recon(XiStarRecon XiStar,Recon Xi);
};

#endif
