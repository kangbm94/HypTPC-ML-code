#ifndef HelixTrack_h
#define HelixTrack_h
#include "Math.hh"
#include "Kinematics.hh"
namespace{
	int nitr = 0;
	vector<TVector3> m_pos;//In Target coordinate	
	vector<TVector3> m_res;//In Target coordinate	
	std::string s_hel="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
	TF1 fhel  = TF1("fhel",s_hel.c_str(),-3.5,7.);
	double m_par[8];
		bool crossing = false;
		bool positive = false;
		bool negative = false;
	void HelixChi2(int& npar, double* gin, double &f, double *par,int iflag){
		int n = m_pos.size();

		vector<double>mint;
		double chisqr=0.0;
		int dof = 0;
		double f_par[8];
		for(int i=0;i<5;++i){
			f_par[i]= par[i];
		}
		for(int ip=0;ip<n;++ip){
			auto pos = m_pos.at(ip);
			auto res = m_res.at(ip);
			f_par[5] =pos.x(); 
			f_par[6] =pos.y(); 
			f_par[7] =pos.z(); 
			fhel.SetParameters(f_par);
//			double min_t = atan2(pos.y()-f_par[1],pos.x()-f_par[0]);
//			if(crossing and min_t<0) min_t+=2*acos(-1); 	
			double min_t;
			if(crossing)min_t= fhel.GetMinimumX(0,2*acos(-1));
			else min_t = fhel.GetMinimumX(-acos(-1),acos(-1));
			double xCal = f_par[0]+cos(min_t)*f_par[3];
			double yCal = f_par[1]+sin(min_t)*f_par[3];
			double zCal = f_par[2]+f_par[4]*(min_t)*f_par[3];
			TVector3 posCal(xCal,yCal,zCal);
			auto d = pos-posCal;
			auto pul = d;
			chisqr += pul.Mag();
			dof+=2;
		}
		f = chisqr /(double)(dof-5);
//		TString parl = Form(", pars: cx,cy,z0,r,dz = (%f,%f,%f,%f,%f)",par[0],par[1],par[2],par[3],par[4]);
		nitr++;
	};
}
class HelixTrack{
	protected:
		double par_fit[5];
		TMinuit *minimizer;
		double mom0=0;
		double mom0Fit=0;
		vector<int> m_l;
		int charge;
	public:
		HelixTrack(vector<TVector3> posarr,double* param);
		HelixTrack(double* param,int q=0);
		HelixTrack(){}
		~HelixTrack(){}

		void SetResolution(vector<TVector3> resarr){
			int n = resarr.size();
			if(n != m_res.size()){
				cout<<"Warning! Nres != Npos"<<endl;
				return;
			}
			for(int i = 0;i<n;++i){
				m_res.at(i) = resarr.at(i);
			}
		}
		bool DoHelixFit();
		void GetFitPar(double* par){
			for(int i=0;i<5;++i){
				par[i] = par_fit[i];
			}
		}
		void Initialize(){
			m_pos.clear();
			m_res.clear();
			crossing = false;
			negative = false;
			positive = false;
			for(int i=0;i<8;++i){
				m_par[i]=0;
			}
		}
		double CircDist(TVector3 point){
			auto pt = GlobalToTarget(point);
			double x = pt.x(),y = pt.y();
			double dist = hypot(y - m_par[1],x-m_par[0]) - m_par[3];
			return dist;	
		}
		double HelixDist(TVector3 pos){
			auto t = GetTcal(m_par,pos);
			auto v = HelixPos(m_par,t);
			double dist =(pos-v).Mag();
//			cout<<dist<<endl;
			return dist;
		}
		void AddHit(TVector3 point,int l=0){
			m_pos.push_back(GlobalToTarget(point));
			m_l.push_back(l);
			m_res.push_back(TVector3(0.3,0.3,0.3));
		}
		double GetMom0(){
			return mom0;
		}
		double GetMom0Fit(){
			return mom0Fit;
		}
		int GetNHits(){
			return m_pos.size();
		}
		TVector3 GetPosition(int i){
			return m_pos.at(i);
		}
		int GetLayer(int i){
			return m_l.at(i);
		}
};

HelixTrack::HelixTrack(double* param,int q){
	Initialize();
	charge = q;
	for(int i = 0;i<5;++i){
		m_par[i]=param[i];
	}
	mom0 = RadToMom(m_par[3],m_par[4]); 
}
HelixTrack::HelixTrack(vector<TVector3> posarr,double* param){
	Initialize();
	for(int i = 0;i<5;++i){
		m_par[i]=param[i];
	}
	for(auto pos:posarr){
		auto pos_tgt = GlobalToTarget(pos);
		m_pos.push_back(pos_tgt);
		m_res.push_back(TVector3(0.3,0.3,0.3));
	}
}




bool HelixTrack::DoHelixFit(){
	int n = m_pos.size();
//	cout<<"HelixFit:: npoints = "<<n<<endl;
	if(n<8) return false;
	double cir_pars[5];
	CircleFitWithRadius(m_pos,m_par[3],cir_pars,charge);
	m_par[0] = cir_pars[0];	
	m_par[1] = cir_pars[1];	
	for(auto pos:m_pos){
		double tmpt = atan2(pos.y()-m_par[1],pos.x()-m_par[0]);
		if(tmpt > acos(-1) * 3/4) positive = true;
		if(tmpt < acos(-1) *-3/4) negative = true;
	}
	if(positive and negative) crossing = true;
	if(crossing) fhel.SetRange( 0,2*acos(-1));
	else fhel.SetRange(-acos(-1),acos(-1));
	vector<TVector3>YTheta;
	for(int ip=0;ip<n;++ip){
		auto pos = m_pos.at(ip);
		double z = pos.z();
		double tmp_t = atan2(pos.y()-m_par[1],pos.x()-m_par[0]);
		if(crossing and tmp_t < 0 ) tmp_t+=2*acos(-1);
		double t = tmp_t;
		TVector3 YT(t,z,0);
		YTheta.push_back(YT);
	}	
	double parl[2];

	LinearFit(YTheta,parl);
	m_par[2]=parl[0];
	m_par[4]=parl[1]/m_par[3];
	

//	double  FitStep[5] = { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-6 };
	double  FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
	double  LowLimit[5] = { -10000., -5000., -7000., 0., -10. };
	double  UpLimit[5] = { 10000., 5000., 7000., 10000., 10. };
	minimizer = new TMinuit(5);
	minimizer->SetPrintLevel(-1);
	minimizer->SetFCN(HelixChi2);
	nitr=0;
	int ierflg = 0;
	double arglist[10];
	double err[5] = {-999., -999., -999., -999., -999.};
	arglist[0] = 5.89;
	minimizer->mnexcm("SET ERR", arglist,1,ierflg); //Num of parameter
	arglist[0] = 1;
	minimizer->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
	TString name[5] = {"cx", "cy", "z0", "r", "dz"};
	for(int i=0; i<5; i++){
		minimizer->mnparm(i, name[i], m_par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
	}
	minimizer->mnparm(3,name[3],m_par[3],FitStep[3],0.7 * m_par[3],1.3*m_par[3],ierflg);
	minimizer->mnparm(4,name[4],m_par[4],FitStep[4], m_par[4]-0.1,m_par[4]+0.1,ierflg);
	minimizer->Command("SET STRategy 2");
	arglist[0] = 100;
	arglist[1] = 0.01;
//	arglist[0] = arglist[0]*5*5*5;
//	arglist[1] = arglist[1]/(10.*10.*10.);
	int Err;
	double bnd1, bnd2;
	int itry=0;
	minimizer->mnexcm("MIGRAD", arglist, 2, ierflg);
	minimizer->mnimpr();
	for(int i=0; i<5; i++){
		minimizer->mnpout(i, name[i], m_par[i], err[i], bnd1, bnd2, Err);
	}
	for(int i=0;i<5;++i){
		par_fit[i] = m_par[i];
	}
	mom0Fit = RadToMom(par_fit[3],par_fit[4]); 
	delete minimizer;
	return true;
}

#endif
