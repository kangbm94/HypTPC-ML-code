
#ifndef HelixFitter_h
#define HelixFitter_h
class HelixFitter{
	protected:
		vector<TVector3> m_pos;//In Target coordinate	
		vector<TVector3> m_res;//In Target coordinate	
		double m_par[8];
		double par_fit[5];
		bool crossing = false;
		bool positive = false;
		bool negative = false;
		std::string s_hel="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
		TF1 fhel;
		TMinuit *minimizer;
	public:
		HelixFitter(vector<TVector3> posarr,double* param);
//		HelixFitter(){}
		~HelixFitter(){delete minimizer;}

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

		void HelixChi2(int& npar, double* gin, double &f, double *par,int iflag);
		void DoHelixFit();
		void GetFitPar(double* par){
			for(int i=0;i<5;++i){
				par[i] = par_fit[i];
			}
		}
};

HelixFitter::HelixFitter(vector<TVector3> posarr,double* param){
	fhel = TF1("fhel",s_hel.c_str(),-3.5,7.);
	int n = posarr.size();
	for(int i = 0;i<5;++i){
		m_par[i]=param[i];
	}
	for(int ip=0;ip<n;++ip){
		auto pos = posarr.at(ip);
		m_pos.push_back(pos);
		m_res.push_back(TVector3(1,1,1));
		double tmpt = atan2(pos.y(),pos.x());
		if(tmpt > acos(-1) * 3/4) positive = true;
		if(tmpt < acos(-1) *-3/4) negative = true;
	}
	if(positive and negative) crossing = true;
	if(crossing) fhel.SetRange( 0,2*acos(-1));
	else fhel.SetRange(-acos(-1),acos(-1));
	vector<TVector3>YTheta;
	for(int ip=0;ip<n;++ip){
		auto pos = posarr.at(ip);
		double z = pos.z();
		double tmp_t = atan2(pos.y(),pos.x());
		if(crossing and tmp_t < 0 ) tmp_t+=2*acos(-1);
		double t = tmp_t;
		TVector3 YT(t,z,0);
		YTheta.push_back(YT);
	}	
	double parl[2];
//	LinearFit(YTheta,parl);
	m_par[2]=parl[0];
	m_par[5]=parl[1];
}




static void 
HelixFitter::HelixChi2(int& npar, double* gin, double &f, double *par,int iflag){
//	int nd = m_dum.size();
	int n = m_pos.size();

	vector<double>mint;
	double chisqr=0.0;
	int dof = 0;
	double f_par[8];
	for(int i=0;i<5;++i){
		f_par[i]= m_par[i];
	}
	for(int ip=0;ip<n;++ip){
		auto pos = m_pos.at(ip);
		auto res = m_res.at(ip);
		f_par[5] =pos.x(); 
		f_par[6] =pos.y(); 
		f_par[7] =pos.z(); 
		fhel.SetParameters(f_par);
		double min_t;
		min_t = fhel.GetMinimumX();
		double xCal = f_par[0]+cos(min_t)*f_par[3];
		double yCal = f_par[1]+sin(min_t)*f_par[3];
		double zCal = f_par[2]+f_par[4]*(min_t)*f_par[3];
		TVector3 posCal(xCal,yCal,zCal);
		auto d = pos-posCal;
		auto pul = TVector3(d.x()/res.x(),d.y()/res.y(),d.z()/res.z());
		chisqr += pul.Mag();
		dof+=2;
	}
	f = chisqr /(double)(dof-5);

	}


void HelixFitter::DoHelixFit(){
	double par[8];
	double  FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
	double  LowLimit[5] = { -10000., -5000., -7000., 0., -10. };
	double  UpLimit[5] = { 10000., 5000., 7000., 10000., 10. };
	minimizer = new TMinuit(5);
	minimizer->SetPrintLevel(-1);
	minimizer->SetFCN(HelixChi2);

	for(int i=0;i<5;++i){
		par[i] = m_par[i];
	}
	int ierflg = 0;
	double arglist[10];
	double err[5] = {-999., -999., -999., -999., -999.};
	arglist[0] = 5.89;
	minimizer->mnexcm("SET ERR", arglist,1,ierflg); //Num of parameter
	arglist[0] = 1;
	minimizer->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
	TString name[5] = {"cx", "cy", "z0", "r", "dz"};
	for(int i=0; i<5; i++){
		minimizer->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
	}
	minimizer->Command("SET STRategy 0");
	arglist[0] = 1000;
	arglist[1] = 0.1;
	arglist[0] = arglist[0]*5*5*5;
	arglist[1] = arglist[1]/(10.*10.*10.);
	int Err;
	double bnd1, bnd2;
	int itry=0;
	minimizer->mnexcm("MIGRAD", arglist, 2, ierflg);
	minimizer->mnimpr();
	for(int i=0; i<5; i++){
		minimizer->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
	}
	for(int i=0;i<5;++i){
		par_fit[i] = par[i];
	}

}

#endif
