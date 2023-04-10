#include "../include/KinFit.hh"
#ifndef KinFit_cc
#define KinFit_cc
void
KinematicFitter::AssignLorentzVector(TLorentzVector P_,TLorentzVector Q_){
	P=P_;
	Q=Q_;
};

void
KinematicFitter::SetResolution(TVector3 Pres_,TVector3 Qres_){
	for(int i = 0;i<3;++i){
		Pres(i)=abs(Pres_(i));
		Qres(i)=abs(Qres_(i));
	}
};
vector<TLorentzVector>
KinematicFitter::GetFittedLV(){
	vector<TLorentzVector> Vs = {PCor,QCor};
	return Vs;
}


void KinematicFitter::DoKinematicFit(){
	int npar = 7;
	double px0 = P.X(),py0=P.Y(),pz0=P.Z();
	double qx0 = Q.X(),qy0=Q.Y(),qz0=Q.Z();
	double psx = Pres.X(),psy=Pres.Y(),psz=Pres.Z(); 
	double qsx = Qres.X(),qsy=Qres.Y(),qsz=Qres.Z(); 
	double MP=P.Mag(),MQ=Q.Mag();
	double par[15]={px0,py0,pz0,qx0,qy0,qz0,
							psx,psy,psz,qsx,qsy,qsz,
							MP,MQ,InvMass};
	for(int i=0;i<15;++i){
		gpar[i] = par[i];
	}
	minimizer->SetFCN(fcn);
	minimizer->SetPrintLevel(-1);
	double arglist[30];
	int ierflg = 0;
	arglist[0]=1;
	double step = 0.001;
//	minimizer->mnexcm("SET ERR",arglist, 1, ierflg);
	static double par_inni[7] = {0,
							px0,py0,pz0,qx0,qy0,qz0};//Initial value 
	static double par_step[7] = {1e-4,
		step*psx,step*psy,step*psz,step*qsx,step*qsy,step*qsz};
	TString parname [7] = {"lambda","px","py","pz","qx","qy","qz"};
	for(int i=0;i<npar;++i){
		minimizer->mnparm(i,parname[i],par_inni[i],par_step[i],0,0,ierflg);
	}
	arglist[0]=1000;//maxcalls
	arglist[1]=0.01;//tolerance
	minimizer->mnexcm("MIGRAD",arglist,2,ierflg);
	double px,py,pz,qx,qy,qz;
	double par_cor[7];
	double par_err[7];
	for(int i=0;i<7;++i){
		minimizer->GetParameter(i,par_cor[i],par_err[i]);
	}
	lambda = par_cor[0];
	px = par_cor[1];
	py = par_cor[2];
	pz = par_cor[3];
	qx = par_cor[4];
	qy = par_cor[5];
	qz = par_cor[6];
//	cout<<Form("Corrected: labmda = %f",lambda)<<endl;
//	cout<<Form("Corrected: P -> PCor = (%f,%f,%f)->(%f,%f,%f)",px0,py0,pz0,px,py,pz)<<endl;
//	cout<<Form("Corrected: Q -> QCor = (%f,%f,%f)->(%f,%f,%f)",qx0,qy0,qz0,qx,qy,qz)<<endl;

	PCor.SetXYZM(px,py,pz, MP);
	QCor.SetXYZM(qx,qy,qz, MQ);
}
#endif
