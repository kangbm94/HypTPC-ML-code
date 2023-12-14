#include "KinFit2.cc"
#include "../include/FourVectorFitter.hh"
#ifndef LVFit_cc
#define LVFit_cc
void FourVectorFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	MeasDir = status;
	Clear();
	if(MeasDir){//Production and Decay vertex of R is known, hence the direction is known. However, the momentum magnitude is not directly measured. I used the sum of P and Q as an initial value.
		auto dir = Vert1 - Vert2;	
		auto th = dir.Theta();
		auto ph = dir.Phi();
		R.SetTheta(th);
		R.SetPhi(ph);
	}
	Initialize();
}
void FourVectorFitter::Initialize(){
	Clear();
	if(MeasDir){
		nMeas = 8;nUnkn = 1; nConst = 4; 
	}
	else{
		nMeas = 6;nUnkn = 3; nConst = 4; //For 1-C fit, when Vertex information is useless.
	}
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double p_P = TV_P.Mag();	
	double th_P = TV_P.Theta();
	double ph_P = TV_P.Phi();
	
	mQ = Q.Mag();
	TVector3 TV_Q = Q.Vect(); 
	double p_Q = TV_Q.Mag();	
	double th_Q = TV_Q.Theta();
	double ph_Q = TV_Q.Phi();

	TVector3 TV_R = R.Vect(); 
	double p_R = TV_R.Mag();	
	double th_R = TV_R.Theta();
	double ph_R = TV_R.Phi();

	double meas[8];
	double unkn[3];
	if(MeasDir){
		double temp[] = {th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_R,0,0};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	else{
		double temp[] = {p_P,th_P,ph_P,p_Q,th_Q,ph_Q,0,0};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_R,th_R,ph_R};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
	vector<double>Pull;
	Pull.resize(nMeas);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	Chi2s.push_back(-1);
	MassDiffs.push_back(1e9);
}
void FourVectorFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Q= Meas(5,0); 
		th_Q= Meas(6,0); 
		ph_Q= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 
	}
	
	
	double f1 = 
		-p_R*sin(th_R)*cos(ph_R) 
		+p_P*sin(th_P)*cos(ph_P) 
		+p_Q*sin(th_Q)*cos(ph_Q) ;//Constraint on x momentum
	double f2 = 
		-p_R*sin(th_R)*sin(ph_R) 
		+p_P*sin(th_P)*sin(ph_P) 
		+p_Q*sin(th_Q)*sin(ph_Q) ;//Constraint on y momentum
	double f3 =  
		-p_R*cos(th_R) 
		+p_P*cos(th_P)
		+p_Q*cos(th_Q);//Constraint on z momentum 
	double f4	=
		- sqrt(p_R*p_R+mR*mR)
		+ sqrt(p_P*p_P+mP*mP)
		+ sqrt(p_Q*p_Q+mQ*mQ);//Constraint on Energy

	double df1du1 = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
	double df1du2 = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1, in case of 3-C fit.However, I didnt want to change the token name... Mathematically it should be df1 / d (Th_R)
	double df1du3 = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
	double df1dm1 = sin(th_P)*cos(ph_P);//...
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df1dm4 = sin(th_Q)*cos(ph_Q);
	double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
	double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	
	double df2du1 = -sin(th_R)*sin(ph_R);
	double df2du2 = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
	double df2du3 = -p_R*sin(th_R)*cos(ph_R);
	double df2dm1 = sin(th_P)*sin(ph_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	double df2dm4 = sin(th_Q)*sin(ph_Q);
	double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
	double df2dm6 = p_Q*sin(th_Q)*cos(ph_Q);

	double df3du1 = -cos(th_R);
	double df3du2 = p_R*sin(th_R);
	double df3du3 = 0;
	double df3dm1 = cos(th_P);
	double df3dm2 = -p_Q*sin(th_P);
	double df3dm3 = 0;
	double df3dm4 = cos(th_Q);
	double df3dm5 = -p_Q*sin(th_Q);
	double df3dm6 = 0;

	double df4du1 = -p_R/sqrt(p_R*p_R+mR*mR);
	double df4du2 = 0;
	double df4du3 = 0;
	double df4dm1 = p_P/sqrt(p_P*p_P+mP*mP);
	double df4dm2 = 0;
	double df4dm3 = 0;
	double df4dm4 = p_Q/sqrt(p_Q*p_Q+mQ*mQ);
	double df4dm5 = 0;
	double df4dm6 = 0;

	double fs[]={f1,f2,f3,f4};
	double dfdms[200] ;
	double dfdus[200] ;
	if(MeasDir){
		double temp[] = {
			df1du2,df1du3,df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,
			df2du2,df2du3,df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,
			df3du2,df3du3,df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,
			df4du2,df4du3,df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double temp2[] = {
			df1du1,
			df2du1,
			df3du1,
			df4du1
		};	
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=temp2[i];
		};
	}
	else{
		double temp[] = {
			df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,
			df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,
			df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,
			df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double temp2[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=temp2[i];
		};
	}
	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void FourVectorFitter::SampleStepPoint(int steps){
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_R,th_R,ph_R;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Q= Meas(5,0); 
		th_Q= Meas(6,0); 
		ph_Q= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 
	}
	double Ppx = p_P*sin(th_P)*cos(ph_P);
	double Ppy = p_P*sin(th_P)*sin(ph_P);
	double Ppz = p_P*cos(th_P);
	double Qpx = p_Q*sin(th_Q)*cos(ph_Q);
	double Qpy = p_Q*sin(th_Q)*sin(ph_Q);
	double Qpz = p_Q*cos(th_Q);
	double Rpx = p_R*sin(th_R)*cos(ph_R);
	double Rpy = p_R*sin(th_R)*sin(ph_R);
	double Rpz = p_R*cos(th_R);
	TLorentzVector PP(Ppx,Ppy,Ppz,hypot(mP,p_P));
	TLorentzVector QQ(Qpx,Qpy,Qpz,hypot(mQ,p_Q));
	TLorentzVector RR(Rpx,Rpy,Rpz,hypot(mR,p_R));
	auto V = PP + QQ;
	double MassDiff = V.Mag()-mR;
	PCor = PP;
	QCor = QQ;
	RCor = RR;
	MassDiffs.push_back(MassDiff);
}
FourVectorFitter::FourVectorFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	//Does Kinematic fitting for R -> P+Q decay.
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
void
FourVectorFitter::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
FourVectorFitter::ToDecayPlane(){
	auto Zaxis =(P + Q).Vect();
	auto vP = P.Vect();
	auto vQ = Q.Vect();
	auto Yaxis = vP.Cross(vQ);
//	double YNorm = 1./(Yaxis.Mag());
//	Yaxis = YNorm * Yaxis;
	double Th_F = Zaxis.Theta();
	double Ph_F = Zaxis.Phi();
	double RotZ[9] ={
		cos(-Ph_F),	-sin(-Ph_F),	0,	
		sin(-Ph_F),	cos(-Ph_F),		0,
		0,					0,						1
	};
	double RotY[9] ={
		cos(Th_F),	0,					-sin(Th_F),
		0,					-1,					0,
		sin(Th_F),	0,					cos(Th_F)
	};
	TMatrixD RZ(3,3,RotZ);
	TMatrixD RY(3,3,RotY);
	TMatrixD R_F = RY * RZ;
	Yaxis = R_F * Yaxis;
	double Th_Y = Yaxis.Theta();
	double Ph_Y = Yaxis.Phi();
	double RotX[9] ={
		1,				0,					0,
		0,				cos(-Ph_Y),	-sin(-Ph_Y),
		0,				sin(-Ph_Y),	cos(-Ph_Y)
	};
	TMatrixD RX(3,3,RotX);
	Yaxis = RX * Yaxis;

}
#endif
