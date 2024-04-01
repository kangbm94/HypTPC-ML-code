#include "KinFit.cc"
#include "../include/XiFitter.hh"
#ifndef XiFitter_cc
#define XiFitter_cc
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void XiFitter::Initialize(){
	Clear();
	nMeas = 9;nUnkn = 3; nConst = 5; 
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double p_P = TV_P.Mag();	
	double th_P = TV_P.Theta();
	double ph_P = TV_P.Phi();
	
	mPi1 = Pi1.Mag();
	TVector3 TV_Pi1 = Pi1.Vect(); 
	double p_Pi1 = TV_Pi1.Mag();	
	double th_Pi1 = TV_Pi1.Theta();
	double ph_Pi1 = TV_Pi1.Phi();

	TVector3 TV_Pi2 = Pi2.Vect(); 
	double p_Pi2 = TV_Pi2.Mag();	
	double th_Pi2 = TV_Pi2.Theta();
	double ph_Pi2 = TV_Pi2.Phi();

	double meas[9];
	double unkn[3];
	
	

	double temp[] = {p_P,th_P,ph_P,p_Pi1,th_Pi1,ph_Pi1,p_Pi2,th_Pi2,ph_Pi2};
	for(int i=0;i<nMeas;++i)meas[i]=temp[i];
	double temp2[] = {p_Xi,th_Xi,ph_Xi};
	for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
	vector<double>Pull;
	Pull.resize(nMeas);
	vector<double>UPull;
	UPull.resize(nUnkn);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
	Chi2s.push_back(-1);
	MassDiffs.push_back(1e9);
}
void XiFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double p_Xi,th_Xi,ph_Xi,p_P,th_P,ph_P,p_Pi1,th_Pi1,ph_Pi1,p_Pi2,th_Pi2,ph_Pi2;
	p_Xi= Unkn(0,0); 
	th_Xi= Unkn(1,0); 
	ph_Xi= Unkn(2,0); 
	p_P= Meas(0,0); 
	th_P= Meas(1,0); 
	ph_P= Meas(2,0);
	p_Pi1= Meas(3,0); 
	th_Pi1= Meas(4,0); 
	ph_Pi1= Meas(5,0); 
	p_Pi2= Meas(6,0); 
	th_Pi2= Meas(7,0); 
	ph_Pi2= Meas(8,0); 


	TVector3 TV_P,TV_Pi1,TV_Pi2,TV_Xi;
	TV_P.SetMagThetaPhi(p_P,th_P,ph_P);
	TV_Pi1.SetMagThetaPhi(p_Pi1,th_Pi1,ph_Pi1);
	TV_Pi2.SetMagThetaPhi(p_Pi2,th_Pi2,ph_Pi2);
	TV_Xi.SetMagThetaPhi(p_Xi,th_Xi,ph_Xi);




	auto TV_Ld = TV_P+TV_Pi1;
	double p_Ld = TV_Ld.Mag();	
	double th_Ld = TV_Ld.Theta();
	double ph_Ld = TV_Ld.Phi();
	
	double p_Xi = TV_Xi.Mag();	
	double th_Xi = TV_Xi.Theta();
	double ph_Xi = TV_Xi.Phi();


	double E_P = hypot(mP,p_P);
	double E_Pi1 = hypot(mPi,p_Pi1);
	double E_Pi2 = hypot(mPi,p_Pi2);
	double E_Ld = hypot(mLd,p_Ld);
	double E_Xi = hypot(mXi,p_Xi);

	double px_Xi = 	p_Xi*sin(th_Xi)*cos(ph_Xi) 
	double px_P =	p_P*sin(th_P)*cos(ph_P) 
	double px_Pi1 =	p_Pi1*sin(th_Pi1)*cos(ph_Pi1)
	double px_Pi2 =	p_Pi2*sin(th_Pi2)*cos(ph_Pi2);//Constraint on x momentum
	double py_Xi =	-p_Xi*sin(th_Xi)*sin(ph_Xi) 
	double py_P = p_P*sin(th_P)*sin(ph_P) 
	double py_Pi1 =	p_Pi1*sin(th_Pi1)*sin(ph_Pi1)
	double py_Pi2 =	p_Pi2*sin(th_Pi2)*sin(ph_Pi2);//Constraint on y momentum
	double pz_Xi =	-p_Xi*cos(th_Xi) 
	double pz_P =	p_P*cos(th_P)
	double pz_Pi1 =	p_Pi1*cos(th_Pi1)
	double pz_Pi2 = p_Pi2*cos(th_Pi2);//Constraint on z momentum 

	double f1 = -px_Xi + px_P+px_Pi1+px+Pi2; 
	double f2 = -py_Xi + py_P+py_Pi1+py+Pi2; 
	double f3 = -pz_Xi + pz_P+pz_Pi1+pz+Pi2; 
	double f4	= -E_Ld + E_P + E_Pi1;
	double f5	= -E_Xi + E_Ld + E_Pi2;

	double df1du1 = -sin(th_Xi)*cos(ph_Xi);
	double df1du2 = -p_Xi*cos(th_Xi)*cos(ph_Xi);
	double df1du3 = p_Xi*sin(th_Xi)*sin(ph_Xi);
	
	double df1dm1 = sin(th_P)*cos(ph_P);
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	
	double df1dm4 = sin(th_Pi1)*cos(ph_Pi1);
	double df1dm5 = p_Pi1*cos(th_Pi1)*cos(ph_Pi1);
	double df1dm6 = -p_Pi1*sin(th_Pi1)*sin(ph_Pi1);
	
	double df1dm7 = sin(th_Pi2)*cos(ph_Pi2);
	double df1dm8 = p_Pi2*cos(th_Pi2)*cos(ph_Pi2);
	double df1dm9 = -p_Pi2*sin(th_Pi2)*sin(ph_Pi2);


	double df2du1 = -sin(th_Xi)*sin(ph_Xi);
	double df2du2 = -p_Xi*cos(th_Xi)*sin(ph_Xi);
	double df2du2 = -p_Xi*sin(th_Xi)*cos(ph_Xi);
	
	double df2dm1 = sin(th_P)*sin(ph_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	
	double df2dm4 = sin(th_Pi1)*sin(ph_Pi1);
	double df2dm5 = p_Pi1*cos(th_Pi1)*sin(ph_Pi1);
	double df2dm6 = p_Pi1*sin(th_Pi1)*cos(ph_Pi1);
	
	double df2dm7 = sin(th_Pi2i1)*sin(ph_Pi2i1);
	double df2dm8 = p_Pi2*cos(th_Pi2)*sin(ph_Pi2);
	double df2dm9 = p_Pi2*sin(th_Pi2)*cos(ph_Pi2);
	

	double df3du1 = -cos(th_Xi);
	double df3du2 = p_Xi;
	double df3du3 = 0;

	double df3dm1 = cos(th_P);
	double df3dm2 = -p_P;
	double df3dm3 = 0;
	
	double df3dm4 = cos(th_Pi1);
	double df3dm5 = -p_Pi1;
	double df3dm6 = 0;
	
	double df3dm7 = cos(th_Pi2);
	double df3dm8 = -p_Pi2;
	double df3dm9 = 0;

	double PPi1 = TV_P*TV_Pi1;
	double dp_Lddp_P = (p_P + PPi1/p_P) / p_Ld;
	double dp_Lddth_P =(
		(px_P+px_Pi1)*cos(th_P)*cos(ph_P)
		+(py_P+py_Pi1)*cos(th_P)*sin(ph_P)
		+(pz_P+pz_Pi1)*(-sin(th_P))
			)  /p_Ld;// dp_Ld / dth_P =( px_Ld*dpx_Ld/dth_P + ...)  /p_Ld ...; dpx_Ld/dth_P = dpx_P/dth_P ...

	double dp_Lddp_Pi1 = (p_Pi1 + PPi1/p_Pi1) / p_Ld;





	double df4du1 = 0; 
	double df4du2 = 0; 
	double df4du3 = 0; 

	double df4dm1 = -p_Ld/E_Ld *  


	double fs[]={f1,f2,f3,f4,f5};
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
		double tempu[] = {
			df1du1,
			df2du1,
			df3du1,
			df4du1
		};	
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[]={
			df1du1du1
		};
		double temp2[]={
			df2du1du1
		};
		double temp3[]={
			df3du1du1
		};
		double temp4[]={
			df4du1du1
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);
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
		double tempu[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[] = {
			df1du1du1,df1du1du2,df1du1du3,
			df1du2du1,df1du2du2,df1du2du3,
			df1du3du1,df1du3du2,df1du3du3
		};
		double temp2[] = {
			df2du1du1,df2du1du2,df2du1du3,
			df2du2du1,df2du2du2,df2du2du3,
			df2du3du1,df2du3du2,df2du3du3
		};
		double temp3[] = {
			df3du1du1,df3du1du2,df3du1du3,
			df3du2du1,df3du2du2,df3du2du3,
			df3du3du1,df3du3du2,df3du3du3
		};
		double temp4[] = {
			df4du1du1,df4du1du2,df4du1du3,
			df4du2du1,df4du2du2,df4du2du3,
			df4du3du1,df4du3du2,df4du3du3
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);


	}

	



	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void XiFitter::SampleStepPoint(int steps){
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double p_P,th_P,ph_P,p_Pi1,th_Pi1,ph_Pi1,p_R,th_R,ph_R;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Pi1= Meas(5,0); 
		th_Pi1= Meas(6,0); 
		ph_Pi1= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Pi1= Meas(3,0); 
		th_Pi1= Meas(4,0); 
		ph_Pi1= Meas(5,0); 
	}
	double Ppx = p_P*sin(th_P)*cos(ph_P);
	double Ppy = p_P*sin(th_P)*sin(ph_P);
	double Ppz = p_P*cos(th_P);
	double Pi1px = p_Pi1*sin(th_Pi1)*cos(ph_Pi1);
	double Pi1py = p_Pi1*sin(th_Pi1)*sin(ph_Pi1);
	double Pi1pz = p_Pi1*cos(th_Pi1);
	double Rpx = p_R*sin(th_R)*cos(ph_R);
	double Rpy = p_R*sin(th_R)*sin(ph_R);
	double Rpz = p_R*cos(th_R);
	TLorentzVector PP(Ppx,Ppy,Ppz,hypot(mP,p_P));
	TLorentzVector Pi1Pi1(Pi1px,Pi1py,Pi1pz,hypot(mPi1,p_Pi1));
	TLorentzVector RR(Rpx,Rpy,Rpz,hypot(mR,p_R));
	auto V = PP + Pi1Pi1;
	double MassDiff = V.Mag()-mR;
	PCor = PP;
	Pi1Cor = Pi1Pi1;
	RCor = RR;
	MassDiffs.push_back(MassDiff);
}
XiFitter::XiFitter(TLorentzVector P_,TLorentzVector Pi1_, TLorentzVector R_){
	//Does Kinematic fitting for R -> P+Pi1 decay.
	P=P_;
	Pi1=Pi1_;
	R=R_;
	Initialize();
};
TMatrixD
XiFitter::JacobianSphToCart(double p, double th, double ph){
	// x = p sin(th) cos(ph)
	// y = p sin(th) sin(ph)
	// z = p cos(th)
	//V_c = J^T V J|->
	//		dxdp, dxdth,dxdph
	//J	=	dydp, dydth,dydph
	//		dzdp, dzdth,dzdph

	double dxdp = sin(th)*cos(ph);
	double dydp = sin(th)*sin(ph);
	double dzdp = cos(th);

	double dxdth = p*cos(th)*cos(ph);
	double dydth = p*cos(th)*sin(ph);
	double dzdth = -p*sin(th);

	double dxdph = -p*sin(th)*sin(ph);
	double dydph = p*sin(th)*cos(ph);
	double dzdph = 0;
	double mat[9] = 
	{ dxdp, dxdth, dxdph,
		dydp, dydth, dydph,
		dzdp, dzdth, dzdph
	};
	/*
	double mat[9] = 
	{ dxdp, dydp, dzdp,
		dxdth, dydth, dzdth,
		dxdph, dydph, dzdph
	};
	*/
	return TMatrixD(3,3,mat);

}
void
XiFitter::CalcVariance(int istep){
	auto Meas = Measurements.at(istep); 
	auto Unkn = Unknowns.at(istep);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Pi1,th_Pi1,ph_Pi1;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Pi1= Meas(5,0); 
		th_Pi1= Meas(6,0); 
		ph_Pi1= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Pi1= Meas(3,0); 
		th_Pi1= Meas(4,0); 
		ph_Pi1= Meas(5,0); 
	}
	TMatrixD Jsc_P = JacobianSphToCart(p_P,th_P,ph_P);
	TMatrixD Jsc_Pi1 = JacobianSphToCart(p_Pi1,th_Pi1,ph_Pi1);
	
	double El_Jsc_PPi1[6*6]= {0};
	for(int ic =0;ic<3;++ic){
	for(int ir =0;ir<3;++ir){
		int col_P = ic, row_P = ir;
		int col_Pi1 = ic+3, row_Pi1 = ir+3;
		El_Jsc_PPi1[row_P+6*col_P] = Jsc_P(ic,ir);
		El_Jsc_PPi1[row_Pi1+6*col_Pi1] = Jsc_Pi1(ic,ir);
	}
	}
	TMatrixD Jsc_PPi1 = TMatrixD(6,6,El_Jsc_PPi1);
	/*
	Jsc_P.Print();
	Jsc_Pi1.Print();
	Jsc_PPi1.Print();
	cin.ignore();
	*/
	TMatrixD Jsc_PPi1_T = TransposeMatrix(Jsc_PPi1);
	TMatrixD VMat = Variancies.at(istep);
	TMatrixD dV = dVMats.at(istep);
	TMatrixD VMat_C = Jsc_PPi1_T*(VMat-dV)*Jsc_PPi1;
//	TMatrixD VMat_C = Jsc_PPi1_T*(VMat)*Jsc_PPi1;
	double El_contract[18]={
		1,0,0,1,0,0,
		0,1,0,0,1,0,
		0,0,1,0,0,1
	};
	TMatrixD ContT(3,6,El_contract);
	TMatrixD Cont = TransposeMatrix(ContT);
	TMatrixD UVMat_C = ContT*VMat_C*Cont;
	TMatrixD Jcs_R = JacobianSphToCart(p_R,th_R,ph_R);
	Jcs_R.Invert();
	auto Jcs_RT = TransposeMatrix(Jcs_R);
	TMatrixD UVMat = Jcs_RT*UVMat_C*Jcs_R;
//	VarianciesU.push_back(UVMat);
}















void
XiFitter::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
XiFitter::ToDecayPlane(){
	auto Zaxis =(P + Pi1).Vect();
	auto vP = P.Vect();
	auto vPi1 = Pi1.Vect();
	auto Yaxis = vP.Cross(vPi1);
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
