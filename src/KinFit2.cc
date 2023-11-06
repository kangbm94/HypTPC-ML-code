#include "../include/KinFit2.hh"
#ifndef KinFit_cc
#define KinFit_cc
#define Debug 0
#define ShowChi2 0
void KinematicFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	MeasDir = status;
	Measurements.clear();
	Unknowns.clear();
	if(MeasDir){
		auto dir = Vert1 - Vert2;	
		auto th = dir.Theta();
		auto ph = dir.Phi();
		R.SetTheta(th);
		R.SetPhi(ph);
	}
	Initialize();
}
void KinematicFitter::Initialize(){
	if(MeasDir){
		nMeas = 8;nUnkn = 1; nConst = 4; 
	}
	else{
		nMeas = 6;nUnkn = 3; nConst = 4; 
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
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
}

void KinematicFitter::Finalize(){
	for(int is=0;is<step+1;++is){
		double chi2 = Chi2s.at(is);
#if ShowChi2
		cout<<Form("Step %d, chisqr %g",is,Chi2s.at(is))<<endl;
#endif
		if(chi2<Best_Chi2){
			Best_Chi2 = chi2;
			best_step = is;
		}
	}
	auto Meas = Measurements.at(best_step); 
	auto Unkn = Unknowns.at(best_step); 
	double Chi2 = Chi2s.at(best_step);
		
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
	TVector3 TV_P = TVector3(p_P*sin(th_P)*cos(ph_P),
			p_P*sin(th_P)*sin(ph_P),
			p_P*cos(th_P));
	PCor = TLorentzVector(TV_P.x(),TV_P.y(),TV_P.z(),hypot(mP,p_P));

	TVector3 TV_Q = TVector3(p_Q*sin(th_Q)*cos(ph_Q),
			p_Q*sin(th_Q)*sin(ph_Q),
			p_Q*cos(th_Q));
	QCor = TLorentzVector(TV_Q.x(),TV_Q.y(),TV_Q.z(),hypot(mQ,p_Q));
	
	TVector3 TV_R = TVector3(p_R*sin(th_R)*cos(ph_R),
			p_R*sin(th_R)*sin(ph_R),
			p_R*cos(th_R));
	RCor = TLorentzVector(TV_R.x(),TV_R.y(),TV_R.z(),hypot(mR,p_R));
}
void KinematicFitter::SetConstraints(){
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
		+p_Q*sin(th_Q)*cos(ph_Q) ;//px
	double f2 = 
		-p_R*sin(th_R)*sin(ph_R) 
		+p_P*sin(th_P)*sin(ph_P) 
		+p_Q*sin(th_Q)*sin(ph_Q) ;// py
	double f3 =  
		-p_R*cos(th_R) 
		+p_P*cos(th_P)
		+p_Q*cos(th_Q);// pz
	double f4	=
		- sqrt(p_R*p_R+mR*mR)
		+ sqrt(p_P*p_P+mP*mP)
		+ sqrt(p_Q*p_Q+mQ*mQ);//E0
	
	double df1du1 = -sin(th_R)*cos(th_R);
	double df1du2 = -p_R*cos(th_R)*cos(ph_R);// d/ dth_R
	double df1du3 = p_R*sin(th_R)*sin(ph_R);
	double df1dm1 = sin(th_P)*cos(th_P);
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df1dm4 = sin(th_Q)*cos(th_Q);
	double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
	double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	
	double df2du1 = -sin(th_R)*sin(th_R);
	double df2du2 = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
	double df2du3 = -p_R*sin(th_R)*cos(ph_R);
	double df2dm1 = sin(th_P)*sin(th_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	double df2dm4 = sin(th_Q)*sin(th_Q);
	double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
	double df2dm6 = -p_Q*sin(th_Q)*cos(ph_Q);

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
	FMats.push_back(FMat);
	dFdMs.push_back(dFdM);
	dFdUs.push_back(dFdU);
}



KinematicFitter::KinematicFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
void
KinematicFitter::SetVariance(double* var){
	double variance[200];
	for(int i = 0;i<nMeas*nMeas;++i){
		variance[i]=0;
	}
	for(int i = 0;i<nMeas;++i){
		variance[i+i*nMeas] = var[i];
	}
	TMatrixD Variance(nMeas,nMeas,variance);//Actually, the inverse of the variance
	Variancies.push_back(Variance);
};
void KinematicFitter::ProcessStep(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step); 
	auto Meas0 = Measurements.at(0);
	auto Unkn0 = Unknowns.at(0);
	auto VMat = Variancies.at(0);
#if Debug
	cout<<"Step: "<<step<<endl;
	cout<<"Meas Mat"<<endl;
	Meas.Print();
	cout<<"Unkn Mat"<<endl;
	Unkn.Print();
	cout<<"V Mat"<<endl;
	VMat.Print();
#endif
	SetConstraints();

	TMatrixD FMat = FMats.at(step);
	TMatrixD dFdM = dFdMs.at(step);
	TMatrixD dFdU = dFdUs.at(step);
	double det_dFdM=0,det_dFdU=0;
	auto dFdMT = TransposeMatrix(dFdM);
	auto dFdUT = TransposeMatrix(dFdU);
	auto rMat = FMat + dFdM*(Meas0-Meas);
	auto sMat = dFdM*VMat*dFdMT;
#if Debug
	cout<<"F Mat"<<endl;
	FMat.Print();
	cout<<"dFdM Mat"<<endl;
	dFdM.Print();
	cout<<"dFdM Transpose"<<endl;
	dFdMT.Print();
	cout<<"dFdU Mat"<<endl;
	dFdU.Print();
	cout<<"dFdU Transpose"<<endl;
	dFdUT.Print();
	cout<<"r Mat"<<endl;
	rMat.Print();
	cout<<"s Mat"<<endl;
	sMat.Print();
#endif
	
	auto Unkn_next = Unkn - (dFdUT*(sMat.Invert())*dFdU).Invert() *dFdUT* sMat.Invert() * rMat;

#if Debug
	cout<<"Unkn Mat"<<endl;
	Unkn.Print();
	cout<<"UnknNext Mat"<<endl;
	Unkn_next.Print();
#endif
	auto Lambda = (sMat.Invert()* (rMat+ dFdU*(Unkn_next - Unkn)));
	auto LambdaT = TransposeMatrix(Lambda);
	
	auto VFL =  VMat*dFdMT*Lambda;
	auto Meas_next = Meas0 - VFL; 
	
#if Debug
	cout<<"Meas0 Mat"<<endl;
	Meas0.Print();
	cout<<"Lambda Mat"<<endl;
	Lambda.Print();
	cout<<"VFL Mat"<<endl;
	VFL.Print();
	cout<<"MeasNext Mat"<<endl;
	Meas_next.Print();
#endif
	TMatrixD dM = Meas0-Meas;
	auto dMT = TransposeMatrix(dM);
	double Chi2 = (dMT* (VMat.Invert())*dM)(0,0) + 2 * (LambdaT * FMat )(0,0);
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(Meas_next);	
#if Debug
	cout<<"dM "<<endl;
	dM.Print();
	cout<<"dM Transpose "<<endl;
	dMT.Print();
#endif
	rMats.push_back(rMat);
	sMats.push_back(sMat);
//	Lambdas.push_back(Lambda);
	Chi2s.push_back(Chi2);
}


double KinematicFitter::DoKinematicFit(){
	int Cnt = 0;
	while(1){
		ProcessStep();
		double Chi2 = Chi2s.at(step);
		double Chi2_prev = -1;
		if(step > 0 )Chi2_prev = Chi2s.at(step-1);
		if(Cnt > 3) break;
		if(abs(Chi2_prev - Chi2) < 0.1){
			Cnt++;
		}
		else{
			Cnt= 0;
		}
		if(step > MaxStep-1){
			break;
		}
		step++;
	}
	Finalize();
//	return Best_Chi2;
	return Chi2s.at(0);
}
void
KinematicFitter::Clear(){
}

TMatrixD
KinematicFitter::TransposeMatrix(TMatrixD M){
	int row = M.GetNrows();
	int col = M.GetNcols();
	double elem[500];
	for(int r=0;r<row;++r){
	for(int c=0;c<col;++c){
		elem[c*row + r]= M(r,c); 
	}	
	}
	return TMatrixD(col,row,elem);
}
#endif
