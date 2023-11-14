#include "../include/KinFit2.hh"
#ifndef KinFit_cc
#define KinFit_cc
#define Debug 0
#define ShowChi2 0
void
KinematicFitter::SetVariance(double* var){
	double variance[200];
	for(int i = 0;i<nMeas*nMeas;++i){
		variance[i]=0;
	}
	for(int i = 0;i<nMeas;++i){
		variance[i+i*nMeas] = var[i];
	}
	TMatrixD Variance(nMeas,nMeas,variance);
	Variancies.push_back(Variance);
};
void KinematicFitter::ProcessStep(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step); 
	auto Meas0 = Measurements.at(0);
	auto Unkn0 = Unknowns.at(0);
	auto VMat = Variancies.at(step);
	auto VInv = VMat;
	VInv.Invert();
#if Debug
	cout<<"Step: "<<step<<endl;
	cout<<"Meas Mat";
	Meas.Print();
	cout<<"Unkn Mat";
	Unkn.Print();
	cout<<"V Mat";
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
	auto sInv = sMat;
	sInv.Invert();
#if Debug
	cout<<"F Mat";
	FMat.Print();
	cout<<"dFdM Mat";
	dFdM.Print();
	cout<<"dFdU Mat";
	dFdU.Print();
	cout<<"r Mat";
	rMat.Print();
	cout<<"s Mat";
	sMat.Print();
	cout<<"s Inv";
	sInv.Print();
#endif
	auto FuSIFu = dFdUT*sInv*dFdU;
	FuSIFu.Invert();
	auto dU = (FuSIFu) * (dFdUT* (sInv) * rMat) ;
	dU = dU -dU - dU;
	auto Unkn_next = Unkn + dU;

#if Debug
#endif
	auto Lambda = (sInv* (rMat+ dFdU*dU));
	auto LambdaT = TransposeMatrix(Lambda);
	
	auto VFL =  VMat*dFdMT*Lambda;
	auto Meas_next = Meas0 - VFL; 
	
#if Debug
	cout<<"VFL Mat";
	VFL.Print();
	cout<<"dU Mat";
	dU.Print();
	cout<<"UnknNext Mat";
	Unkn_next.Print();
	cout<<"MeasNext Mat";
	Meas_next.Print();
#endif
	TMatrixD dM = Meas0-Meas;
	auto dMT = TransposeMatrix(dM);
	double Chi2 = (dMT* (VInv)*dM)(0,0) + 2 * (LambdaT * FMat )(0,0);
	
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(Meas_next);	
#if Debug
	cout<<"dM ";
	dM.Print();
	cout<<"Lambda Mat";
	Lambda.Print();
#endif
	rMats.push_back(rMat);
	sMats.push_back(sMat);
//	Lambdas.push_back(Lambda);
	Chi2s.push_back(Chi2);
	auto VMat_next = VMat;
	Variancies.push_back(VMat_next);
	step++;
	SampleStepPoint();
}


double KinematicFitter::DoKinematicFit(bool Do = true){
	int Cnt = 0;
	if(!Do){
		Finalize();
		return -1;
	}
	while(1){
		ProcessStep();
		double Chi2 = Chi2s.at(step);
		double Chi2_prev = Chi2s.at(step-1);
		if(Cnt > 3){
			break;
		}
		if(abs(Chi2_prev - Chi2 )< 0.1 and step > 1){
			Cnt++;
		}
		else{
			Cnt= 0;
		}
		if(step > MaxStep){
			break;
		}
	}
	Finalize();
	return Best_Chi2;
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
