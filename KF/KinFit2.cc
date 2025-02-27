#include "../include/KinFit2.hh"
#ifndef KinFit_cc
#define KinFit_cc
#define Debug 0
#define ShowChi2 0
void
KinematicFitter::SetVariance(double* var){
	double variance[200];
	double varianceInv[200];
	TMatrixD ScaleUp(nMeas,nMeas);
	TMatrixD ScaleDn(nMeas,nMeas);
	for(int i = 0;i<nMeas*nMeas;++i){
		variance[i]=0;
		varianceInv[i]=0;
	}
	for(int i = 0;i<nMeas;++i){
		variance[i+i*nMeas] = var[i];
		varianceInv[i+i*nMeas] = 1./var[i];
		if(ScaleParams){
			ScaleUp(i,i) = sqrt(1./var[i]);
			ScaleDn(i,i) = sqrt(var[i]);
			variance[i+i*nMeas] = 1.;
			varianceInv[i+i*nMeas] = 1.;
		}
		else{
			ScaleUp(i,i) = 1.; 
			ScaleDn(i,i) = 1.; 
		}
	}

	TMatrixD Variance(nMeas,nMeas,variance);
	TMatrixD VarianceInv(nMeas,nMeas,varianceInv);
	TMatrixD ZeroMat(nMeas,nMeas);
	Variancies.push_back(Variance);
	VarianceInvs.push_back(VarianceInv);
	dVMats.push_back(ZeroMat);
	ScalingMats.push_back(ScaleUp);
	ScalingMats.push_back(ScaleDn);
};
void
KinematicFitter::AddDiagonals(TMatrixD Cov){
	if(ScaleParams){
		auto ScaleUp = ScalingMats.at(0);
		Cov = Cov * ScaleUp;
		Cov = ScaleUp * Cov;
	}
	Variancies.at(0)+= Cov;
//	cout<<"CovarianceMat : ";
}
void KinematicFitter::ProcessStep(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step); 
	auto Meas0 = Measurements.at(0);
	auto Unkn0 = Unknowns.at(0);
	auto ScaleUp = ScalingMats.at(0);
	auto ScaleDn = ScalingMats.at(1);
	auto MS0 = ScaleUp * Meas0;// MS = M/sigma 
	auto MS = ScaleUp * Meas;
	/*
	if(UpdateVariancies and step ){
		Meas0 = Measurements.at(step-1);	
		Unkn0 = Unknowns.at(step-1);	
	}
	*/
	auto VMat = Variancies.at(step);
	auto VInv = VMat;
	VInv.SetTol(1e-26);
	VInv.Invert();
	SetConstraints();

	TMatrixD FMat = FMats.at(step);
	TMatrixD dFdM = dFdMs.at(step);
	TMatrixD dFdMS = dFdM*ScaleDn;
#if Debug
	cout<<"dFdM";
	dFdM.Print();
	cout<<"dFdMScaled";
	dFdMS.Print();
#endif
	TMatrixD dFdU = dFdUs.at(step);
	double det_dFdM=0,det_dFdU=0;
	auto dFdMT = TransposeMatrix(dFdMS);
	auto dFdUT = TransposeMatrix(dFdU);
	auto rMat = FMat + dFdMS*(MS0-MS);
	auto sMat =dFdMS*VMat*dFdMT;
	auto sInv = sMat;
	sInv.Invert();
	auto FuSIFu =	dFdUT*sInv*dFdU;
	FuSIFu.Invert();
	auto dU = (FuSIFu) * (dFdUT* (sInv) * rMat) ;
	dU = dU -dU - dU;
	auto Unkn_next = Unkn + dU;

	auto Lambda = (sInv* (rMat+ dFdU*dU));
	auto LambdaT = TransposeMatrix(Lambda);
		
	auto dM =  VMat*dFdMT*Lambda;
	auto Meas_next = MS0 - dM; 
	
	auto dMT = TransposeMatrix(dM);
	
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(ScaleDn*Meas_next);	
	rMats.push_back(rMat);
	sMats.push_back(sMat);
	auto GMat = dFdMT*sInv*dFdMS;
	auto HMat = dFdMT*sInv*dFdU;
	auto UMat = dFdUT*sInv*dFdU;
	UMat.Invert();
	auto HMatT = TransposeMatrix(HMat);
	auto HUH = (HMat*UMat*HMatT);
	auto dV = VMat*(GMat -(HUH) )*VMat;
	dVMats.push_back(dV);
	auto VMat_next = VMat-dV;
	auto VInv_next = VMat_next;
	VInv_next.SetTol(1e-26);
	VInv_next.Invert();
	int ip = 0;
	vector<double>Pull;
	for(int i = 0; i<nMeas;++i){
		double dm = dM(i,0); 	
		double dv = dV(i,i);
		dv = sqrt(dv);
		Pull.push_back( dm / dv);
	}
	Pulls.push_back(Pull);
	double Chi2 = (dMT* (VInv)*dM)(0,0) + 2 * (LambdaT * FMat )(0,0);
	Chi2s.push_back(Chi2);


#if Debug>1
	TString StepIndi = Form("[Step::%d]",step);
	cout<<StepIndi<<"##########Variance Matrix##########"<<endl;
	cout<<"V Mat : Determinant = "<<VMat.Determinant();
	VMat.Print();
	cout<<"V Mat_next : Determinant = "<<VMat_next.Determinant();
	VMat_next.Print();
	cout<<"V*VInv : Determinant = "<<(VMat*VInv).Determinant();
	(VMat*VInv).Print();
	cout<<"V*VInv_next : Determinant = "<<(VMat_next*VInv_next).Determinant();
	(VMat_next*VInv_next).Print();
	cout<<"dV : Determinant "<<dV.Determinant();
	dV.Print();
	TCanvas*canv = new TCanvas("canv","canv",600,600);
	gStyle->SetOptStat(0);
	TH2D* HdV = new TH2D("HDV","HDV",8,0,8,8,0,8);
	for(int ir=0;ir<nMeas;++ir){
		for(int il=0;il<nMeas;++il){
			HdV->SetBinContent(ir+1,8-il,dV(ir,il));
		}
	}
	cout<<"Chi2 = "<<Chi2<<endl;
	for(auto pul:Pull){
		cout<<Form("Pull %d = %g ,",ip,pul);
		ip++;
	}
	cout<<endl;
		

	cout<<StepIndi<<"##########Constraint Matrix##########"<<endl;
	cout<<"F Mat";
	FMat.Print();
	cout<<"dFdM Mat";
	dFdMS.Print();
	cout<<"dFdU Mat";
	dFdU.Print();
	cout<<"r Mat";
	rMat.Print();
	cout<<"s Mat, Det = "<<sMat.Determinant();
	sMat.Print();
	cout<<"s Inv";
	sInv.Print();
	

	cout<<StepIndi<<"##########Updated Matrix##########"<<endl;
	cout<<"dU Mat";
	dU.Print();
	cout<<"UnknNext Mat";
	Unkn_next.Print();
	cout<<"dM Mat";
	dM.Print();
	cout<<"MeasNext Mat";
	Meas_next.Print();
	cout<<"Lambda Mat";
	Lambda.Print();
	cout<<"GMat ";
	GMat.Print();
	cout<<"HMat";
	HMat.Print();
	cout<<"UMat";
	UMat.Print();
	cout<<"HUHTMat";
	(HMat*UMat*HMatT).Print();
	
	HdV->Draw("colz");
	gPad->SetMargin(0.1,0.2,0.1,0.1);
	canv->Modified();
	canv->Update();
	gSystem->ProcessEvents();
	cin.ignore();
#endif
	if(UpdateVariancies and step < 1){
		if(abs((VMat_next*VInv_next).Determinant()-1)>0.1){
			VMat_next = VMat;
			VInv_next = VInv;
		}
	}else{
		VMat_next = VMat;
		VInv_next = VInv;
	}
	Variancies.push_back(VMat_next);
	VarianceInvs.push_back(VInv_next);
	SampleStepPoint(step);
	step++;
}
void KinematicFitter::Finalize(){

	for(int is=0;is<step;++is){
		double chi2 = Chi2s.at(is);
#if ShowChi2
		cout<<Form("Step %d, chisqr %g",is,chi2)<<endl;
#endif
		if(chi2<Best_Chi2 and chi2>0){
			Best_Chi2 = chi2;
			best_step = is;
			best_pull = Pulls.at(is);
			auto BestF = FMats.at(is);
			for(int i=0;i<nConst;++i){
				best_constraints.push_back(BestF(i,0));
			}
		}
	}
	if(best_step == 0){
		best_step = step-1;
		Best_Chi2 = Chi2s.at(step-1);
		best_pull = Pulls.at(step-1); 
		auto BestF = FMats.at(step-1);
		for(int i=0;i<nConst;++i){
			best_constraints.push_back(BestF(i,0));
		}
	}
	auto InitialF = FMats.at(0);
	for(int i=0;i<nConst;++i){
		initial_constraints.push_back(InitialF(i,0));
	}
	SampleStepPoint(best_step);
	//Restore LVs from fitted parameters//	
}


double KinematicFitter::DoKinematicFit(bool Do = true){
	int Cnt = 0;
	if(!Do){
		Finalize();
		return -1;
	}
	while(1){
		ProcessStep();
#if Debug
		cout<<"Processing step "<<step<<endl;
#endif
		double Chi2 = Chi2s.at(step);
		double Chi2_prev = Chi2s.at(step-1);
		if(Cnt > 3){
			break;
		}
		if(step and Chi2<0) break;
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
	step = 0;
	best_step = 0;
	Best_Chi2 = 1e18;
	Variancies.clear();
	VarianceInvs.clear();
	Measurements.clear();
	Unknowns.clear();
	Pulls.clear();
	Chi2s.clear();
	dVMats.clear();
	FMats.clear();
	dFdMs.clear();
	dFdUs.clear();
	rMats.clear();
	sMats.clear();
	best_constraints.clear();
	initial_constraints.clear();
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
void
KinematicFitter::RotateVariance(TMatrixD J){
	auto VMat = Variancies.at(0);
	Variancies.clear();
	VMat = TransposeMatrix(J)*VMat*J;
	Variancies.push_back(VMat);
}







#endif
