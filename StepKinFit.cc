#include "src/FourVectorFitter.cc"
#include "TestKinfit.hh"
int nev = 10000;
void StepKinFit(){
	gStyle->SetOptFit(11);
	double MomScale = 1.;	
	
	double pK18 = 1.8*MomScale;
	mK*=MomScale;
	mXi*=MomScale;
	mL*=MomScale;
	mP*=MomScale;
	mPi*=MomScale;
	TLorentzVector KM(0,0,pK18,hypot(mK,pK18));
	TLorentzVector PT(0,0,0,mP);
	TLorentzVector Vertex = KM + PT;
	double VertMass[2] = {mXi,mK};
	double XiDecayMass[2] = {mL,mPi};
	double LdDecayMass[2] = {mP,mPi};
	TGenPhaseSpace EvVert;
	TGenPhaseSpace EvXi;
	TGenPhaseSpace EvLd;
	for(int i=0;i<nev;++i){
		cout<<Form("Event %d",i)<<endl;
		EvVert.SetDecay(Vertex,2,VertMass);
		EvVert.Generate();
		iev = i;
		auto Xi = *EvVert.GetDecay(0);
		auto KP = *EvVert.GetDecay(1);
		EvXi.SetDecay(Xi,2,XiDecayMass);
		EvXi.Generate();
		auto Ld = *EvXi.GetDecay(0);
		auto Pi2 = *EvXi.GetDecay(1);
		EvLd.SetDecay(Ld,2,LdDecayMass);
		EvLd.Generate();
		auto P = *EvLd.GetDecay(0);
		auto Pi1 = *EvLd.GetDecay(1);
		
		auto TVLd = Ld.Vect();
		double ThLd = Ld.Theta();
		double PhLd = Ld.Phi();

		auto TVP = P.Vect();
		auto TVPi1 = Pi1.Vect();
		auto TVPi2 = Pi2.Vect();
		PP = TVP.Mag();
		ThP = TVP.Theta();
		PhP = TVP.Phi();
		PPi1 = TVPi1.Mag();
		ThPi1 = TVPi1.Theta();
		PhPi1 = TVPi1.Phi();
		PPi2 = TVPi2.Mag();
		ThPi2 = TVPi2.Theta();
		PhPi2 = TVPi2.Phi();

		ThLdMeas = gRandom->Gaus(ThLd,ResThV);
		PhLdMeas = gRandom->Gaus(PhLd,ResPhV);
		TVector3 V1 = TVLd;
		V1.SetTheta(ThLdMeas);
		V1.SetPhi(PhLdMeas);
		TVector3 V2(0,0,0);
		PPMeas = gRandom->Gaus(PP,PP*ResP);
		ThPMeas = gRandom->Gaus(ThP,ResTh);
		PhPMeas = gRandom->Gaus(PhP,ResPh);
		PPi1Meas = gRandom->Gaus(PPi1,PPi1*ResPi1);
		ThPi1Meas = gRandom->Gaus(ThPi1,ResTh);
		PhPi1Meas = gRandom->Gaus(PhPi1,ResPh);
		PPi2Meas = gRandom->Gaus(PPi2,PPi2*ResPi2);
		ThPi2Meas = gRandom->Gaus(ThPi2,ResTh);
		PhPi2Meas = gRandom->Gaus(PhPi2,ResPh);
		
		TVector3 TVPMeas(0,0,PPMeas);
		TVPMeas.SetTheta(ThPMeas);
		TVPMeas.SetPhi(PhPMeas);
		TVector3 TVPi1Meas(0,0,PPi1Meas);
		TVPi1Meas.SetTheta(ThPi1Meas);
		TVPi1Meas.SetPhi(PhPi1Meas);
		TVector3 TVPi2Meas(0,0,PPi2Meas);
		TVPi2Meas.SetTheta(ThPi2Meas);
		TVPi2Meas.SetPhi(PhPi2Meas);
		TLorentzVector PMeas(TVPMeas,hypot(PPMeas,mP));
		TLorentzVector Pi1Meas(TVPi1Meas,hypot(PPi1Meas,mPi));
		TLorentzVector Pi2Meas(TVPi2Meas,hypot(PPi2Meas,mPi));
		double rp = PPMeas*ResP;
		double rpi1 = PPi1Meas*ResPi1;
		double rpi2 = PPi2Meas*ResPi2;

		auto LdRecon = PMeas + Pi1Meas;
		auto LdReconMeas = LdRecon;
		LdReconMeas.SetTheta(ThLdMeas);
		LdReconMeas.SetPhi(PhLdMeas);

		auto XiRecon = LdRecon + Pi2Meas;
		InvMLd = LdRecon.Mag();
		InvMXi = XiRecon.Mag();
		

		FourVectorFitter KFLd(PMeas,Pi1Meas,LdReconMeas);
		KFLd.UseVertex(true,V1,V2);
		double Variance[8] = {ResThV*ResThV,ResPhV*ResPhV,rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};

/*
		KinematicFitter KFLd(PMeas,Pi1Meas,LdRecon);
		double Variance[8] = {rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh,0,0};
*/		
		KFLd.SetVariance(Variance);
		KFLd.SetInvMass(mL);
		KFLd.SetMaximumStep(300);
		Chi2 = KFLd.DoKinematicFit();
		for(int k=0;k<100;++k){
			cout<<"..."<<endl;
		}
		Chi2 = KFLd.DoSecondaryFit();
		NStep = KFLd.GetNStep();
		BestStep = KFLd.GetBestStep();
		auto stepChi2 = KFLd.GetStepChi2();
		auto stepMassDiff = KFLd.GetStepMassDiff();
		for(int i=0;i<200;++i){
			StepChi2[i]=0;
			StepMassDiff[i]=0;
			Step[i]=-1;
		}
		for(int i=0;i<NStep;i++){
			Step[i]=i;
			StepChi2[i]=stepChi2.at(i);
			StepMassDiff[i]=stepMassDiff.at(i);
		}
		auto cont = KFLd.GetFittedLV();
		auto PCor = cont.at(0);
		auto Pi1Cor = cont.at(1);
		auto LdFit = PCor+Pi1Cor;
		InvMLdCor = LdFit.Mag();	
		auto TVPCor = PCor.Vect();
		PPCor = TVPCor.Mag(); ThPCor = TVPCor.Theta();
		PhPCor = TVPCor.Phi();
		
		auto TVPi1Cor = Pi1Cor.Vect();
		PPi1Cor = TVPi1Cor.Mag();
		ThPi1Cor = TVPi1Cor.Theta();
		PhPi1Cor = TVPi1Cor.Phi();
		auto XiFit = LdFit + Pi2Meas;
		InvMXiCor = XiFit.Mag();	
		cout<<Form("MassLd = %g,MassXi=%g",InvMLdCor,InvMXiCor)<<endl;
		cin.ignore();
	}
}
