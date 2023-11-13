#include "src/KinFit2.cc"
#include "TestKinfit.hh"
double mXi = 1.321;
double mL = 1.115;
double mP = 0.938;
double mK = 0.493;
double mPi = 0.139;
int nev = 10000;
void GenerateKinFit(){
	gStyle->SetOptFit(11);
	double pK18 = 1.8;
	TLorentzVector KM(0,0,pK18,hypot(mK,pK18));
	TLorentzVector PT(0,0,0,mP);
	TLorentzVector Vertex = KM + PT;
	double VertMass[2] = {mXi,mK};
	double XiDecayMass[2] = {mL,mPi};
	double LdDecayMass[2] = {mP,mPi};
	TGenPhaseSpace EvVert;
	TGenPhaseSpace EvXi;
	TGenPhaseSpace EvLd;
	TFile* file = new TFile("KinFitTest.root","recreate");
	TTree* tree = new TTree("tree","tree");
	SetBranches(tree);
	TH1D* HistLd = new TH1D("InvMLd","InvMLd",100,1.08,1.14);
	TH1D* HistXi = new TH1D("InvMXi","InvMXi",100,1.25,1.4);
	TH1D* HistLdFit = new TH1D("InvMLdFit","InvMLdFit",1000,1.08,1.14);
	TH1D* HistXiFit = new TH1D("InvMXiFit","InvMXiFit",100,1.25,1.4);
	TH1D* HistLdP = new TH1D("HistLdP","HistLdP",100,-0.3,0.3);
	TH1D* HistLdTh = new TH1D("HistLdTh","HistLdTh",100,-0.3,0.3);
	TH1D* HistLdPh = new TH1D("HistLdPh","HistLdPh",100,-0.3,0.3);
	TH1D* HistPP = new TH1D("HistPP","HistPP",100,-0.3,0.3);
	TH1D* HistPTh = new TH1D("HistPTh","HistPTh",100,-0.3,0.3);
	TH1D* HistPPh = new TH1D("HistPPh","HistPPh",100,-0.3,0.3);
	TH1D* HistPi1P = new TH1D("HistPi1P","HistPi1P",100,-0.1,0.1);
	TH1D* HistPi1Th = new TH1D("HistPi1Th","HistPi1Th",100,-0.3,0.3);
	TH1D* HistPi1Ph = new TH1D("HistPi1Ph","HistPi1Ph",100,-0.3,0.3);
	TH1D* HistXiP = new TH1D("HistXiP","HistXiP",100,-0.3,0.3);
	TH1D* HistXiTh = new TH1D("HistXiTh","HistXiTh",100,-0.3,0.3);
	TH1D* HistXiPh = new TH1D("HistXiPh","HistXiPh",100,-0.3,0.3);
	TH1D* HistLdCorP = new TH1D("HistLdCorP","HistLdCorP",100,-0.3,0.3);
	TH1D* HistLdCorTh = new TH1D("HistLdCorTh","HistLdCorTh",100,-0.3,0.3);
	TH1D* HistLdCorPh = new TH1D("HistLdCorPh","HistLdCorPh",100,-0.3,0.3);
	
	TH1D* HistPi2P = new TH1D("HistPi2P","HistPi2P",100,-0.1,0.1);
	TH1D* HistPi2Th = new TH1D("HistPi2Th","HistPi2Th",100,-0.3,0.3);
	TH1D* HistPi2Ph = new TH1D("HistPi2Ph","HistPi2Ph",100,-0.3,0.3);
	
	TH1D* HistDLdP = new TH1D("HistDLdP","HistDLdP",100,-0.3,0.3);
	TH1D* HistDLdTh = new TH1D("HistDLdTh","HistDLdTh",100,-0.3,0.3);
	TH1D* HistDLdPh = new TH1D("HistDLdPh","HistDLdPh",100,-0.3,0.3);
	TH1D* HistDPP = new TH1D("HistDPP","HistDPP",100,-0.3,0.3);
	TH1D* HistDPTh = new TH1D("HistDPTh","HistDPTh",100,-0.3,0.3);
	TH1D* HistDPPh = new TH1D("HistDPPh","HistDPPh",100,-0.3,0.3);
	TH1D* HistDPi1P = new TH1D("HistDPi1P","HistDPi1P",100,-0.1,0.1);
	TH1D* HistDPi1Th = new TH1D("HistDPi1Th","HistDPi1Th",100,-0.3,0.3);
	TH1D* HistDPi1Ph = new TH1D("HistDPi1Ph","HistDPi1Ph",100,-0.3,0.3);

	TH1D* HistDXiP = new TH1D("HistDXiP","HistDXiP",100,-0.3,0.3);
	TH1D* HistDXiTh = new TH1D("HistDXiTh","HistDXiTh",100,-0.3,0.3);
	TH1D* HistDXiPh = new TH1D("HistDXiPh","HistDXiPh",100,-0.3,0.3);
	TH1D* HistDLdCorP = new TH1D("HistDLdCorP","HistDLdCorP",100,-0.3,0.3);
	TH1D* HistDLdCorTh = new TH1D("HistDLdCorTh","HistDLdCorTh",100,-0.3,0.3);
	TH1D* HistDLdCorPh = new TH1D("HistDLdCorPh","HistDLdCorPh",100,-0.3,0.3);
	TH1D* HistDPi2P = new TH1D("HistDPi2P","HistDPi2P",100,-0.1,0.1);
	TH1D* HistDPi2Th = new TH1D("HistDPi2Th","HistDPi2Th",100,-0.3,0.3);
	TH1D* HistDPi2Ph = new TH1D("HistDPi2Ph","HistDPi2Ph",100,-0.3,0.3);



	TH1D* HistDifLdP = new TH1D("HistDifLdP","HistDifLdP",100,-0.1,0.1);
	TH1D* HistDifLdTh = new TH1D("HistDifLdTh","HistDifLdTh",100,-0.3,0.3);
	TH1D* HistDifLdPh = new TH1D("HistDifLdPh","HistDifLdPh",100,-0.3,0.3);
	TH1D* HistDifPP = new TH1D("HistDifPP","HistDifPP",100,-0.3,0.3);
	TH1D* HistDifPTh = new TH1D("HistDifPTh","HistDifPTh",100,-0.3,0.3);
	TH1D* HistDifPPh = new TH1D("HistDifPPh","HistDifPPh",100,-0.3,0.3);
	TH1D* HistDifPi1P = new TH1D("HistDifPi1P","HistDifPi1P",100,-0.1,0.1);
	TH1D* HistDifPi1Th = new TH1D("HistDifPi1Th","HistDifPi1Th",100,-0.3,0.3);
	TH1D* HistDifPi1Ph = new TH1D("HistDifPi1Ph","HistDifPi1Ph",100,-0.3,0.3);

	TH1D* HistChi2 = new TH1D("Chi2","Chi2",500,0,50);
	TH2D* HistChi22D = new TH2D("Chi2:InvM","Chi2:InvM",100,0,20,100,1.1,1.14);
	TH1D* HistProb = new TH1D("Prob","Prob",100,0,1.1);
	TH2D* HistProb2D = new TH2D("Prob:InvM","Prob:InvM",100,0,1.1,100,1.1,1.14);
	TH1D* HistMM = new TH1D("MM","MM",100,-0.3,0.3);
	TH1D* HistMMCor = new TH1D("MMCor","MMCor",100,-0.3,0.3);
	for(int i=0;i<nev;++i){
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
		PLd = TVLd.Mag();
		ThLd = Ld.Theta();
		PhLd = Ld.Phi();

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
		auto TVXiMeas = XiRecon.Vect();
		auto TVXi = Xi.Vect();
		PXiMeas = TVXiMeas.Mag();
		ThXiMeas = TVXiMeas.Theta();
		PhXiMeas = TVXiMeas.Phi();
		PXi = TVXi.Mag();
		ThXi = TVXi.Theta();
		PhXi = TVXi.Phi();
		HistLd->Fill(InvMLd);
		HistXi->Fill(InvMXi);
	
		KinematicFitter KFLd(PMeas,Pi1Meas,LdReconMeas);
		KFLd.UseVertex(true,V1,V2);
		double Variance[8] = {ResThV*ResThV,ResPhV*ResPhV,rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};
/*

		KinematicFitter KFLd(PMeas,Pi1Meas,LdRecon);
		double Variance[8] = {rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh,0,0};
*/	
		KFLd.SetVariance(Variance);
		KFLd.SetInvMass(mL);
		KFLd.SetMaximumStep(300);
		Chi2=		KFLd.DoKinematicFit();
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
		HistChi2->Fill(Chi2);
		auto cont = KFLd.GetFittedLV();
		auto PCor = cont.at(0);
		auto Pi1Cor = cont.at(1);
		auto LdCor = cont.at(2);
		auto TVLdCor = LdCor.Vect();
		PLdCor = TVLdCor.Mag();
		ThLdCor = TVLdCor.Theta();
		PhLdCor = TVLdCor.Phi();
		auto LdFit = PCor+Pi1Cor;
		auto TVLdFit = LdFit.Vect();
		auto TVLdMeas = TVPMeas + TVPi1Meas;
		PLdMeas = TVLdMeas.Mag();
		InvMLdCor = LdFit.Mag();	
		auto TVPCor = PCor.Vect();
		PPCor = TVPCor.Mag();
		ThPCor = TVPCor.Theta();
		PhPCor = TVPCor.Phi();
		
		auto TVPi1Cor = Pi1Cor.Vect();
		PPi1Cor = TVPi1Cor.Mag();
		ThPi1Cor = TVPi1Cor.Theta();
		PhPi1Cor = TVPi1Cor.Phi();
		auto XiFit = LdCor + Pi2Meas;
		KinematicFitter KFXi(LdCor,Pi2Meas,XiFit);
		double rL = ResPLd*PLdCor;
		double VarianceXi[8] = {rp*rp,ResThLd*ResThLd,ResPhLd*ResPhLd,rpi2*rpi2,ResTh*ResTh,ResPh*ResPh,0,0};
		KFXi.SetVariance(VarianceXi);
		KFXi.SetInvMass(mXi);	
		KFXi.SetMaximumStep(300);
		Chi2Xi=		KFXi.DoKinematicFit(false);
		auto contXi = KFXi.GetFittedLV();
		auto LdCorCor = contXi.at(0);
		auto Pi2Cor = contXi.at(1);
		auto XiCor = LdCorCor + Pi2Cor;
		auto TVPi2Cor = Pi2Cor.Vect();
		PPi2Cor = TVPi2Cor.Mag();
		ThPi2Cor = TVPi2Cor.Theta();
		PhPi2Cor = TVPi2Cor.Phi();
		auto TVLdCorCor = LdCorCor.Vect();
		PLdCorCor = TVLdCorCor.Mag();
		ThLdCorCor = TVLdCorCor.Theta();
		PhLdCorCor = TVLdCorCor.Phi();
		auto TVXiCor = XiCor.Vect();
		PXiCor = TVXiCor.Mag();
		ThXiCor = TVXiCor.Theta();
		PhXiCor = TVXiCor.Phi();
		//		auto XiFit = LdFit + Pi2Meas;
//		auto XiFit = LdCor + Pi2Meas;
		InvMXiCor = XiFit.Mag();	
		HistDLdP->Fill(PLd-PLdMeas);
		HistDLdTh->Fill(ThLd-ThLdMeas);
		HistDLdPh->Fill(PhLd - PhLdMeas);
		HistDPP->Fill(PP-PPMeas);
		HistDPTh->Fill(ThP-ThPMeas);
		HistDPPh->Fill(PhP - PhPMeas);
		HistDPi1P->Fill(PPi1-PPi1Meas);
		HistDPi1Th->Fill(ThPi1-ThPi1Meas);
		HistDPi1Ph->Fill(PhPi1 - PhPi1Meas);
		

		HistDXiP->Fill(PXi-PXiMeas);
		HistDXiTh->Fill(ThXi-ThXiMeas);
		HistDXiPh->Fill(PhXi - PhXiMeas);
		HistDLdCorP->Fill(PLd-PLdCor);
		HistDLdCorTh->Fill(ThLd-ThLdCor);
		HistDLdCorPh->Fill(PhLd - PhLdCor);
		HistDPi2P->Fill(PPi2-PPi2Meas);
		HistDPi2Th->Fill(ThPi2-ThPi2Meas);
		HistDPi2Ph->Fill(PhPi2 - PhPi2Meas);
		
		
		
		HistDifPP->Fill(PPCor-PPMeas);
		HistDifPTh->Fill(ThPCor-ThPMeas);
		HistDifPPh->Fill(PhPCor - PhPMeas);
		HistDifPi1P->Fill(PPi1Cor-PPi1Meas);
		HistDifPi1Th->Fill(ThPi1Cor-ThPi1Meas);
		HistDifPi1Ph->Fill(PhPi1Cor - PhPi1Meas);
		
		HistLdFit->Fill(InvMLdCor);
		HistXiFit->Fill(InvMXiCor);
		HistLdP->Fill(PLd-PLdCor);
		HistLdTh->Fill(ThLd-ThLdCor);
		HistLdPh->Fill(PhLd - PhLdCor);
		HistPP->Fill(PP-PPCor);
		HistPTh->Fill(ThP-ThPCor);
		HistPPh->Fill(PhP - PhPCor);
		HistPi1P->Fill(PPi1-PPi1Cor);
		HistPi1Th->Fill(ThPi1-ThPi1Cor);
		HistPi1Ph->Fill(PhPi1 - PhPi1Cor);


		HistXiP->Fill(PXi-PXiCor);
		HistXiTh->Fill(ThXi-ThXiCor);
		HistXiPh->Fill(PhXi - PhXiCor);
		HistLdCorP->Fill(PLd-PLdCorCor);
		HistLdCorTh->Fill(ThLd-ThLdCorCor);
		HistLdCorPh->Fill(PhLd - PhLdCorCor);
		HistPi2P->Fill(PPi2-PPi2Cor);
		HistPi2Th->Fill(ThPi2-ThPi2Cor);
		HistPi2Ph->Fill(PhPi2 - PhPi2Cor);



		double prob = ROOT::Math::chisquared_cdf(Chi2,KFLd.GetNDF());
		HistProb->Fill(prob);
		HistProb2D->Fill(prob,LdFit.Mag());
		HistChi22D->Fill(Chi2,LdFit.Mag());
		auto MM = Xi - XiRecon;
		auto MMCor = Xi - XiCor;
		HistMM->Fill(MM.Mag());
		HistMMCor->Fill(MMCor.Mag());
		tree->Fill();	
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	c1->cd(1);
	HistLd->Draw();
	HistLd->Fit("gaus","Q");
	c1->cd(2);
	HistXi->Draw();
	HistXi->Fit("gaus","Q");
	c1->cd(3);
	HistLdFit->Draw();
	HistLdFit->Fit("gaus","Q");
	c1->cd(4);
	HistXiFit->Draw();
	HistXiFit->Fit("gaus","Q");
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(3,3);
	c2->cd(1);
	HistLdP->Draw();
	HistLdP->SetLineColor(kGreen);
	HistLdP->Fit("gaus","Q");
	c2->cd(2);
	HistLdTh->Draw();
	HistLdTh->SetLineColor(kGreen);
	HistLdTh->Fit("gaus","Q");
	c2->cd(3);
	HistLdPh->Draw();
	HistLdPh->SetLineColor(kGreen);
	HistLdPh->Fit("gaus","Q");
	c2->cd(4);
	HistPP->Draw();
	HistPP->SetLineColor(kGreen);
	HistPP->Fit("gaus","Q");
	c2->cd(5);
	HistPTh->Draw();
	HistPTh->SetLineColor(kGreen);
	HistPTh->Fit("gaus","Q");
	c2->cd(6);
	HistPPh->Draw();
	HistPPh->SetLineColor(kGreen);
	HistPPh->Fit("gaus","Q");
	c2->cd(7);
	HistPi1P->Draw();
	HistPi1P->SetLineColor(kGreen);
	HistPi1P->Fit("gaus","Q");
	c2->cd(8);
	HistPi1Th->Draw();
	HistPi1Th->SetLineColor(kGreen);
	HistPi1Th->Fit("gaus","Q");
	c2->cd(9);
	HistPi1Ph->Draw();
	HistPi1Ph->SetLineColor(kGreen);
	HistPi1Ph->Fit("gaus","Q");
	TCanvas* c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(3,3);
	c3->cd(1);
	HistDLdP->Draw();
	HistDLdP->Fit("gaus","Q");
	c3->cd(2);
	HistDLdTh->Draw();
	HistDLdTh->Fit("gaus","Q");
	c3->cd(3);
	HistDLdPh->Draw();
	HistDLdPh->Fit("gaus","Q");
	c3->cd(4);
	HistDPP->Draw();
	HistDPP->Fit("gaus","Q");
	c3->cd(5);
	HistDPTh->Draw();
	HistDPTh->Fit("gaus","Q");
	c3->cd(6);
	HistDPPh->Draw();
	HistDPPh->Fit("gaus","Q");
	c3->cd(7);
	HistDPi1P->Draw();
	HistDPi1P->Fit("gaus","Q");
	c3->cd(8);
	HistDPi1Th->Draw();
	HistDPi1Th->Fit("gaus","Q");
	c3->cd(9);
	HistDPi1Ph->Draw();
	HistDPi1Ph->Fit("gaus","Q");
	TCanvas* c4 = new TCanvas("c4","c4",1200,600);
	c4->Divide(3,3);
	c4->cd(1);
	HistDifPP->Draw();
	HistDifPP->Fit("gaus","Q");
	c4->cd(2);
	HistDifPTh->Draw();
	HistDifPTh->Fit("gaus","Q");
	c4->cd(3);
	HistDifPPh->Draw();
	HistDifPPh->Fit("gaus","Q");
	c4->cd(4);
	HistDifPi1P->Draw();
	HistDifPi1P->Fit("gaus","Q");
	c4->cd(5);
	HistDifPi1Th->Draw();
	HistDifPi1Th->Fit("gaus","Q");
	c4->cd(6);
	HistDifPi1Ph->Draw();
	HistDifPi1Ph->Fit("gaus","Q");
	

	TCanvas* c5 = new TCanvas("c5","c5",600,600);
	c5->Divide(2,1);
	c5->cd(1);
	HistChi2->Draw();
	c5->cd(2);
	HistChi22D->Draw("colz");
	file->Write();
	TCanvas* c6 = new TCanvas("c6","c6",1200,600);
	c6->Divide(3,3);
	c6->cd(1);
	HistLdP->Draw();
	HistDLdP->Draw("same");
	c6->cd(2);
	HistLdTh->Draw();
	HistDLdTh->Draw("same");
	c6->cd(3);
	HistLdPh->Draw();
	HistDLdPh->Draw("same");
	c6->cd(4);
	HistPP->Draw();
	HistDPP->Draw("same");
	c6->cd(5);
	HistPTh->Draw();
	HistDPTh->Draw("same");
	c6->cd(6);
	HistPPh->Draw();
	HistDPPh->Draw("same");
	c6->cd(7);
	HistPi1P->Draw();
	HistDPi1P->Draw("same");
	c6->cd(8);
	HistPi1Th->Draw();
	HistDPi1Th->Draw("same");
	c6->cd(9);
	HistPi1Ph->Draw();
	HistDPi1Ph->Draw("same");
	TCanvas* c7 = new TCanvas("c7","c7",600,600);
	c7->Divide(2,1);
	c7->cd(1);
	HistProb->Draw();
	c7->cd(2);
	HistProb2D->Draw("colz");

	TCanvas* c8 = new TCanvas("c8","c8",1200,600);
	c8->Divide(3,3);
	c8->cd(1);
	HistXiP->Draw();
	HistXiP->SetLineColor(kGreen);
	HistDXiP->Draw("same");
	c8->cd(2);
	HistXiTh->Draw();
	HistXiTh->SetLineColor(kGreen);
	HistDXiTh->Draw("same");
	c8->cd(3);
	HistXiPh->Draw();
	HistXiPh->SetLineColor(kGreen);
	HistDXiPh->Draw("same");
	c8->cd(4);
	HistLdCorP->Draw();
	HistLdCorP->SetLineColor(kGreen);
	HistDLdCorP->Draw("same");
	c8->cd(5);
	HistLdCorTh->Draw();
	HistLdCorTh->SetLineColor(kGreen);
	HistDLdCorTh->Draw("same");
	c8->cd(6);
	HistLdCorPh->Draw();
	HistLdCorPh->SetLineColor(kGreen);
	HistDLdCorPh->Draw("same");
	c8->cd(7);
	HistPi2P->Draw();
	HistPi2P->SetLineColor(kGreen);
	HistDPi2P->Draw("same");
	c8->cd(8);
	HistPi2Th->Draw();
	HistPi2Th->SetLineColor(kGreen);
	HistDPi2Th->Draw("same");
	c8->cd(9);
	HistPi2Ph->Draw();
	HistPi2Ph->SetLineColor(kGreen);
	HistDPi2Ph->Draw("same");
	TCanvas* c9 = new TCanvas("c9","c9",600,600);
	HistMM->Draw();
	HistMMCor->SetLineColor(kGreen);
	HistMMCor->Draw("same");
}
void TestKinFit(){
	gStyle->SetOptFit(11);
	double pK18 = 1.8;
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
		

		KinematicFitter KFLd(PMeas,Pi1Meas,LdReconMeas);
		KFLd.UseVertex(true,V1,V2);
		double Variance[8] = {ResThV*ResThV,ResPhV*ResPhV,rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};

/*
		KinematicFitter KFLd(PMeas,Pi1Meas,LdRecon);
		double Variance[8] = {rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh,0,0};
*/		
		KFLd.SetVariance(Variance);
		KFLd.SetInvMass(mL);
		KFLd.SetMaximumStep(300);
		Chi2=		KFLd.DoKinematicFit();
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
