#include "src/KinFit2.cc"
double mXi = 1.321;
double mL = 1.115;
double mP = 0.938;
double mK = 0.493;
double mPi = 0.139;
int nev = 10000;
void TestKinfit(){
	double pK18 = 1.8;
	TLorentzVector KM(0,0,pK18,hypot(mK,pK18));
	TLorentzVector PT(0,0,0,mP);
	TLorentzVector Vertex = KM + PT;
	cout<<Vertex.Mag()<<endl;
	double VertMass[2] = {mXi,mK};
	double XiDecayMass[2] = {mL,mPi};
	double LdDecayMass[2] = {mP,mPi};
	TGenPhaseSpace EvVert;
	TGenPhaseSpace EvXi;
	TGenPhaseSpace EvLd;
	TH1D* HistLd = new TH1D("InvMLd","InvMLd",100,1.04,1.24);
	TH1D* HistXi = new TH1D("InvMXi","InvMXi",100,1.2,1.45);
	TH1D* HistLdFit = new TH1D("InvMLdFit","InvMLdFit",100,1.10,1.14);
	TH1D* HistXiFit = new TH1D("InvMXiFit","InvMXiFit",100,1.2,1.45);
	TH1D* HistPP = new TH1D("HistPP","HistPP",100,-0.3,0.3);
	TH1D* HistPTh = new TH1D("HistPTh","HistPTh",100,-1.7,1.7);
	TH1D* HistPPh = new TH1D("HistPPh","HistPPh",100,-1.7,1.7);
	TH1D* HistPi1P = new TH1D("HistPi1P","HistPi1P",100,-0.3,0.3);
	TH1D* HistPi1Th = new TH1D("HistPi1Th","HistPi1Th",100,-1.7,1.7);
	TH1D* HistPi1Ph = new TH1D("HistPi1Ph","HistPi1Ph",100,-1.7,1.7);
	
	TH1D* HistDPP = new TH1D("HistDPP","HistDPP",100,-0.3,0.3);
	TH1D* HistDPTh = new TH1D("HistDPTh","HistDPTh",100,-1.7,1.7);
	TH1D* HistDPPh = new TH1D("HistDPPh","HistDPPh",100,-1.7,1.7);
	TH1D* HistDPi1P = new TH1D("HistDPi1P","HistDPi1P",100,-0.3,0.3);
	TH1D* HistDPi1Th = new TH1D("HistDPi1Th","HistDPi1Th",100,-1.7,1.7);
	TH1D* HistDPi1Ph = new TH1D("HistDPi1Ph","HistDPi1Ph",100,-1.7,1.7);


	double ResP = 0.05;
	double ResPi1 = 0.03;
	double ResPi2 = 0.03;
	double ResTh = 0.02; 
	double ResPh = 0.03; 


	for(int i=0;i<nev;++i){
		EvVert.SetDecay(Vertex,2,VertMass);
		EvVert.Generate();
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
		
		auto TVP = P.Vect();
		auto TVPi1 = Pi1.Vect();
		auto TVPi2 = Pi2.Vect();
		double PP = TVP.Mag();
		double ThP = TVP.Theta();
		double PhP = TVP.Phi();
		double PPi1 = TVPi1.Mag();
		double ThPi1 = TVPi1.Theta();
		double PhPi1 = TVPi1.Phi();
		double PPi2 = TVPi2.Mag();
		double ThPi2 = TVPi2.Theta();
		double PhPi2 = TVPi2.Phi();

		double PPMeas = gRandom->Gaus(PP,PP*ResP);
		double ThPMeas = gRandom->Gaus(ThP,ResTh);
		double PhPMeas = gRandom->Gaus(PhP,ResPh);
		double PPi1Meas = gRandom->Gaus(PPi1,PPi1*ResPi1);
		double ThPi1Meas = gRandom->Gaus(ThPi1,ResTh);
		double PhPi1Meas = gRandom->Gaus(PhPi1,ResPh);
		double PPi2Meas = gRandom->Gaus(PPi2,PPi2*ResPi2);
		double ThPi2Meas = gRandom->Gaus(ThPi2,ResTh);
		double PhPi2Meas = gRandom->Gaus(PhPi2,ResPh);
		
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


		auto XiRecon = LdRecon + Pi2Meas;
		HistLd->Fill(LdRecon.Mag());
		HistXi->Fill(XiRecon.Mag());
		KinematicFitter KF(PMeas,Pi1Meas,LdRecon);
		double Variance[6] = {rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};
		KF.SetVariance(Variance);
		KF.SetInvMass(mL);
//		KF.SetMaximumStep(1);
		KF.DoKinematicFit();
		auto cont = KF.GetFittedLV();
		auto PCor = cont.at(0);
		auto Pi1Cor = cont.at(1);
		auto LdFit = PCor+Pi1Cor;
		
		auto TVPCor = PCor.Vect();
		auto PPCor = TVPCor.Mag();
		auto ThPCor = TVPCor.Theta();
		auto PhPCor = TVPCor.Phi();
		
		auto TVPi1Cor = Pi1Cor.Vect();
		auto PPi1Cor = TVPi1Cor.Mag();
		auto ThPi1Cor = TVPi1Cor.Theta();
		auto PhPi1Cor = TVPi1Cor.Phi();
		HistLdFit->Fill(LdFit.Mag());
		auto XiFit = LdFit + Pi2Meas;
		HistXiFit->Fill(XiFit.Mag());
		HistDPP->Fill(PP-PPMeas);
		HistDPTh->Fill(ThP-ThPMeas);
		HistDPPh->Fill(PhP - PhPMeas);
		HistDPi1P->Fill(PPi1-PPi1Meas);
		HistDPi1Th->Fill(ThPi1-ThPi1Meas);
		HistDPi1Ph->Fill(PhPi1 - PhPi1Meas);
		HistPP->Fill(PP-PPCor);
		HistPTh->Fill(ThP-ThPCor);
		HistPPh->Fill(PhP - PhPCor);
		HistPi1P->Fill(PPi1-PPi1Cor);
		HistPi1Th->Fill(ThPi1-ThPi1Cor);
		HistPi1Ph->Fill(PhPi1 - PhPi1Cor);


	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	c1->cd(1);
	HistLd->Draw();
	HistLd->Fit("gaus");
	c1->cd(2);
	HistXi->Draw();
	HistXi->Fit("gaus");
	c1->cd(3);
	HistLdFit->Draw();
	HistLdFit->Fit("gaus");
	c1->cd(4);
	HistXiFit->Draw();
	HistXiFit->Fit("gaus");
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(3,2);
	c2->cd(1);
	HistPP->Draw();
	HistPP->Fit("gaus");
	c2->cd(2);
	HistPTh->Draw();
	HistPTh->Fit("gaus");
	c2->cd(3);
	HistPPh->Draw();
	HistPPh->Fit("gaus");
	c2->cd(4);
	HistPi1P->Draw();
	HistPi1P->Fit("gaus");
	c2->cd(5);
	HistPi1Th->Draw();
	HistPi1Th->Fit("gaus");
	c2->cd(6);
	HistPi1Ph->Draw();
	HistPi1Ph->Fit("gaus");
	TCanvas* c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(3,2);
	c3->cd(1);
	HistDPP->Draw();
	HistDPP->Fit("gaus");
	c3->cd(2);
	HistDPTh->Draw();
	HistDPTh->Fit("gaus");
	c3->cd(3);
	HistDPPh->Draw();
	HistDPPh->Fit("gaus");
	c3->cd(4);
	HistDPi1P->Draw();
	HistDPi1P->Fit("gaus");
	c3->cd(5);
	HistDPi1Th->Draw();
	HistDPi1Th->Fit("gaus");
	c3->cd(6);
	HistDPi1Ph->Draw();
	HistDPi1Ph->Fit("gaus");

}
