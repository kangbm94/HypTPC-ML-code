#include "src/FourVectorXYZFitter.cc"
#include "TestKinfit.hh"
int nev = 10000;
double res_x = 0.07;
double res_y = 0.07;
double res_z = 0.07;
void TestKinfit1CXYZ(){
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


	TH2D* HistCor[81];
	for(int i=0;i<81 ;++i){
		TString title = Form("Hist%d",i);
		HistCor[i] = new TH2D(title,title,300,-10,10,300,-10,10);
	}

	for(int i=0;i<nev;++i){
		EvVert.SetDecay(Vertex,2,VertMass);
		EvVert.Generate();
		iev = i;
		if(iev%1000==0){
			cout<<Form("Event %d",iev)<<endl;
		}
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
		double PxLd = TVLd.x();
		double PyLd = TVLd.y();
		double PzLd = TVLd.z();

		auto TVP = P.Vect();
		auto TVPi1 = Pi1.Vect();
		auto TVPi2 = Pi2.Vect();
		
		PP = TVP.Mag();
		double PxP = TVP.x();
		double PyP = TVP.y();
		double PzP = TVP.z();
		PPi1 = TVPi1.Mag();
		double PxPi1 = TVPi1.x();
		double PyPi1 = TVPi1.y();
		double PzPi1 = TVPi1.z();
		PPi2 = TVPi2.Mag();
		double PxPi2 = TVPi2.x();
		double PyPi2 = TVPi2.y();
		double PzPi2 = TVPi2.z();

//		ThLdMeas = (ThLd);
//		PhLdMeas = (PhLd);
		TVector3 V1 = TVLd;
//		V1.SetTheta(ThLdMeas);
//		V1.SetPhi(PhLdMeas);
		TVector3 V2(0,0,0);
		double PxPMeas = gRandom->Gaus(PxP,PxP*res_x);
		double PyPMeas = gRandom->Gaus(PyP,PyP*res_y);
		double PzPMeas = gRandom->Gaus(PzP,PzP*res_z);
		double PxPi1Meas = gRandom->Gaus(PxPi1,PxPi1*res_x);
		double PyPi1Meas = gRandom->Gaus(PyPi1,PyPi1*res_y);
		double PzPi1Meas = gRandom->Gaus(PzPi1,PzPi1*res_z);
		double PxPi2Meas = gRandom->Gaus(PxPi2,PxPi2*res_x);
		double PyPi2Meas = gRandom->Gaus(PyPi2,PyPi2*res_y);
		double PzPi2Meas = gRandom->Gaus(PzPi2,PzPi2*res_z);


		TVector3 TVPMeas(PxPMeas,PyPMeas,PzPMeas);
		TVector3 TVPi1Meas(PxPi1Meas,PyPi1Meas,PzPi1Meas);
		TVector3 TVPi2Meas(PxPi2Meas,PyPi2Meas,PzPi2Meas);
		double PPMeas = TVPMeas.Mag();
		double PPi1Meas = TVPi1Meas.Mag();
		double PPi2Meas = TVPi2Meas.Mag();
		TLorentzVector PMeas(TVPMeas,hypot(PPMeas,mP));
		TLorentzVector Pi1Meas(TVPi1Meas,hypot(PPi1Meas,mPi));
		TLorentzVector Pi2Meas(TVPi2Meas,hypot(PPi2Meas,mPi));
		double rxp = PxPMeas*res_x;
		double ryp = PyPMeas*res_y;
		double rzp = PzPMeas*res_z;
		double rxpi1 = PxPi1Meas*res_x;
		double rypi1 = PyPi1Meas*res_y;
		double rzpi1 = PzPi1Meas*res_z;
		double rxpi2 = PxPi2Meas*res_x;
		double rypi2 = PyPi2Meas*res_y;
		double rzpi2 = PzPi2Meas*res_z;

		auto LdRecon = PMeas + Pi1Meas;
		auto LdReconMeas = LdRecon;
		




		auto XiRecon = LdRecon + Pi2Meas;
		InvMLd = LdRecon.Mag();
		InvMXi = (LdRecon+Pi2).Mag();
		auto TVXiMeas = XiRecon.Vect();
		auto TVXi = Xi.Vect();
		PXiMeas = TVXiMeas.Mag();
		PXi = TVXi.Mag();
		HistLd->Fill(InvMLd);
		HistXi->Fill(InvMXi);
	


		FourVectorXYZFitter KFLd(PMeas,Pi1Meas,LdRecon);
		double Variance[6] = {rxp*rxp,ryp*ryp,rzp*rzp,rxpi1*rxpi1,rypi1*rypi1,rzpi1*rzpi1};
	
		KFLd.SetInvMass(mL);
		KFLd.SetMaximumStep(10);
//		KFLd.UpdateVariance();
		KFLd.ScaleParameters(false);
		KFLd.SetVariance(Variance);
		Chi2=		KFLd.DoKinematicFit();
		NStep = KFLd.GetNStep();
		BestStep = KFLd.GetBestStep();
		auto stepChi2 = KFLd.GetStepChi2();
		auto PullLd = KFLd.GetPull();
		auto UPullLd = KFLd.GetUPull();
		auto UV = KFLd.GetUnmeasuredCovariance();
		auto cont = KFLd.GetFittedLV();
		auto PCor = cont.at(0);
		auto Pi1Cor = cont.at(1);
//		auto LdCor = cont.at(2);
		auto LdCor = PCor+Pi1Cor;
		auto TVLdCor = LdCor.Vect();
		PLdCor = TVLdCor.Mag();
		ThLdCor = TVLdCor.Theta();
		PhLdCor = TVLdCor.Phi();
		vector<double>tempPull;
		auto sigpx = sqrt(abs(UV(0,0)));
		auto sigpy = sqrt(abs(UV(1,1)));
		auto sigpz = sqrt(abs(UV(2,2)));
		auto dPxLd = (Ld.X() - LdCor.X())/sigpx;
		auto dPyLd = (Ld.Y() - LdCor.Y())/sigpy;
		auto dPzLd = (Ld.Z() - LdCor.Z())/sigpz;
//		auto dPLd = (PLd - PLdCor);
//		auto dThLd = (ThLd - ThLdCor);
//		auto dPhLd = (PhLd - PhLdCor);	
//		UV.Print();
		tempPull.push_back(dPxLd);
		tempPull.push_back(dPyLd);
		tempPull.push_back(dPzLd);
		for(auto p:UPullLd){
//			tempPull.push_back(p);
		}
		for(auto p:PullLd){
			tempPull.push_back(p);
		}
		PullLd = tempPull;
		for(int col=0;col<9;++col){
			for(int row=0;row<col;++row){
				if(row == col) continue;
				HistCor[9*col + row] ->Fill(PullLd.at(col),PullLd.at(row));
			}
		}
		PullPLd=PullLd.at(0);	
		PullThLd=PullLd.at(1);	
		PullPhLd=PullLd.at(2);	
		PullPP=PullLd.at(3);	
		PullThP=PullLd.at(4);	
		PullPhP=PullLd.at(5);	
		PullPPi1=PullLd.at(6);	
		PullThPi1=PullLd.at(7);	
		PullPhPi1=PullLd.at(8);	
		HistPullPLd->Fill(PullPLd);
		HistPullThLd->Fill(PullThLd);
		HistPullPhLd->Fill(PullPhLd);
		HistPullPP->Fill(PullPP);
		HistPullThP->Fill(PullThP);
		HistPullPhP->Fill(PullPhP);
		HistPullPPi1->Fill(PullPPi1);
		HistPullThPi1->Fill(PullThPi1);
		HistPullPhPi1->Fill(PullPhPi1);
		for(int i=0;i<200;++i){
			StepChi2[i]=0;
			Step[i]=-1;
		}
		for(int i=0;i<NStep;i++){
			Step[i]=i;
			StepChi2[i]=stepChi2.at(i);
		}
		HistChi2->Fill(Chi2);
		/*
		auto LdFit = PCor+Pi1Cor;
		auto TVLdFit = LdFit.Vect();
		auto TVLdMeas = TVPMeas + TVPi1Meas;
		PLdMeas = TVLdMeas.Mag();
		InvMLdCor = LdFit.Mag();	
		PLdFit = TVLdFit.Mag();
		ThLdFit = TVLdFit.Theta();
		PhLdFit = TVLdFit.Phi();
		auto TVPCor = PCor.Vect();
		PPCor = TVPCor.Mag();
		ThPCor = TVPCor.Theta();
		PhPCor = TVPCor.Phi();
		
		auto TVPi1Cor = Pi1Cor.Vect();
		PPi1Cor = TVPi1Cor.Mag();
		ThPi1Cor = TVPi1Cor.Theta();
		PhPi1Cor = TVPi1Cor.Phi();
		auto XiFit = LdCor + Pi2Meas;
		FourVectorFitter KFXi(LdCor,Pi2Meas,XiFit);
		double rL = ResPLd*PLdCor;
		double VarianceXi[8] = {rp*rp,ResThLd*ResThLd,ResPhLd*ResPhLd,rpi2*rpi2,ResTh*ResTh,ResPh*ResPh,0,0};
		KFXi.SetVariance(VarianceXi);
		KFXi.SetInvMass(mXi);	
		KFXi.SetMaximumStep(800);
		Chi2Xi=		KFXi.DoKinematicFit();
		
		auto PullXi = KFXi.GetPull();
		PullPLdCor = PullXi.at(0);
		PullThLdCor = PullXi.at(1);
		PullPhLdCor = PullXi.at(2);
		PullPPi2 = PullXi.at(3);
		PullThPi2 = PullXi.at(4);
		PullPhPi2 = PullXi.at(5);
		HistPullPLdCor->Fill(PullPLdCor);
		HistPullThLdCor->Fill(PullThLdCor);
		HistPullPhLdCor->Fill(PullPhLdCor);
		HistPullPPi2->Fill(PullPPi2);
		HistPullThPi2->Fill(PullThPi2);
		HistPullPhPi2->Fill(PullPhPi2);

		auto contXi = KFXi.GetFittedLV();
		auto LdCorCor = contXi.at(0);
		auto Pi2Cor = contXi.at(1);
		auto XiCor = LdCorCor + Pi2Cor;
		double probXi = ROOT::Math::chisquared_cdf(Chi2Xi,KFXi.GetNDF());
		HistChi2Xi->Fill(Chi2Xi);
		HistProbXi->Fill(probXi);
		HistChi2Xi2D->Fill(Chi2Xi,XiCor.Mag());
		HistProbXi2D->Fill(probXi,XiCor.Mag());
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
		InvMXiCor = (LdCor + Pi2).Mag();	
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
		
		
		
		HistDifLdP->Fill(PLdCor-PLdFit);
		HistDifLdTh->Fill(ThLdCor-ThLdFit);
		HistDifLdPh->Fill(PhLdCor - PhLdFit);
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
	*/



		double prob = ROOT::Math::chisquared_cdf(Chi2,KFLd.GetNDF());
		HistProb->Fill(prob);
		HistProb2D->Fill(prob,LdCor.Mag());
		HistChi22D->Fill(Chi2,LdCor.Mag());
		auto MM = Xi - XiRecon;
//		auto MMCor = Xi - XiCor;
		auto MMCor = MM;
		HistMM->Fill(MM.Mag());
		HistMMCor->Fill(MMCor.Mag());
		tree->Fill();	
		auto FIni = KFLd.GetInitialConstraints();
		auto FAfter = KFLd.GetKinematicConstraints();
		double Ipx=FIni.at(0);
		double Ipy=FIni.at(1);
		double Ipz=FIni.at(2);
		double IE=FIni.at(3);
		double Apx=FAfter.at(0);
		double Apy=FAfter.at(1);
		double Apz=FAfter.at(2);
		double AE=FAfter.at(3);
		HistPx->Fill(Ipx);
		HistPy->Fill(Ipy);
		HistPz->Fill(Ipz);
		HistE->Fill(IE);
		HistDPx->Fill(Apx);
		HistDPy->Fill(Apy);
		HistDPz->Fill(Apz);
		HistDE->Fill(AE);
	}
	HistLdP->SetLineColor(kRed);
	HistLdTh->SetLineColor(kRed);
	HistLdPh->SetLineColor(kRed);
	HistPP->SetLineColor(kRed);
	HistPTh->SetLineColor(kRed);
	HistPPh->SetLineColor(kRed);
	HistPi1P->SetLineColor(kRed);
	HistPi1Th->SetLineColor(kRed);
	HistPi1Ph->SetLineColor(kRed);
	HistPTh->Fit("gaus","Q0");
	HistLdP->Fit("gaus","Q0");
	/*
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	c1->cd(1);
	HistLd->Draw();
	HistLd->Fit("gaus","Q0");
	c1->cd(2);
	HistXi->Draw();
	HistXi->Fit("gaus","Q0");
	c1->cd(3);
	HistLdFit->Draw();
	HistLdFit->Fit("gaus","Q0");
	c1->cd(4);
	HistXiFit->Draw();
	HistXiFit->Fit("gaus","Q0");
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(3,3);
	c2->cd(1);
	HistLdP->Draw();
	HistLdP->SetLineColor(kGreen);
	c2->cd(2);
	HistLdTh->Draw();
	HistLdTh->Fit("gaus","Q0");
	c2->cd(3);
	HistLdPh->Draw();
	HistLdPh->Fit("gaus","Q0");
	c2->cd(4);
	HistPP->Draw();
	HistPP->Fit("gaus","Q0");
	c2->cd(5);
	HistPTh->Draw();
	c2->cd(6);
	HistPPh->Draw();
	HistPPh->Fit("gaus","Q0");
	c2->cd(7);
	HistPi1P->Draw();
	HistPi1P->Fit("gaus","Q0");
	c2->cd(8);
	HistPi1Th->Draw();
	HistPi1Th->Fit("gaus","Q0");
	c2->cd(9);
	HistPi1Ph->Draw();
	HistPi1Ph->Fit("gaus","Q0");
*/
	/*
	TCanvas* c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(3,3);
	c3->cd(1);
	HistDLdP->Draw();
	HistDLdP->Fit("gaus","Q0");
	c3->cd(2);
	HistDLdTh->Draw();
	HistDLdTh->Fit("gaus","Q0");
	c3->cd(3);
	HistDLdPh->Draw();
	HistDLdPh->Fit("gaus","Q0");
	c3->cd(4);
	HistDPP->Draw();
	HistDPP->Fit("gaus","Q0");
	c3->cd(5);
	HistDPTh->Draw();
	HistDPTh->Fit("gaus","Q0");
	c3->cd(6);
	HistDPPh->Draw();
	HistDPPh->Fit("gaus","Q0");
	c3->cd(7);
	HistDPi1P->Draw();
	HistDPi1P->Fit("gaus","Q0");
	c3->cd(8);
	HistDPi1Th->Draw();
	HistDPi1Th->Fit("gaus","Q0");
	c3->cd(9);
	HistDPi1Ph->Draw();
	HistDPi1Ph->Fit("gaus","Q0");
	*/
	TCanvas* c4 = new TCanvas("c4","c4",1200,600);
	c4->Divide(3,3);
	c4->cd(1);
	HistDifLdP->Draw();
	HistDifLdP->Fit("gaus","Q0");
	c4->cd(2);
	HistDifLdTh->Draw();
	HistDifLdTh->Fit("gaus","Q0");
	c4->cd(3);
	HistDifLdPh->Draw();
	HistDifLdPh->Fit("gaus","Q0");
	c4->cd(4);
	HistDifPP->Draw();
	HistDifPP->Fit("gaus","Q0");
	c4->cd(5);
	HistDifPTh->Draw();
	HistDifPTh->Fit("gaus","Q0");
	c4->cd(6);
	HistDifPPh->Draw();
	HistDifPPh->Fit("gaus","Q0");
	c4->cd(7);
	HistDifPi1P->Draw();
	HistDifPi1P->Fit("gaus","Q0");
	c4->cd(8);
	HistDifPi1Th->Draw();
	HistDifPi1Th->Fit("gaus","Q0");
	c4->cd(9);
	HistDifPi1Ph->Draw();
	HistDifPi1Ph->Fit("gaus","Q0");
	

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
	HistXiP->SetLineColor(kRed);
	HistDXiP->Draw("same");
	c8->cd(2);
	HistXiTh->Draw();
	HistXiTh->SetLineColor(kRed);
	HistDXiTh->Draw("same");
	c8->cd(3);
	HistXiPh->Draw();
	HistXiPh->SetLineColor(kRed);
	HistDXiPh->Draw("same");
	c8->cd(4);
	HistLdCorP->Draw();
	HistLdCorP->SetLineColor(kRed);
	HistDLdCorP->Draw("same");
	c8->cd(5);
	HistLdCorTh->Draw();
	HistLdCorTh->SetLineColor(kRed);
	HistDLdCorTh->Draw("same");
	c8->cd(6);
	HistLdCorPh->Draw();
	HistLdCorPh->SetLineColor(kRed);
	HistDLdCorPh->Draw("same");
	c8->cd(7);
	HistPi2P->Draw();
	HistPi2P->SetLineColor(kRed);
	HistDPi2P->Draw("same");
	c8->cd(8);
	HistPi2Th->Draw();
	HistPi2Th->SetLineColor(kRed);
	HistDPi2Th->Draw("same");
	c8->cd(9);
	HistPi2Ph->Draw();
	HistPi2Ph->SetLineColor(kRed);
	HistDPi2Ph->Draw("same");
	TCanvas* c9 = new TCanvas("c9","c9",600,600);
	HistMMCor->SetLineColor(kRed);
	HistMMCor->Draw();
	HistMM->Draw("same");
	TCanvas* c10 = new TCanvas("c10","c10",600,600);
	c10->Divide(2,1);
	c10->cd(1);
	HistChi2Xi->Draw();
	c10->cd(2);
	HistChi2Xi2D->Draw("colz");
	TCanvas* c11 = new TCanvas("c11","c11",600,600);
	c11->Divide(2,1);
	c11->cd(1);
	HistProbXi->Draw();
	c11->cd(2);
	HistProbXi2D->Draw("colz");
	TCanvas* c12 = new TCanvas("c12","c12",1200,600);
	c12->Divide(3,3);
	c12->cd(1);
	HistPullPLd->Draw();
	HistPullPLd->Fit("gaus");
	c12->cd(2);
	HistPullThLd->Draw();
	HistPullThLd->Fit("gaus");
	c12->cd(3);
	HistPullPhLd->Draw();
	HistPullPhLd->Fit("gaus");
	c12->cd(4);
	HistPullPP->Draw();
	HistPullPP->Fit("gaus");
	c12->cd(5);
	HistPullThP->Draw();
	HistPullThP->Fit("gaus");
	c12->cd(6);
	HistPullPhP->Draw();
	HistPullPhP->Fit("gaus");
	c12->cd(7);
	HistPullPPi1->Draw();
	HistPullPPi1->Fit("gaus");
	c12->cd(8);
	HistPullThPi1->Draw();
	HistPullThPi1->Fit("gaus");
	c12->cd(9);
	HistPullPhPi1->Draw();
	HistPullPhPi1->Fit("gaus");
	TCanvas* c13 = new TCanvas("c13","c13",1200,600);
	c13->Divide(3,2);
	c13->cd(1);
	HistPullPLdCor->Draw();
	HistPullPLdCor->Fit("gaus");
	c13->cd(2);
	HistPullThLdCor->Draw();
	HistPullThLdCor->Fit("gaus");
	c13->cd(3);
	HistPullPhLdCor->Draw();
	HistPullPhLdCor->Fit("gaus");
	c13->cd(4);
	HistPullPPi2->Draw();
	HistPullPPi2->Fit("gaus");
	c13->cd(5);
	HistPullThPi2->Draw();
	HistPullThPi2->Fit("gaus");
	c13->cd(6);
	HistPullPhPi2->Draw();
	HistPullPhPi2->Fit("gaus");
	TCanvas* c14 = new TCanvas("c14","c14",1200,600);
	c14->Divide(9,9);
		for(int col=0;col<9;++col){
			for(int row=0;row<col;++row){
				c14->cd(9*col + row + 1);
				HistCor[9*col + row]->Draw("colz");
			}
		}
	TCanvas* c15 = new TCanvas("c15","c15",1200,800);
	c15->Divide(2,2);
	c15->cd(1);
	HistDPx->Draw();
	HistDPx->SetLineColor(kRed);
	HistPx->Draw("SAME");
	c15->cd(2);
	HistDPy->Draw();
	HistDPy->SetLineColor(kRed);
	HistPy->Draw("SAME");
	c15->cd(3);
	HistDPz->Draw();
	HistDPz->SetLineColor(kRed);
	HistPz->Draw("SAME");
	c15->cd(4);
	HistDE->Draw();
	HistDE->SetLineColor(kRed);
	HistE->Draw("SAME");
}

