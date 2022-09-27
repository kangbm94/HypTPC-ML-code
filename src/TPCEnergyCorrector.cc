#include "TPCManager.cc"
void TPCEnergyCorrector(){
	cout<<"PadEnergy(TString filename)"<<endl;
	cout<<"ClusterPos(TString filename)"<<endl;
}
void ClusterPos(TString filename){
	gTPCManager.LoadFile(filename);
	gTPCManager.LoadBcOut();
	gTPCManager.LoadClusterChain();
	gTPCManager.InitializeHistograms();
	int nev = gTPCManager.GetNEvent();
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	//	c1->DrawFrame(-250,250,-250,250);
	//	c1->SetGrid();
	TGraph* BcGr;
	TGraph* ClGr;
	TH1I* hist[6]; 
	for(int i=0;i<6;++i){
		TString title = Form("Run05756_L%d",i+27);
		hist[i]=new TH1I(title,title,11,-5,6);
	}
	for(int i=0;i<nev;++i){
		if(i%10000==0){
			cout<<i<<" th"<<endl;
		}
		gTPCManager.SetEvent(i);
		int nh = gTPCManager.GetNhits(0);
		int ncl = gTPCManager.GetNhits(1);
		int ncl2=ncl;
		if(gTPCManager.GetBCnt()!=1||ncl>35){
			continue;
		}
		Track Bc = gTPCManager.GetTrack();
		double bcx=Bc.GetPosition(-K18HS).X();
		int cnt =0;
		for(int j=0;j<ncl2;++j){
			bool ct = false;
			TVector3 ClPos = gTPCManager.GetClusterPosition(cnt);
			double clx=ClPos.X(),clz=ClPos.Z();
			double bcx=Bc.GetPosition(clz).X();
			int l,r;
			int bcpad = tpc::findPadID(clz,bcx);
			GetLayerRow(bcpad,l,r);
			if(l>25){
				for(int k=0;k<nh;++k){
					int padID=gTPCManager.GetPadID(k);
					int lay,row;
					GetLayerRow(padID,lay,row);
					int cs = gTPCManager.GetClSize(cnt);
					if(lay==l&&cs<2){
						hist[l-26]->Fill(row-r);
					}
				}
			}
			cnt++;
		}
	}
	TCanvas* c2 = new TCanvas("c2","c2",1500,800);
	c2->Divide(3,2);
	for(int i=0;i<6;++i){
		c2->cd(i+1);
		hist[i]->Draw();
	}
}

void PadEnergy(TString filename){
	gTPCManager.LoadFile(filename);
	gTPCManager.LoadBcOut();
	gTPCManager.LoadClusterChain();
	gTPCManager.InitializeHistograms();
	int nev = gTPCManager.GetNEvent();
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->DrawFrame(-250,250,-250,250);
	c1->SetGrid();
	TGraph* BcGr;
	TGraph* ClGr;
	for(int i=0;i<nev;++i){
		if(i%10000==0){
			cout<<i<<" th"<<endl;
		}
		gTPCManager.SetEvent(i);
		int nh = gTPCManager.GetNhits(0);
		int ncl = gTPCManager.GetNhits(1);
		int ncl2=ncl;
		if(gTPCManager.GetBCnt()!=1||ncl>35){
			continue;
		}
//		Track Bc = gTPCManager.GetTrack();
//		double bcx=Bc.GetPosition(-K18HS).X();
//		cout<<"nh = "<<nh<<endl;
		for(int j=0;j<nh;++j){
//			int padID=gTPCManager.GetPadID(j);
//			double de = gTPCManager.GetDE(j);
//				double de=1;
			gTPCManager.FillHist(j);
//			gTPCManager.SetPadContent(padID,de);
		}

		int cnt =0;
		/*
		TF1* Bcf = new TF1("Bcf","pol1",-300,300);	
		TF1* Clf = new TF1("Clf","pol1",-300,300);
		Clf->SetLineColor(kViolet);
		Bcf->SetLineWidth(3);
		Clf->SetLineWidth(3);
		auto* h = gTPCManager.GetPadHistogram();
		c1->cd();
		h->SetTitle(Form("Run05754,Event%d",i));
		h->Draw("col");
		//		gTPCManager.ClearHistogram();
		gPad->SetLogz();
		//		h->GetYaxis()->SetRange(49,51);
		BcGr->Draw("Psame");
		ClGr->Draw("Psame");
		BcGr->SetMarkerStyle(1);
		//		BcGr->SetMarkerColor(kMagenta);
		BcGr->SetMarkerSize(0.1);
		BcGr->Fit("Bcf");
		ClGr->SetMarkerStyle(37);
		ClGr->SetMarkerColor(kGreen);
		ClGr->SetMarkerSize(0.3);
		ClGr->Fit("Clf");
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
		gTPCManager.ClearHistogram();
		delete BcGr;
		delete ClGr;
		cout<<"GrDeleted"<<endl;*/
	}
		auto* h = gTPCManager.GetPadHistogram();
		h->Draw("colz");
}
