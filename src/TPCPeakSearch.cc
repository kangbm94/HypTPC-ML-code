#include "TPCManager.cc"
int runnum = 5000;
TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
TString tpcfile = dir + Form("SelectedHelix.root");
//TString tpcfile = dir + Form("run0%d_DstSelectedTPCHelixTracking.root ",runnum);
int SearchPeaks(TH1D* hist,vector<int> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.1;
	int npeaks = spec.Search(hist,sig,"",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}










void TPCEventDisplay(){
}
void TPCEventDisplayHelix(){
	TF1* fpol = new TF1("fpol2","pol2",-250,250);
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile,"tree");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	auto h3 = gTPCManager.GetYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",1000,700);
	c1->Divide(2,1);
	TCanvas* c2 = new TCanvas("c2","c2",800,500);
	c2->Divide(2,1);
	bool Single = true;
	//	Single = false;
	TFile* file = new TFile("SelectedEvents.root");
	TTree* tree = (TTree*)file->Get("tree");
	int xirn,xiev;
	double xim2;
	tree->SetBranchAddress("runnum",&xirn);
	tree->SetBranchAddress("evnum",&xiev);
	tree->SetBranchAddress("XiM2",&xim2);
	int xient = tree->GetEntries();
	TH2D* hist_beam = new TH2D("bh","bh",100,-250,250,100,-100,100);
	for(int i=0;i<ent;++i){
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int runnum = gTPCManager.GetRunnum();
		bool go = false;
		for(int j=0;j<xient;j++){
			tree->GetEntry(j);
			if(runnum == xirn and evnum == xiev and abs(xim2-1314)<0.05){cout<<xirn<<endl; go = false;break;}
			else{
				go = true;
			}
		}
		//		if(go) continue;
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		gTPCManager.ReconEvent();

		if(!gTPCManager.XiEvent()) continue;
		gTPCManager.InitializeHelix();
		vector<int>peaks;

		if(Single){
			for(int itr=0;itr<nh;++itr){
				gTPCManager.FillHist(itr);
			}
			int np = SearchPeaks(h3,peaks); 
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
			double width =  5;
			TMinuit *minuit;
			double arglist[5];
			arglist[0]=5.89;
			TString name[5] = {"cx", "cy", "z0", "r", "dz"};
			double err[5]={-999.,-999.,-999.,-999.,-999.};
			int Err;
			double bnd1,bnd2;
			double cx,cy,z0,r,dz;
			double par[5] = {6000,-6000,0,0,0};
			Int_t ierflg = 0;
			c1->cd(1);
			h->Draw();
			h->SetTitle("CircleFit");
			//	h->SetTitle(Form("MissMass = %f",xim2));
			//	h->SetTitle(Form("MissMass = %f",xim2));
//			gTPCManager.DrawHelix();
			gTPCManager.DrawVertex();
			c1->cd(2);
			h2->Draw("col");
//			gTPCManager.DrawHelixZY();
			gTPCManager.DrawVertexZY();
			c2->cd(1);
			h3->Draw();
			fpol->SetParLimits(2,1/5000,1/7000);
			for(int ipk=0;ipk<np;++ipk){
				gHitPos.clear();
				gRes.clear();
				if(abs(peaks[ipk])<60) continue;
				gTPCManager.FillgHitPos(peaks[ipk],width);	
				if(gHitPos.size()<10) continue;
				double p0,p1,p2;
				for(auto g:gHitPos){
					hist_beam->Fill(g.Z(),g.X());
				}
				hist_beam->Fit("fpol2","QR0");
				hist_beam->Reset("ICES");
				p0=fpol->GetParameter(0);
				p1=fpol->GetParameter(1);
				p2=fpol->GetParameter(2);
				for(auto g:gHitPos){
					double val = p0+p1*g.Z()+p2*g.Z()*g.Z();
					if(abs(val-g.X())<20)hist_beam->Fill(g.Z(),g.X());
				}
				hist_beam->Fit("fpol2","QR0");
				hist_beam->Reset("ICES");
				p0=fpol->GetParameter(0);
				p1=fpol->GetParameter(1);
				p2=fpol->GetParameter(2);
				for(auto g:gHitPos){
					double val = p0+p1*g.Z()+p2*g.Z()*g.Z();
					if(abs(val-g.X())<8)hist_beam->Fill(g.Z(),g.X());
				}

				break;
			}
			c2->cd(2);
			hist_beam->Draw("colz");
			hist_beam->Fit("fpol2","R");
			hist_beam->Reset("ICES");
			for(auto g:gHitPos){
					hist_beam->Fill(g.Z(),g.X());
			}
			c1->Modified();
			c1->Update();
			c2->Modified();
			c2->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"Searching..."<<endl;
			hist_beam->Reset("ICES");
			gTPCManager.ClearHistogram();
		}
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
}
void DisplayY(){
	gTPCManager.LoadClusterFile(tpcfile);
	auto h2 = gTPCManager.GetZYHistogram();
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);



}
