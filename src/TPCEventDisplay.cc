#include "TPCManager.cc"
int runnum = 5000;
TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
TString tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
//TString tpcfile = dir + Form("run0%d_DstSelectedTPCHelixTracking.root ",runnum);
void TPCEventDisplay(){
}
void TPCEventDisplayRaw(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",800,1200);
	c1->Divide(2,1);
	bool Single = true;
//	Single = false;
	for(int i=0;i<ent;++i){
		if(i%10000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		if(Single){
			cout<<"Drawing"<<endl;
			c1->cd(1);
			h->Draw();
			h->SetTitle(Form("Event%d",evnum));
			c1->cd(2);
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
		}
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
}
void TPCEventDisplayHelix(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
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
		gTPCManager.InitializeHelix();
		gTPCManager.ReconEvent();
		if(!gTPCManager.XiEvent()) continue;
		if(Single){
			for(int itr=0;itr<nh;++itr){
				gTPCManager.FillHist(itr);
			}
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
			c1->cd(1);
			h->Draw();
			h->SetTitle("CircleFit");
		//	h->SetTitle(Form("MissMass = %f",xim2));
		//	h->SetTitle(Form("MissMass = %f",xim2));
			gTPCManager.DrawHelix();
			gTPCManager.DrawVertex();
			c1->cd(2);
			h2->Draw("col");
			gTPCManager.DrawHelixZY();
			gTPCManager.DrawVertexZY();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"Searching..."<<endl;
			gTPCManager.ClearHistogram();
		}
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
}
