#include "TPCManager.cc"
int runnum = 6384;
TString dir = base_dir+"MayRun/rootfiles/Cosmic/";
TString tpcfile = dir + Form("run0%d_TPCHit.root",runnum);
void TPCEventDisplay(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	bool Single = true;
	Single = false;
	for(int i=0;i<ent;++i){
		if(i%10000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int nh=0;
		nh=gTPCManager.GetNhits(0);
		if(!nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		if(Single){
			h->Draw();
			h->SetTitle(Form("Event%d",i));
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
