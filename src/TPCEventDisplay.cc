#include "TPCManager.cc"
int runnum = 5000;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
TString dir = base_dir;
//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
//TString tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
TString tpcfile = dir + Form("run05641_DstBeamRemover.root");
//TString tpcfile = dir + Form("run0%d_DstSelectedTPCHelixTracking.root ",runnum);
void TPCEventDisplay(){
}
void TPCEventDisplayRaw(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",800,1200);
	c1->Divide(2,1);
	bool Single = true;
//	Single = false;
	cout<<ent<<endl;
	for(int i=0;i<ent;++i){
		if(i%1==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		if(Single){
			cout<<"Drawing..."<<i<<endl;
			c1->cd(1);
			h->Draw("colz");
			h->SetTitle(Form("Event%d",evnum));
			c1->cd(2);
			h2->Draw("colz");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
		}
	}
//	h->Draw("colz");
//	h->SetTitle(Form("Run%d",runnum));
}
void TPCEventDisplayHelix(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile,"tree");
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
void DisplayBeamRemver(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	auto hf = gTPCManager.GetPadHistogramf();
	auto hf2 = gTPCManager.GetZYHistogramf();
	auto hb = gTPCManager.GetPadHistogramb();
	auto hb2 = gTPCManager.GetZYHistogramb();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
	c1->Divide(2,2);
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,1);
	bool Single = true;
//	Single = false;
	cout<<ent<<endl;
	for(int i=0;i<ent;++i){
		if(i%1==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=gTPCManager.GetNhits();
		int nhf=gTPCManager.GetNhitsf();
		int nhb=gTPCManager.GetNhitsb();
		if(!nhf or !nhb) continue;
//		if(nhf+ nhb ==nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		for(int itr=0;itr<nhf;++itr){
			gTPCManager.FillHistf(itr);
		}
		for(int itr=0;itr<nhb;++itr){
			gTPCManager.FillHistb(itr);
		}
		if(Single){
			cout<<"Drawing..."<<i<<endl;
			auto by = gTPCManager.GetBeamY();
			auto p0 = gTPCManager.GetBeamP0();
			auto p1 = gTPCManager.GetBeamP1();
			auto p2 = gTPCManager.GetBeamP2();
			c1->cd(1);
			hf->Draw("colz");
			hf->SetTitle(Form("Event%d",evnum));
			for(int i=0;i<by.size();++i){
				TF1* f = new TF1(Form("f%d",i),"pol2",-250,250);
				f->SetParameter(0,p0[i]);
				f->SetParameter(1,p1[i]);
				f->SetParameter(2,p2[i]);
				f->SetLineColor(i+2);
				f->SetLineWidth(3);
				f->Draw("same");
			}
			c1->cd(2);
			hf2->Draw("colz");
			for(int i=0;i<by.size();++i){
				TF1* f2 = new TF1(Form("f2%d",i),"pol1",-250,250);
				f2->SetParameter(0,by[i]);
				f2->SetParameter(1,0);
				f2->SetLineColor(i+2);
				f2->SetLineWidth(3);
				f2->Draw("same");
			}
			c1->cd(3);
			hb->Draw("colz");
			hb->SetTitle(Form("Event%d",evnum));
			c1->cd(4);
			hb2->Draw("colz");
			c1->Modified();
			c1->Update();
			c2->cd(1);
			h->Draw("colz");
			h->SetTitle(Form("Event%d",evnum));
			for(int i=0;i<by.size();++i){
				TF1* f = new TF1(Form("f%d",i),"pol2",-250,250);
				f->SetParameter(0,p0[i]);
				f->SetParameter(1,p1[i]);
				f->SetParameter(2,p2[i]);
				f->SetLineColor(i+2);
				f->SetLineWidth(3);
				f->Draw("same");
			}
			c2->cd(2);
			h2->Draw("colz");
			for(int i=0;i<by.size();++i){
				TF1* f2 = new TF1(Form("f2%d",i),"pol1",-250,250);
				f2->SetParameter(0,by[i]);
				f2->SetParameter(1,0);
				f2->SetLineColor(i+2);
				f2->SetLineWidth(3);
				f2->Draw("same");
			}
			
			c2->Modified();
			c2->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
		}
	}
//	h->Draw("colz");
//	h->SetTitle(Form("Run%d",runnum));
}
