#include "TPCManager.cc"
#include "../include/TPCEventDisplay.hh"
void Display(){
//	int runnum = 5641;
	int runnum = 5721;
//	TString tpcdir = base_dir+"MayRun/rootfiles/CH2/TPC/";
	TString tpcdir = base_dir+"MayRun/rootfiles/Callibration/";
	TString tpcfilename = Form("run0%d_TPCHit.root",runnum);
	gTPCManager.LoadFile(tpcdir+tpcfilename);
	gStyle->SetPalette(kCividis);
	gStyle->SetNumberContours(255);
	gTPCManager.InitializeTPC();
	gTPCManager.DrawTPC();
	for(int i=0;i<10;++i){
		gTPCManager.SetEvent(i);
		gTPCManager.LoadTPC3D();
		gSystem->ProcessEvents();
		cin.ignore();
	}
	//	new TPCMainFrame(gClient->GetRoot(),1600,800);
}
