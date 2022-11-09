#include "../include/RiemanAnalysis.hh"
#include "TPCManager.cc"
void Geometry(){
	int runnum = 5641;
	RiemanManager RS;
	TString tpcdir = base_dir+"MayRun/rootfiles/CH2/TPC/";
	TString tpcfilename = Form("run0%d_TPCHit.root",runnum);
	RS.LoadFile(tpcdir+tpcfilename);
	for(int i=0;i<10;++i){
		RS.SetEvent(i);
		RS.ProjectPoints();
		//		gTPCManager.LoadTPC3D();
		gSystem->ProcessEvents();
		cin.ignore();
	}
}
