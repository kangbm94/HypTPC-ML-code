#include "src/TPCManager.cc"
#include "src/TPCCorrectionMapMaker.cc"

void TPC(){
	gTPCManager.InitializeHistograms();
	cout<<"ViewTPC(TString Filename)"<<endl;
	cout<<"TagTPCEvents(TString Filename, TString OutFileName)"<<endl;
	cout<<"TagTPCTracks(TString Filename)"<<endl;
	cout<<"TestML(TString Filename)"<<endl;
	cout<<"ConvertRealTPC(int runnum)"<<endl;
	cout<<"CheckRealTPCTraining(TString Filename)"<<endl;
	TPCCorrectionMapMaker();
	gStyle->SetPalette(22);
}
