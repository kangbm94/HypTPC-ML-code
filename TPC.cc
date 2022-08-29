#include "src/TPCManager.cc"
#include "src/TPCCorrectionMapMaker.cc"

void TPC(){
	T.InitializeHistograms();
	cout<<"ViewTPC(TString Filename)"<<endl;
	cout<<"TagTPCEvents(TString Filename, TString OutFileName)"<<endl;
	cout<<"TagTPCTracks(TString Filename)"<<endl;
	cout<<"TestML(TString Filename)"<<endl;
	cout<<"ConvertRealTPC(int runnum)"<<endl;
	cout<<"CheckRealTPCTraining(TString Filename)"<<endl;
	gStyle->SetPalette(22);
}
