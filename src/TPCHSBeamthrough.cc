#include "TPCManager.cc"
void TPCHSBeamthrough(){
	TString tpcfile;
	TString dir = base_dir+ "Cor1st/";
	TString name = "HSpi500.root"; 
	tpcfile = dir + name;
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	TFile* file = new TFile(dir+"Fit"+name,"recreate");
	double mom0,mom0Fit;
	TTree* tree = new TTree("tree","tree");
	tree->Branch("mom0",&mom0);
	tree->Branch("mom0Fit",&mom0Fit);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();
	gTPCManager.ClearHistogram();
	for(int i=0;i<ent;++i){
		gTPCManager.SetEvent(i);
//		if(!gTPCManager.Acpt()) continue;
		gTPCManager.MakeHSTrack();
		auto HS = gTPCManager.GetHSTracks();
		if(HS.size()==1){
			auto HST = HS.at(0);
			mom0 = HST.GetMom0();
			mom0Fit = HST.GetMom0Fit();
			tree->Fill();
		}
		if(i%1000==0){
			cout<<"Event : "<<i<<" HSsize = "<<HS.size()<<endl;
		}
	}
	file->Write();

}
