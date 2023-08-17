#include "TPCManager.cc"
int runnum = 5856;
//int runnum = 58**;//0.4G pi/P:44/47, 0.5 Pi/P:64/66, 0.6 Pi/P:58,60,0.8 Pi/P:55/56
//int runnum = 58**;//-0.3G Pi 21, -0.4G 24, -0.5G: 28, -0.8G: 18
//int runnum = 5764;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString dir = base_dir+"/AllPosCor/";
TString dir = base_dir;
//TString dir = base_dir;//+"before_mod/";
static int gevnum = 0;
TString tpcfile;
void TPCBeamthroughSelection(){
	gStyle->SetOptStat(0);
	tpcfile = dir +"NoCorL29/"+ Form("run0%d_DstTPCHelixTracking.root",runnum);
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	cout<<"SelectBeam()"<<endl;
}

void SelectBeam(){
	int ent = gTPCManager.GetEntries();
	TFile* file = new TFile(Form("run0%d_Beamthrough.root",runnum),"recreate");
	TTree* tree = new TTree("tpc","tpc");
	double mom0;
	vector<double>layer;	
	vector<double>row;	
	vector<double>hitpos_x;	
	vector<double>hitpos_y;	
	vector<double>hitpos_z;	
	vector<double>residual_x;	
	vector<double>residual_y;	
	vector<double>residual_z;	
	tree->Branch("mom0",&mom0);
	tree->Branch("layer",&layer);
	tree->Branch("row",&row);
	tree->Branch("hitpos_x",&hitpos_x);
	tree->Branch("hitpos_y",&hitpos_y);
	tree->Branch("hitpos_z",&hitpos_z);
	tree->Branch("residual_x",&residual_x);
	tree->Branch("residual_y",&residual_y);
	tree->Branch("residual_z",&residual_z);
	for(int iev=0; iev<ent;++iev){
		gTPCManager.SetEvent(iev);
		gTPCManager.InitializeHelix();
		if(!gTPCManager.SelectBeamthrough()) continue;	
		layer.clear();
		row.clear();
		hitpos_x.clear();
		hitpos_y.clear();
		hitpos_z.clear();
		residual_x.clear();
		residual_y.clear();
		residual_z.clear();
		mom0 = gTPCManager.Getmom0(0);
		auto hl = gTPCManager.HitLayer(0);
		auto hr = gTPCManager.HitRow(0);
		auto hx = gTPCManager.HitX(0);
		auto hy = gTPCManager.HitY(0);
		auto hz = gTPCManager.HitZ(0);
		auto rx = gTPCManager.ResidualX(0);
		auto ry = gTPCManager.ResidualY(0);
		auto rz = gTPCManager.ResidualZ(0);
		int nh = hl.size();
		for(int ih = 0; ih < nh; ++ih){
			layer.push_back(hl[ih]);	
			row.push_back(hr[ih]);	
			hitpos_x.push_back(hx[ih]);	
			hitpos_y.push_back(hy[ih]);	
			hitpos_z.push_back(hz[ih]);	
			residual_x.push_back(rx[ih]);	
			residual_y.push_back(ry[ih]);	
			residual_z.push_back(rz[ih]);	
		}
		tree->Fill();
	}
	file->Write();
}
