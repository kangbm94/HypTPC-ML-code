#include "TPCManager.cc"
int runnum = 5641;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString dir = base_dir+"before_mod/";
TString dir = base_dir;//+"before_mod/";

TString tpcfile;
void TPCHitAngle(){
	tpcfile = dir + Form("run0%d_DstTPCHSKuramaHelixTracking.root",runnum);
//	tpcfile = "test.root"; 
	gTPCManager.LoadClusterFile(tpcfile,"tpc");

//	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
//	gTPCManager.InitializeTPC();
}
void MemoryTest(){
	int ent = gTPCManager.GetEntries();
//	ent = 100000;
	for(int i=0;i<ent;++i){
		if(i%10000==0){
			cout<<i<<endl;
		}
		gTPCManager.SetEvent(i);	
		gTPCManager.InitializeHelix();
	}
}
void DrawAngle(){
	for(int i=0;i<max_pad;++i){
		cout<<i<<endl;
		auto pos = gTPCManager.GetPadPosition(i);
		cout<<"pos: "<<pos.x()<<","<<pos.y()<<","<<pos.z()<<endl;
		double t = gTPCManager.GetHitAngle(pos);
		cout<<"t: "<<t<<endl;
		gTPCManager.SetPadContent(i,t);
	}
	gTPCManager.DrawHist();
}
void ResVsTheta(){
	TFile* file = new TFile(Form("run0%d_ResVsTheta.root",runnum),"recreate");
	TTree* tree = new TTree("tree","tree");
	vector<vector<double>> ang,padang,Res,ResX,ResY,ResZ,HitX,HitY,HitZ,CalX,CalY,CalZ;
	tree->Branch("TrackAngle",&ang);
	tree->Branch("PadAngle",&padang);
	tree->Branch("Residual",&Res);
	tree->Branch("ResidualX",&ResX);
	tree->Branch("ResidualY",&ResY);
	tree->Branch("ResidualZ",&ResZ);
	tree->Branch("HitX",&HitX);
	tree->Branch("HitY",&HitY);
	tree->Branch("HitZ",&HitZ);
	tree->Branch("CalX",&CalX);
	tree->Branch("CalY",&CalY);
	tree->Branch("CalZ",&CalZ);
	TH2D* TrackAngleR = new TH2D("angleVRes","angleVRes",1000,-1.5,1.5,1000,-10,10);
	TH2D* TrackAngleRX = new TH2D("angleVResX","angleVResX",1000,-1.5,1.5,1000,-10,10);
	TH2D* TrackAngleRY = new TH2D("angleVResY","angleVResY",1000,-1.5,1.5,1000,-10,10);
	TH2D* TrackAngleRZ = new TH2D("angleVResZ","angleVResZ",1000,-1.5,1.5,1000,-10,10);
	int ent = gTPCManager.GetEntries();
//	ent = 100000;
	for(int i=0;i<ent;++i){
		if(i%10000==0){
			cout<<i<<endl;
		}
		gTPCManager.SetEvent(i);	
		gTPCManager.InitializeHelix();
		int nt = gTPCManager.GetNHelixTrack();
		ang.clear();
		padang.clear();
		Res.clear();
		ResX.clear();
		ResY.clear();
		ResZ.clear();
		HitX.clear();
		HitY.clear();
		HitZ.clear();
		CalX.clear();
		CalY.clear();
		CalZ.clear();
		for(int it=0;it<nt;++it){
			ang.push_back(gTPCManager.TrackHitAngle2(it));
			padang.push_back(gTPCManager.PadAngle(it));
			Res.push_back(gTPCManager.Residual(it));
			ResX.push_back(gTPCManager.ResidualX(it));
			ResY.push_back(gTPCManager.ResidualY(it));
			ResZ.push_back(gTPCManager.ResidualZ(it));
			HitX.push_back(gTPCManager.HitX(it));
			HitY.push_back(gTPCManager.HitY(it));
			HitZ.push_back(gTPCManager.HitZ(it));
			CalX.push_back(gTPCManager.CalX(it));
			CalY.push_back(gTPCManager.CalY(it));
			CalZ.push_back(gTPCManager.CalZ(it));
			int nh = ang[it].size();
			for(int ih=0;ih<nh;++ih){
				double angd = padang[it].at(ih);
				double Resd = Res[it].at(ih);
				double ResXd = ResX[it].at(ih);
				double ResYd = ResY[it].at(ih);
				double ResZd = ResZ[it].at(ih);
				TrackAngleR->Fill(angd,Resd);
				TrackAngleRX->Fill(angd,ResXd);
				TrackAngleRY->Fill(angd,ResYd);
				TrackAngleRZ->Fill(angd,ResZd);
			}
		}
		tree->Fill();
	}
	file->Write();
	TCanvas* c1 = new TCanvas("c1","c1",900,900);
	c1->Divide(2,2);
	c1->cd(1);
	TrackAngleR->Draw("colz");
	c1->cd(2);
	TrackAngleRX->Draw("colz");
	c1->cd(3);
	TrackAngleRY->Draw("colz");
	c1->cd(4);
	TrackAngleRZ->Draw("colz");
}
void DrawTrack(int evt_id,int tr_id=0){
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gTPCManager.ClearHistogram();
	gTPCManager.SetEvent(evt_id);	
	gTPCManager.InitializeHelix();
	int nt = gTPCManager.GetNHelixTrack();
	cout<<"NTracks : "<<nt<<endl;
	gTPCManager.DrawHist();
	TH1D* TrackAngle = new TH1D("angle","angle",100,-1,1);
	if(tr_id<nt){
		c1->cd(1);
		cout<<"Drawing..."<<endl;
		gTPCManager.DrawHelixHit(tr_id);
		gTPCManager.DrawHelixDir(tr_id);
		auto t = gTPCManager.TrackHitAngle(tr_id);
		auto t2 = gTPCManager.TrackHitAngle2(tr_id);
		int nhs = t.size();
		for(int ih=0;ih<nhs;++ih){
			double tt = t.at(ih)-t2.at(ih);
			TrackAngle->Fill(tt);
		}
	}
	TCanvas* c2 = new TCanvas("c2","c2",100,100,600,600);
	TrackAngle->Draw();
}


