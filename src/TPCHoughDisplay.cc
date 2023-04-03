#include "TPCManager.cc"
int runnum = 5642;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
int SearchPeaks(TH1D* hist,vector<double> &peaks){
	TSpectrum spec(30);
	double sig=0.8,th=0.1;
	//	int npeaks = spec.Search(hist,sig,"goff",th);
	int npeaks = spec.Search(hist,sig,"",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}
TString dir = base_dir;
const Int_t nBin_rdiff = 220;
const Double_t rdiff_min = -220.;
const Double_t rdiff_max = 220.;
const Int_t nBin_theta = 180;
const Double_t theta_min = -0.3*acos(-1);
const Double_t theta_max = 0.3*acos(-1);//Charge < -1
const Int_t nBin_p = 100;
const Double_t p_min = 1600.;//MeV/c
const Double_t p_max = 2000.;//MeV/c
const int    thetaY_ndiv =  180;
const double thetaY_min  =  60.;
const double thetaY_max  =  120.;
const int    r_ndiv =  1000;
const double r_min  = -500.;
const double r_max  =  500.;


TString tpcfile;
void TPCHoughDisplay(){
	cout<<"TPCEventDisplayZY(int ievt)"<<endl;
	
	int runnum = 5641;
	tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
	tpcfile = "test.root"; 
	gTPCManager.LoadClusterFile(tpcfile,"tpc");

	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.InitializeTPC();
	gTPCManager.InitializeHoughHistograms();
}







void TPCHoughDisplayZY(int ievt){
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TH1D* Yhist= new TH1D("histY","histY",80,-500,500);
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
	TCanvas* c3 = new TCanvas("c3","c3",50,50,900,900);
	TCanvas* c4 = new TCanvas("c4","c4",650,650,300,300);
//	TCanvas* c5 = new TCanvas("c5","c5",350,350,600,600);
	gTPCManager.SetEvent(ievt);
	gTPCManager.ClearHistogram();
	int evnum = gTPCManager.GetEvnum();
	bool go = false;
	int nh=0;
	//		if(!gTPCManager.TagTrig(23)) continue;
	nh=gTPCManager.GetNhits(1);
	bool acc = false;
	//		gTPCManager.InitializeAccidental();
	//		accidental_dist = gTPCManager.GetAccidentalDist();	
	for(int itr=0;itr<nh;++itr){
		//				if(gTPCManager.GetClDe(itr)>60)
		gTPCManager.FillHist(itr);
		auto pos = gTPCManager.GetPosition(itr);
		double x = pos.X(),y=pos.Y(),z=pos.Z();
		if(-0 < x and x < 0 and -0 < y and y < 0)continue;
			Yhist->Fill(y);
//		if( gTPCManager.GetHoughFlag()->at(itr)<990)
			//				cout<<"ADist : "<<accidental_dist->at(itr)<<endl;
	}
	gTPCManager.DoZYHough();
	cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
	c1->Clear("D");
	c1->cd(1);
	h->Draw("colz");
	h->SetTitle("CircleFit");
	gTPCManager.GetPadHistogram()->Draw("colz");
	c1->cd(2);
	h2->Draw("colz");
	gTPCManager.GetZYHistogram()->Draw("colz");
	gTPCManager.DrawZYHough();
	c1->Modified();
	c1->Update();
	c3->cd();
	c3->Clear("D");
	gTPCManager.DrawTPC();
	gTPCManager.LoadTPC3D();
	c3->Modified();
	c3->Update();
	c4->cd();
	vector<double>peaks;
	Yhist->Draw();
	int npeaks = SearchPeaks(Yhist,peaks);
	double YWindow = 25;
	double max = Yhist->GetMaximum();
	vector<vector<TVector3>> peakarr;
	peakarr.resize(npeaks);
	for(int ip = 0; ip < npeaks; ++ip){
		double peak = peaks.at(ip);
		for(int itr=0;itr<nh;++itr){
			auto pos = gTPCManager.GetPosition(itr);
			double x = pos.X(),y=pos.Y(),z=pos.Z();
			if(abs(peak - y)<YWindow){
				peakarr.at(ip).push_back(pos);	
			}
		}
		TLine* left = new TLine(peak-YWindow,0,peak-YWindow,max);
		TLine* right = new TLine(peak+YWindow,0,peak+YWindow,max);
		left->Draw("same");
		left->SetLineColor(2);
		right->Draw("same");
		right->SetLineColor(3);
	}
	c4->Modified();
	c4->Update();









	gSystem->ProcessEvents();
}


