#include "TPCManager.cc"
int runnum = 5641;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
int SearchPeaks(TH1D* hist,vector<double> &peaks){
	TSpectrum spec(30);
	double sig=2,th=0.1;
	//	int npeaks = spec.Search(hist,sig,"goff",th);
	int npeaks = spec.Search(hist,sig,"",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}
//TString dir = base_dir+"before_mod/";
TString dir = base_dir;//+"before_mod/";
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
vector<int> Xievnum;


TString tpcfile;
void TPCEventDisplay(){
	cout<<"TPCEventDisplayAccidental(int ievt)"<<endl;
	cout<<dir<<endl;	
	tpcfile = dir + Form("run0%d_DstTPCKuramaSelectedHelixTracking.root",runnum);
//	tpcfile = "test.root"; 
	gTPCManager.LoadClusterFile(tpcfile,"tpc");

	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.InitializeTPC();
	TFile* fileXi = new TFile("Sorted.root");
	TTree* treeXi = (TTree*)fileXi->Get("tree");
	int XiRunnum,XiEvnum;
	treeXi->SetBranchAddress("runnum",&XiRunnum);
	treeXi->SetBranchAddress("evnum",&XiEvnum);
	for(int i = 0; i< treeXi->GetEntries();++i){
		treeXi->GetEntry(i);
		if(XiRunnum == runnum) Xievnum.push_back(XiEvnum);
	}
	cout<<"TPCXiDisplay(int iXi)"<<endl;
	cout<<"NXi = "<<Xievnum.size()<<endl;
}

void TPCEventDisplayHelix(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile,"tree");
	//	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	gTPCManager.SetBetheProton();	
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
	TFile* file = new TFile("SelectedEvents.root");
	TTree* tree = (TTree*)file->Get("tree");
	int xirn,xiev;
	double xim2;
	tree->SetBranchAddress("runnum",&xirn);
	tree->SetBranchAddress("evnum",&xiev);
	tree->SetBranchAddress("XiM2",&xim2);
	int xient = tree->GetEntries();
	TCanvas* c2 = new TCanvas("c2","c1",1200,1200);
	c2->Divide(4,3);
	for(int i=0;i<ent;++i){
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int runnum = gTPCManager.GetRunnum();
		bool go = false;
		for(int j=0;j<xient;j++){
/*
			tree->GetEntry(j);
			if(runnum == xirn and evnum == xiev and abs(xim2-1314)<0.05){cout<<xirn<<endl; go = false;break;}
			else{
				go = true;
			}
			*/
		}
		//		if(go) continue;
		
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		gTPCManager.InitializeHelix();
		gTPCManager.ReconEvent();
		//		if(!gTPCManager.XiEvent()) continue;
		for(int itr=0;itr<nh;++itr){
			//				if(gTPCManager.GetClDe(itr)>60)
			gTPCManager.FillHist(itr);
		}
		if(h->GetEffectiveEntries()==0){
			gTPCManager.ClearHistogram();
			continue;
		}
		cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
		c1->cd(1);
		h->Draw("col");
		h->SetTitle("CircleFit");
		gTPCManager.DrawHelix();
		if(gTPCManager.LambdaEvent())gTPCManager.DrawVertex();
		c1->cd(2);
		h2->Draw("col");
		gTPCManager.DrawHelixZY();
		if(gTPCManager.LambdaEvent())gTPCManager.DrawVertexZY();
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
		cout<<"Searching..."<<endl;
		gTPCManager.ClearHistogram();
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
}






void TPCEventDisplayAccidental(int ievt){
//	tpcfile = dir + Form("run0%d_DstTPCHelixTrackingWTarget.root",runnum);
//	tpcfile = dir + Form("test.root");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TH1D* Yhist= new TH1D("histY","histY",500,-500,500);
	gTPCManager.SetBetheProton();	
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
	TCanvas* c3 = new TCanvas("c3","c3",50,50,900,900);
	TCanvas* c4 = new TCanvas("c4","c4",650,650,300,300);
//	TCanvas* c5 = new TCanvas("c5","c5",350,350,600,600);
	gTPCManager.SetEvent(ievt);
	/*
	while(gTPCManager.TagTrig(20)){
		++ievt;
		gTPCManager.SetEvent(ievt);
	}
	*/
	cout<<"Event = "<<ievt<<endl;
	gTPCManager.ClearHistogram();
	int evnum = gTPCManager.GetEvnum();
	gTPCManager.SetXiStarFlag(true);
	bool go = false;
	int nh=0;
	//		if(!gTPCManager.TagTrig(23)) continue;
	nh=gTPCManager.GetNhits(1);
	bool acc = false;
	gTPCManager.InitializeHelix();
	gTPCManager.ReconEvent();
	//		gTPCManager.InitializeAccidental();
	//		accidental_dist = gTPCManager.GetAccidentalDist();	
	for(int itr=0;itr<nh;++itr){
		gTPCManager.FillHist(itr);
		auto pos = gTPCManager.GetPosition(itr);
		double x = pos.X(),y=pos.Y(),z=pos.Z();
		if(-0 < x and x < 0 and -0 < y and y < 0)continue;
			Yhist->Fill(y);
	}
	c1->cd(1);
	h->Draw("colz");
	h->SetTitle("CircleFit");
	gTPCManager.DrawHelix();
	c1->cd(2);
	h2->Draw("colz");
	gTPCManager.DrawHelixZY();
	c1->Modified();
	c1->Update();
	c3->cd();
	c3->Clear("D");
	gTPCManager.DrawTPC();
	gTPCManager.LoadTPC3D();
	gTPCManager.LoadHelix3D();
	gTPCManager.LoadAccidental3D();
	gTPCManager.DrawHelix3D();
	gTPCManager.DrawAccidental3D();
	gTPCManager.DrawVertex3D();
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
	cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;

	gSystem->ProcessEvents();
}
void TPCXiDisplay(int i){
	int ievt = Xievnum.at(i);	
	TPCEventDisplayAccidental(ievt);
}


void TPCAccidentalCD(){
	int runnum = 5641;
	tpcfile = dir + Form("run0%d_DstTPCHelixTrackingWTarget.root",runnum);
	tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	TFile* file = new TFile("dists.root","recreate");
	TTree* tree = new TTree("tree","tree");
	vector<double>* hough_dist = new vector<double>;
	vector<double>* accidental_dist = new vector<double>;
	vector<int>* hough_flag = new vector<int>;
	tree->Branch("hough_dist",&hough_dist);
	tree->Branch("hough_flag",&hough_flag);
	tree->Branch("accidental_dist",&accidental_dist);
	int ent = gTPCManager.GetEntries();	
	for(int i=0;i<ent;++i){
		gTPCManager.SetEvent(i);
		if(i%1000==0) cout<<i<<" th event"<<endl;
		hough_dist = gTPCManager.GetHoughDist();	
		hough_flag = gTPCManager.GetHoughFlag();	
//		cout<<hough_dist->size()<<endl;
		accidental_dist = gTPCManager.GetAccidentalDist();	
//		cout<<accidental_dist->size()<<endl;
//		cout<<endl;
		tree->Fill();
	}
	file->Write();
}

