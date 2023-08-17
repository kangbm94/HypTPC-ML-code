#include "TPCManager.cc"
//int runnum = 5641;
//int runnum = 58**;//0.4G pi/P:44/47, 0.5 Pi/P:64/66, 0.6 Pi/P:58,60,0.8 Pi/P:55/56
//int runnum = 58**;//-0.3G Pi 21, -0.4G 24, -0.5G: 28, -0.8G: 18
int runnum = 5721;
//int runnum = 5641;
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
//TString dir = base_dir+"/AllPosCor/";
//TString dir = base_dir+"/NoCorL29/";
TString dir = base_dir;//+"before_mod/";
vector<int> Xievnum;
vector<double> XIMM;
static int gevnum = 0;
TString tpcfile;
void TPCEventDisplay(){
//	cout<<"TPCEventDisplayAccidental(int ievt)"<<endl;
//	tpcfile = dir + Form("./dstfiles/run0%d_DstTPCHSKuramaSelectedHelixTracking.root",runnum);
		tpcfile = dir + Form("NoCor/HSp400.root");
//	tpcfile = dir + Form("NoCor/run0%d_DstTPCHelixTrackingHToF.root",runnum);
//	tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
	//tpcfile = dir +"NoCorL29/"+ Form("run0%d_DstTPCHelixTracking.root",runnum);
//	tpcfile = dir +"BeforeUpdate/"+ Form("run0%d_DstTPCHelixTracking.root",runnum);
//	tpcfile = dir +"AfterUpdate/"+ Form("run0%d_DstTPCHelixTracking.root",runnum);
//	tpcfile = dir +"AllCor/"+ Form("run0%d_DstTPCHelixTracking.root",runnum);
//	tpcfile = dir + Form("run0%d_DstTPCTracking.root",runnum);
//	tpcfile = dir + Form("run0%d_DstTPCHSKuramaHelixTracking.root",runnum);
//	tpcfile = "test.root"; 

	gStyle->SetPalette(1);
//	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.InitializeTPC();
	TFile* fileXi = new TFile("Sorted.root");
	TTree* treeXi = (TTree*)fileXi->Get("tree");
	int XiRunnum,XiEvnum;
	double XiMM;
	treeXi->SetBranchAddress("runnum",&XiRunnum);
	treeXi->SetBranchAddress("evnum",&XiEvnum);
	treeXi->SetBranchAddress("MM",&XiMM);
	for(int i = 0; i< treeXi->GetEntries();++i){
		treeXi->GetEntry(i);
		if(XiRunnum == runnum){
			Xievnum.push_back(XiEvnum);
			XIMM.push_back(XiMM);
		}
	}
	cout<<"TPCEventDisplayHelix()"<<endl;
	cout<<"TPCEventDisplayHelix(iev)"<<endl;
	cout<<"TPCEventDisplayLinear()"<<endl;
	cout<<"TPCEventDisplayLinear(iev)"<<endl;
	cout<<"TPCXiDisplay(int iXi)"<<endl;
	cout<<"NXi = "<<Xievnum.size()<<endl;
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
}
void TPCEventDisplayHelix(int iev){
	gStyle->SetOptStat(0);
//	gTPCManager.LoadClusterFile(tpcfile,"tree");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1;
	TCanvas* c2;
	if(!c1){
		c1= new TCanvas("c1","c1",800,600);
		c1->Divide(2,1);
	}
	int xirn,xiev;
	double xim2;
	gTPCManager.ClearHistogram();

	gTPCManager.SetEvent(iev);
	//while(!(gTPCManager.Acpt() and gTPCManager.NearTarget(0,20) and gTPCManager.SelectBeamthrough())){
//	while(!gTPCManager.SelectBeamthrough() and 0){
	while(!gTPCManager.Acpt() and 1){
		iev++;
		gTPCManager.SetEvent(iev);
	}
	gevnum = iev;
	int evnum = gTPCManager.GetEvnum();
	int runnum = gTPCManager.GetRunnum();
	bool go = false;

	int nh=0;
	nh=gTPCManager.GetNhits(0);
//	cout<<"NHits = "<<nh<<endl;
	for(int itr=0;itr<nh;++itr){
		gTPCManager.FillHistHits(itr);
	}
	nh=gTPCManager.GetNhits(1);
	for(int itr=0;itr<nh;++itr){
		gTPCManager.FillHist(itr);
	}
	gTPCManager.InitializeHelix();
	cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
	c1->cd(1);
	h->Draw("col");
	h->SetTitle("CircleFit");
	gTPCManager.DrawHelix();
//	gTPCManager.MakeHSTrack();

//	auto HS = gTPCManager.GetHSTracks();
/*
	if(HS.size()==1){
		auto HSt = HS.at(0);
		cout<<Form("pHS,mom0,mom0fit = (%g,%g,%g)",HSt.GetMom0(),gTPCManager.Getmom0(0),HSt.GetMom0Fit())<<endl;
	}
	gTPCManager.DrawHS();
	*/
	c1->cd(2);
	h2->Draw("col");
	gTPCManager.DrawHelixZY();
	gTPCManager.DrawHSZY();
	auto hists = gTPCManager.GetResHists();
	if(!c2){
		c2 = new TCanvas("c2","c2",300,200,1200,600);
		c2->Divide(2,1);
	}
	if(hists.size()>0){
		c2->cd(1);
		hists.at(0)->Draw("colz");
	}
	if(hists.size()==2){
		c2->cd(2);
		hists.at(1)->Draw("colz");
	}
	
/*
	c2->cd();
	gTPCManager.DrawTPC();
	gTPCManager.LoadHelix3D();
	gTPCManager.LoadAccidental3D();
	gTPCManager.DrawHelix3D();
	gTPCManager.DrawTPCHit3D();
	gTPCManager.DrawAccidental3D();
	*/
	gStyle->SetOptStat(1);
//	gTPCManager.DoCircleHough();
//	gTPCManager.DoYThetaHough();
}

void TPCEventDisplayHelix(){
	TPCEventDisplayHelix(gevnum);
	gevnum++;
}

void TPCEventDisplayLinear(int iev){
	gStyle->SetOptStat(0);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	c1->Divide(2,1);
	int xirn,xiev;
	double xim2;
	TCanvas* c2 = new TCanvas("c2","c2",500,200,600,600);
	gTPCManager.ClearHistogram();

	gTPCManager.SetEvent(iev);
	int evnum = gTPCManager.GetEvnum();
	int runnum = gTPCManager.GetRunnum();
	bool go = false;

	int nh=0;
	nh=gTPCManager.GetNhits(1);
	gTPCManager.InitializeLinear();
	for(int itr=0;itr<nh;++itr){
		gTPCManager.FillHist(itr);
	}
	cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
	c1->cd(1);
	h->Draw("col");
	gTPCManager.DrawLinear();
	gTPCManager.DrawTPCHit3D();
	c1->cd(2);
	h2->Draw("col");
	gTPCManager.DrawLinearZY();
	c2->cd();
	gTPCManager.DrawTPC();
	gTPCManager.LoadLinear3D();
	gTPCManager.DrawLinear3D();
}

void TPCEventDisplayLinear(){
	TPCEventDisplayLinear(gevnum);
	gevnum++;
}



void TPCEventDisplayAccidental(int ievt){
//	tpcfile = dir + Form("run0%d_DstTPCHelixTrackingWTarget.root",runnum);
//	tpcfile = dir + Form("test.root");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TH1D* Yhist= new TH1D("histY","histY",500,-500,500);
//	gTPCManager.SetBetheProton();	
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
	bool acc = false;
	gTPCManager.InitializeHelix();
	gTPCManager.ReconEvent();
	//		gTPCManager.InitializeAccidental();
	//		accidental_dist = gTPCManager.GetAccidentalDist();	
	nh=gTPCManager.GetNhits(0);
	for(int itr=0;itr<nh;++itr){
		gTPCManager.FillHistHits(itr);
	}
	nh=gTPCManager.GetNhits(1);
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
//	gTPCManager.DrawHelix();
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
	c1->cd(1);
	gTPCManager.DrawHelix();
	gSystem->ProcessEvents();
}
void TPCXiDisplay(int i){
	int ievt = Xievnum.at(i);
	double mxi =  XIMM.at(i);
	cout<<"MassXi = "<<mxi<<endl; 
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

void TPCKuramaDisplay(int ievt){
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TH1D* Yhist= new TH1D("histY","histY",500,-500,500);
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
//	TCanvas* c3 = new TCanvas("c3","c3",50,50,900,900);
	gTPCManager.SetEvent(ievt);
	cout<<"Event = "<<ievt<<endl;
	gTPCManager.ClearHistogram();
	int evnum = gTPCManager.GetEvnum();
	gTPCManager.SetXiStarFlag(true);
	bool go = false;
	int nh=0;
	bool acc = false;
	gTPCManager.InitializeHelix();
//	gTPCManager.ReconEvent();
	nh=gTPCManager.GetNhits(1);
	for(int itr=0;itr<nh;++itr){
//		gTPCManager.FillHistHits(itr);
		gTPCManager.FillHist(itr);
	}
	nh=gTPCManager.GetNhits(1);
	c1->cd(1);
	h->Draw("colz");
	h->SetTitle("CircleFit");
	gTPCManager.DrawHelix();
	gTPCManager.MakeHSTrack();
	gTPCManager.MakeKuramaTrack();
	gTPCManager.DrawHS();
	gTPCManager.DrawKurama();
	c1->cd(2);
	h2->Draw("colz");
	gTPCManager.DrawHSZY();
	gTPCManager.DrawKuramaZY();
//	gTPCManager.DrawHelixZY();
	c1->Modified();
	c1->Update();
/*
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
	gTPCManager.DrawHelix();
*/
	gSystem->ProcessEvents();
}
void TPCKuramaXiDisplay(int i){
	int ievt = Xievnum.at(i);
	double mxi =  XIMM.at(i);
	cout<<"MassXi = "<<mxi<<endl; 
	TPCKuramaDisplay(ievt);
}
