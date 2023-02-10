#include "TPCManager.cc"
int runnum = 5641;
//TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
int SearchPeaks(TH1D* hist,vector<double> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.15;
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

//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
//TString tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
//TString tpcfile = dir + Form("run05641_DstBeamRemover.root");
TString tpcfile = dir + Form("SelectedHelix12.root");
//TString tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
//TString tpcfile = dir + Form("run0%d_DstSelectedTPCHelixTracking.root ",runnum);
void TPCEventDisplay(){
}
void TPCEventDisplayAntiProton(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	//	gTPCManager.LoadClusterFile(tpcfile,"tree");
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	gTPCManager.SetBetheProton();	
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
	bool Single = true;
	//	Single = false;
	TFile* file = new TFile("SelectedEvents.root");
	TTree* tree = (TTree*)file->Get("tree");
	int xirn,xiev;
	double xim2;
	tree->SetBranchAddress("runnum",&xirn);
	tree->SetBranchAddress("evnum",&xiev);
	tree->SetBranchAddress("XiM2",&xim2);
	int xient = tree->GetEntries();
	for(int i=0;i<ent;++i){
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int runnum = gTPCManager.GetRunnum();
		bool go = false;
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		gTPCManager.InitializeHelix();
		gTPCManager.ReconEvent();
		//		if(!gTPCManager.XiEvent()) continue;
		if(Single){
			for(int itr=0;itr<nh;++itr){
				//				if(gTPCManager.GetClDe(itr)>60)
				//				gTPCManager.FillHist(itr);
			}
			gTPCManager.FillAntiProtonHist();
			if(h->GetEffectiveEntries()==0){
				gTPCManager.ClearHistogram();
				continue;
			}
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl; c1->cd(1);
			h->Draw("col");
			h->SetTitle("CircleFit");
			//	h->SetTitle(Form("MissMass = %f",xim2));
			//	h->SetTitle(Form("MissMass = %f",xim2));
			for(int i=0;i<gTPCManager.GetNTracks();++i){
				if(gTPCManager.IsAntiProton(i))	gTPCManager.DrawHelix(i);
			}
			//			if(gTPCManager.LambdaEvent())gTPCManager.DrawVertex();
			c1->cd(2);
			h2->Draw("col");
			for(int i=0;i<gTPCManager.GetNTracks();++i){
				if(gTPCManager.IsAntiProton(i))	gTPCManager.DrawHelixZY(i);
			}
			//			if(gTPCManager.LambdaEvent())gTPCManager.DrawVertexZY();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"Searching..."<<endl;
			gTPCManager.ClearHistogram();
		}
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
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






void TPCEventDisplayAccidental(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	int runnum = 5641;
//	tpcfile = dir + Form("run0%d_DstTPCHelixTrackingWTarget.root",runnum);
	tpcfile = dir + Form("run0%d_DstTPCHelixTracking.root",runnum);
	tpcfile = dir + Form("test.root");
	gTPCManager.LoadClusterFile(tpcfile,"tpc");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	gTPCManager.SetBetheProton();	
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);
	TCanvas* c2 = new TCanvas("c2","c2",1600,1200);
	c2->Divide(4,3);
	vector<double>* accidental_dist = new vector<double>;
	for(int i=0;i<ent;++i){
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		gTPCManager.ClearHistogram();
		int evnum = gTPCManager.GetEvnum();
		bool go = false;
		int nh=0;
//		if(!gTPCManager.TagTrig(23)) continue;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		bool acc = false;
		int iacc = 0;
		auto RV = gTPCManager.GetAccidentalR();
		auto Z0V = gTPCManager.GetAccidentalZ0();
		auto DZV = gTPCManager.GetAccidentalDZ();
		int nacc = RV.size();
		for(int ia=0; ia<nacc;++ia){
			double rr = RV[ia];
			double z0v = Z0V[ia];
			double dzv = DZV[ia];
			if(rr<1000){
				cout<<Form("%d : rad %f",iacc,rr)<<endl;
				cout<<Form("z0,dz = (%f, %f)",z0v,dzv)<<endl;
			}
			iacc++;
		}
//		if(!acc) continue;
		gTPCManager.InitializeHelix();
		gTPCManager.InitializeAccidental();
		accidental_dist = gTPCManager.GetAccidentalDist();	
		for(int itr=0;itr<nh;++itr){
			//				if(gTPCManager.GetClDe(itr)>60)
			gTPCManager.FillHist(itr);
			if( gTPCManager.GetHoughFlag()->at(itr)<990)
			cout<<"ADist : "<<accidental_dist->at(itr)<<endl;
		}
		cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
		c1->cd(1);
		h->Draw("colz");
		h->SetTitle("CircleFit");
		gTPCManager.DrawHelix();
		gTPCManager.DrawAccidental();
		c1->cd(2);
		h2->Draw("colz");
		gTPCManager.DrawHelixZY();
		gTPCManager.DrawAccidentalZY();
		c1->Modified();
		c1->Update();
		gTPCManager.FillAccHists();
		int ntAcc = gTPCManager.GetNTracksAcc();
		cout<<ntAcc<<endl;	
		for(int ic = 0 ;ic<6;++ic){
			auto hzy = gTPCManager.GetZYHistAcc(ic);
			auto hcir = gTPCManager.GetCirHistAcc(ic);
			c2->cd(2*ic+1);
			hcir->Draw("col");
			if(ic<ntAcc) gTPCManager.DrawAccidental(ic);
			c2->cd(2*ic+2);
			hzy->Draw("col");
			if(ic<ntAcc) gTPCManager.DrawAccidentalZY(ic);
		}
		c2->Modified();
		c2->Update();
		cout<<"Event "<<endl;
		gSystem->ProcessEvents();
		cin.ignore();
		cout<<"Searching..."<<endl;
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
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

