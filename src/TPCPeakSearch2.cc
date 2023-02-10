#include "TPCManager.cc"
#include "BeamRemover.cc"
int runnum = 5000;
TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
TString tpcfile = base_dir + Form("SelectedHelix12.root");
//TString tpcfile = base_dir + Form("run05146_DstTPCHelixTracking.root");
//TString tpcfile = dir + Form("run0%d_DstSelectedTPCHelixTracking.root ",runnum);
int SearchPeaks(TH1D* hist,vector<double> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.2;
	//	int npeaks = spec.Search(hist,sig,"goff",th);
	int npeaks = spec.Search(hist,sig,"",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}

const double Const = 0.299792458;
const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

double Quadratic(double* x, double* p){
	return p[0]+p[1]*x[0]+p[2]*x[0]*x[0];
}






void TPCPeakSearch2(){
}
void TPCEventDisplayHelix(){
	TF1* fpol = new TF1("fpol2","Quadratic",-250,250,3);
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile,"tree");
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	auto h3 = gTPCManager.GetYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",1000,700);
	c1->Divide(2,1);
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(2,2);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
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
	TH2D* hist_beam = new TH2D("bh","bh",100,-250,250,100,-100,100);
	for(int i=0;i<ent;++i){
		BeamRemover BR(-30,30,-50,50);
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int runnum = gTPCManager.GetRunnum();
		bool go = false;
		/*
		for(int j=0;j<xient;j++){
			tree->GetEntry(j);
			if(runnum == xirn and evnum == xiev and abs(xim2-1314)<0.05){cout<<xirn<<endl; go = false;break;}
			else{
				go = true;
			}
		}
		*/
		//		if(go) continue;
		int nh=0;
		nh=gTPCManager.GetNhits(1);
		if(!nh) continue;
		//		gTPCManager.ReconEvent();

		//		if(!gTPCManager.XiEvent()) continue;
		gTPCManager.InitializeHelix();
		vector<double>peaks;
		vector<double>SelPeaks;
		vector<double>par0;
		vector<double>par1;
		vector<double>par2;
		vector<TH2D> hists;
		if(Single){
			for(int itr=0;itr<nh;++itr){
				gTPCManager.FillHist(itr);
			}
			for(int i=0;i<nh;++i){
				auto pos = gTPCManager.GetPosition(i);
				BR.LoadHit(pos);
				double y = pos.Y();
				if(abs(y)<30)continue;
				else h3->Fill(y);
			}
			auto ydistr = BR.GetYDistrib();
			c3->cd();
			ydistr->Draw();
			int np = BR.SearchPeaks(ydistr,peaks);
			int height = ydistr->GetMaximum();
			double yw = BR.GetYwidth();
			for(auto peak : peaks){
				TLine* left = new TLine(peak-yw,0,peak-yw,height);
				TLine* right = new TLine(peak+yw,0,peak+yw,height);
				left->SetLineWidth(2);
				left->SetLineColor(kRed);
				right->SetLineWidth(2);
				right->SetLineColor(kRed);
				left->Draw("same");
				right->Draw("same");
			}
			c3->Modified();
			c3->Update();
			BR.AssignHits(peaks);
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
			double width =  5;
			double arglist[5];
			arglist[0]=5.89;
			TString name[5] = {"cx", "cy", "z0", "r", "dz"};
			double err[5]={-999.,-999.,-999.,-999.,-999.};
			int Err;
			double bnd1,bnd2;
			double cx,cy,z0,r,dz;
			double par[5] = {6000,-6000,0,0,0};
			Int_t ierflg = 0;
			c1->cd(1);
			h->Draw("colz");
			h->SetTitle("CircleFit");
			//	h->SetTitle(Form("MissMass = %f",xim2));
			//	h->SetTitle(Form("MissMass = %f",xim2));
			//			gTPCManager.DrawHelix();
			//			gTPCManager.DrawVertex();
			c1->cd(2);
			h2->Draw("col");
			c1->Modified();
			c1->Update();
			//			gTPCManager.DrawHelixZY();
			//			gTPCManager.DrawVertexZY();

			fpol->SetParLimits(0,-100,100);
			fpol->SetParLimits(1,-100,100);
			fpol->SetParLimits(2,-1e-4,-1e-7);
			fpol->SetParameter(2,-7.6e-5);
			bool BeamFlg = false;
			for(int ipk=0;ipk<np;++ipk){
				cout<<"Peak: "<<peaks[ipk]<<endl;
				BR.DoCircleHough(ipk);
				BR.DoYThetaHough(ipk);
				int col = 0;
				TLine* up = new TLine(-250,peaks[ipk]+yw,250,peaks[ipk]+yw);
				TLine* dn = new TLine(-250,peaks[ipk]-yw,250,peaks[ipk]-yw);
				c1->cd(2);
				up->SetLineWidth(2);
				dn->SetLineWidth(2);
				up->Draw("same");
				dn->Draw("same");
				c1->Modified();
				c1->Update();
				c2->cd();
				auto h = BR.GetHoughHist();
				h->Draw("colz");
				auto cirs = BR.GetCircle(ipk);
				c2->Modified();
				c2->Update();
				gSystem->ProcessEvents();
				cin.ignore();
				for(auto cir:cirs){
					col++;
					cir->SetLineColor(col);
					cir->SetLineWidth(2);
					cir->Draw("same");
					cir->Draw("same");
					c2->Modified();
					c2->Update();
					gSystem->ProcessEvents();
					cin.ignore();
				}
				cout<<"NBeam : "<<col<<endl;
				h->Reset();
				for(auto cir:cirs){
					delete cir;
				}
				delete up;
				delete dn;
				cout<<"Searching..."<<endl;
			}
			for(auto g:gHitPos){
				hist_beam->Fill(g.Z(),g.X());
			}
			double bend = fpol->GetParameter(2);
			double rad = 1/(2*bend);
			double mom = rad*(Const*dMagneticField);
			cout<<"Mom: "<<mom<<endl;
			
			hist_beam->Reset("ICES");
			gTPCManager.ClearHistogram();
		}
	}
	//	h->Draw("colz");
	//	h->SetTitle(Form("Run%d",runnum));
}
void DisplayY(){
	gTPCManager.LoadClusterFile(tpcfile);
	auto h2 = gTPCManager.GetZYHistogram();
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,1);



}
