#include "TPCManager.cc"
int runnum = 5000;
TString dir = base_dir+"MayRun/rootfiles/CH2/TPC/";
//TString tpcfile = dir + Form("run0%d_TPC_RMSCut.root",runnum);
TString tpcfile = dir + Form("SelectedHelix.root");
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






void TPCEventDisplay(){
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
		if(i%1000==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int runnum = gTPCManager.GetRunnum();
		bool go = false;
		for(int j=0;j<xient;j++){
			tree->GetEntry(j);
			if(runnum == xirn and evnum == xiev and abs(xim2-1314)<0.05){cout<<xirn<<endl; go = false;break;}
			else{
				go = true;
			}
		}
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
				double y = pos.Y();
				if(abs(y)<30)continue;
				else h3->Fill(y);
			}
			c2->cd(1);
			int np = SearchPeaks(h3,peaks); 
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
			//			gTPCManager.DrawHelixZY();
			//			gTPCManager.DrawVertexZY();

			fpol->SetParLimits(0,-100,100);
			fpol->SetParLimits(1,-100,100);
			fpol->SetParLimits(2,-1e-4,-1e-7);
			fpol->SetParameter(2,-7.6e-5);
			bool BeamFlg = false;
			for(int ipk=0;ipk<np;++ipk){
				gHitPos.clear();
				gRes.clear();
				if(abs(peaks[ipk])<60) continue;
				gTPCManager.FillgHitPos(peaks[ipk],width);	
				if(gHitPos.size()<10) continue;
				cout<<"Peak: "<<peaks[ipk]<<endl;
				double p0,p1,p2;
				for(auto g:gHitPos){
					hist_beam->Fill(g.Z(),g.X());
				}
				hist_beam->Fit("fpol2","QR0");
				hist_beam->Reset("ICES");
				p0=fpol->GetParameter(0);
				p1=fpol->GetParameter(1);
				p2=fpol->GetParameter(2);
				for(auto g:gHitPos){
					double xx[2]={g.Z(),0};
					double par[3]={p0,p1,p2};
					double val = Quadratic(xx,par);
					if(abs(val-g.X())<20)hist_beam->Fill(g.Z(),g.X());
				}
				hist_beam->Fit("fpol2","QR0");
				hist_beam->Reset("ICES");
				p0=fpol->GetParameter(0);
				p1=fpol->GetParameter(1);
				p2=fpol->GetParameter(2);
				for(auto g:gHitPos){
					double xx[2]={g.Z(),0};
					double par[3]={p0,p1,p2};
					double val = Quadratic(xx,par);
					if(abs(val-g.X())<10)hist_beam->Fill(g.Z(),g.X());
				}
				hist_beam->Fit("fpol2","QR0");
				hist_beam->Reset("ICES");
				p0=fpol->GetParameter(0);
				p1=fpol->GetParameter(1);
				p2=fpol->GetParameter(2);
				double min_z = 9999;
				for(auto g:gHitPos){
					double xx[2]={g.Z(),0};
					double par[3]={p0,p1,p2};
					double val = Quadratic(xx,par);
					if(abs(val-g.X())<5){
						hist_beam->Fill(g.Z(),g.X());
						min_z = min(g.Z(),min_z);
					}
				}
				if(hist_beam->GetEntries()>10 and min_z<-200){
					BeamFlg=true;
					hist_beam->Fit("fpol2","R0");
					p0=fpol->GetParameter(0);
					p1=fpol->GetParameter(1);
					p2=fpol->GetParameter(2);
					par0.push_back(p0);
					par1.push_back(p1);
					par2.push_back(p2);
					SelPeaks.push_back(peaks[ipk]);
				}
			}
			c2->cd(2);
			for(auto g:gHitPos){
				hist_beam->Fill(g.Z(),g.X());
			}
			double bend = fpol->GetParameter(2);
			double rad = 1/(2*bend);
			double mom = rad*(Const*dMagneticField);
			cout<<"Mom: "<<mom<<endl;
			c2->Modified();
			c2->Update();
			for(int ipk=0;ipk<SelPeaks.size();++ipk){
				double  y = SelPeaks[ipk];
				c1->cd(1);
				TString title = Form("pol2%d",ipk);
				TF1* fquad = new TF1(title,"Quadratic",-250,250,3);
				fquad->SetParameter(0,par0[ipk]);
				fquad->SetParameter(1,par1[ipk]);
				fquad->SetParameter(2,par2[ipk]);
				fquad->SetLineWidth(2);
				fquad->SetLineColor(ipk+1);
				fquad->Draw("same");
				c1->cd(2);
				TLine* line = new TLine(-250,y,250,y);
				line->SetLineWidth(4);
				line->SetLineColor(ipk+1);
				line->Draw("same");
			}
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"Searching..."<<endl;
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
