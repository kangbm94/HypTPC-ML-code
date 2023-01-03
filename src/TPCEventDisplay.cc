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
void TPCEventDisplayRaw(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",800,1200);
	c1->Divide(2,1);
	bool Single = true;
//	Single = false;
	cout<<ent<<endl;
	for(int i=0;i<ent;++i){
		if(i%1==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=0;
		nh=gTPCManager.GetNhits(1);
//		if(!nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		if(Single){
			cout<<"Drawing..."<<i<<endl;
			c1->cd(1);
			h->Draw("colz");
			h->SetTitle(Form("Event%d",evnum));
			c1->cd(2);
			h2->Draw("colz");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
		}
	}
//	h->Draw("colz");
//	h->SetTitle(Form("Run%d",runnum));
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
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
			c1->cd(1);
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
		gTPCManager.InitializeHelix();
		gTPCManager.ReconEvent();
//		if(!gTPCManager.XiEvent()) continue;
		if(Single){
			for(int itr=0;itr<nh;++itr){
//				if(gTPCManager.GetClDe(itr)>60)
				gTPCManager.FillHist(itr);
			}
//			gTPCManager.FillAntiProtonHist();
			if(h->GetEffectiveEntries()==0){
				gTPCManager.ClearHistogram();
				continue;
			}
			cout<<Form("Drawing Event (%d,%d)",runnum,evnum)<<endl;
			c1->cd(1);
			h->Draw("col");
			h->SetTitle("CircleFit");
			gTPCManager.DrawHelix();
		//	h->SetTitle(Form("MissMass = %f",xim2));
		//	h->SetTitle(Form("MissMass = %f",xim2));
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
	}
	h->Draw("colz");
	h->SetTitle(Form("Run%d",runnum));
}







void DisplayBeamRemover(){
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	auto hf = gTPCManager.GetPadHistogramf();
	auto hf2 = gTPCManager.GetZYHistogramf();
	auto hb = gTPCManager.GetPadHistogramb();
	auto hb2 = gTPCManager.GetZYHistogramb();
	int ent = gTPCManager.GetEntries();	
	TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
	c1->Divide(2,2);
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,1);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	TH1D* hist = new TH1D("bhy","bhy",140,-350,350);
	bool Single = true;
//	Single = false;
	cout<<ent<<endl;
	int start = 0;
	for(int i=start;i<ent;++i){
		if(i%1==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=gTPCManager.GetNhits();
		int nhf=gTPCManager.GetNhitsf();
		int nhb=gTPCManager.GetNhitsb();
//		if(!nhf or !nhb) continue;
//		if(nhb<80) continue;
//		if(nhf+ nhb ==nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		for(int itr=0;itr<nh;++itr){
			auto pos = gTPCManager.GetPosition(itr);
			double y = pos.Y();
			if(abs(y)>50) hist->Fill(y);
		}
		for(int itr=0;itr<nhf;++itr){
			gTPCManager.FillHistf(itr);
		}
		for(int itr=0;itr<nhb;++itr){
			gTPCManager.FillHistb(itr);
		}
		if(Single){
			cout<<"Drawing..."<<i<<endl;
			c1->cd(1);
			hf->Draw("colz");
			hf->SetTitle(Form("Event%d",evnum));
			c1->cd(2);
			hf2->Draw("colz");
			c1->cd(3);
			hb->Draw("colz");
			auto by = gTPCManager.GetBeamY();
			auto bv = gTPCManager.GetBeamV();
			auto p0 = gTPCManager.GetBeamP0();
			auto p1 = gTPCManager.GetBeamP1();
			auto p2 = gTPCManager.GetBeamP2();
			for(int i=0;i<by.size();++i){
				TF1* f = new TF1(Form("f%d",i),"pol2",-250,250);
				f->SetParameter(0,p0[i]);
				f->SetParameter(1,p1[i]); f->SetParameter(2,p2[i]);
				f->SetLineColor(i+2);
				f->SetLineWidth(3);
				f->Draw("same");
			}
			hb->SetTitle(Form("Event%d",evnum));
			c1->cd(4);
			hb2->Draw("colz");
			for(int i=0;i<by.size();++i){
				TF1* f2 = new TF1(Form("f2%d",i),"pol1",-250,250);
				f2->SetParameter(0,by[i]);
				f2->SetParameter(1,bv[i]);
				f2->SetLineColor(i+2);
				f2->SetLineWidth(3);
				f2->Draw("same");
			}
			c1->Modified();
			c1->Update();
			c2->cd(1);
			h->Draw("colz");
			h->SetTitle(Form("Event%d",evnum));
			for(int i=0;i<by.size();++i){
				TF1* f = new TF1(Form("f%d",i),"pol2",-250,250);
				f->SetParameter(0,p0[i]);
				f->SetParameter(1,p1[i]);
				f->SetParameter(2,p2[i]);
				f->SetLineColor(i+2);
				f->SetLineWidth(3);
				f->Draw("same");
			}
			c2->cd(2);
			h2->Draw("colz");
			for(int i=0;i<by.size();++i){
				TF1* f2 = new TF1(Form("f2%d",i),"pol1",-250,250);
				f2->SetParameter(0,by[i]);
				f2->SetParameter(1,bv[i]);
				f2->SetLineColor(i+2);
				f2->SetLineWidth(3);
				f2->Draw("same");
			}
			
			c2->Modified();
			c2->Update();
			c3->cd();
			hist->Draw();
			vector<double> dum;
			SearchPeaks(hist,dum);
			c3->Modified();
			c3->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
			hist->Reset("ice");
		}
	}
//	h->Draw("colz");
//	h->SetTitle(Form("Run%d",runnum));
}
void DisplayBeamFit(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	gTPCManager.LoadClusterFile(tpcfile);
	auto h = gTPCManager.GetPadHistogram();
	auto h2 = gTPCManager.GetZYHistogram();
	int ent = gTPCManager.GetEntries();	
	int nc2 = 6;
	TH1D* hist = new TH1D("bhy","bhy",140,-350,350);
	TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
	c1->Divide(2,1);
	TCanvas* c2 = new TCanvas("c2","c2",1200,1200);
	c2->Divide(3,2);
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	bool Single = true;
	TH2D* h_zx[9];
	TH2D* h_zxraw[9];
	TH2D* h_zy[9];
	TF1* quad[9];
	TF1* lin[9];
	for(int i=0;i<nc2;++i){
		TString title = Form("HistZX_%d",i);
		TString titler = Form("HistZXraw_%d",i);
		TString quadt = Form("quad_%d",i);
		TString titley = Form("HistZX_%d",i);
		TString lint = Form("lin_%d",i);
		h_zx[i]=new TH2D(title,title,100,-250,250,120,-120,120);
		h_zxraw[i]=new TH2D(titler,titler,100,-250,250,120,-120,120);
		quad[i] = new TF1(quadt,"pol2",-250,250); 
		h_zy[i]=new TH2D(titley,titley,100,-250,250,700,-350,350);
		lin[i] = new TF1(lint,"pol1",-250,250); 
	}
//	Single = false;
	cout<<ent<<endl;
	int start = 0;
	double y_min = -50,y_max=50;
	for(int i=start;i<ent;++i){
		if(i%1==0) cout<<"Event "<<i<<endl;
		gTPCManager.SetEvent(i);
		int evnum = gTPCManager.GetEvnum();
		int nh=gTPCManager.GetNhits();
//		if(!nhf or !nhb) continue;
//		if(nhb<80) continue;
//		if(nhf+ nhb ==nh) continue;
		for(int itr=0;itr<nh;++itr){
			gTPCManager.FillHist(itr);
		}
		vector<TVector3> posCont;
		posCont.clear();
		for(int itr=0;itr<nh;++itr){
			auto pos = gTPCManager.GetPosition(itr);
			double y = pos.Y();
			if(y_max<y or y<y_min){
				hist->Fill(y);
				posCont.push_back(pos);
			}
		}
		if(Single){
			cout<<"Drawing..."<<i<<endl;
			c3->cd();
			hist->Draw();
			vector<double> peaks;
			int nhc = posCont.size();
			int npeaks = SearchPeaks(hist,peaks);
			if(npeaks>nc2) npeaks=nc2;
			double Ywidth=10;
			double dist_cut[3]={20,10,5};
			double dist_cutY[3]={10,7.5,5};
			for(int ih=0;ih<nhc;++ih){
				auto pos = posCont.at(ih);
				double x = pos.X(),y=pos.Y(),z=pos.Z();
				for(int ip = 0 ; ip<npeaks;++ip){
					double peak = peaks.at(ip);
					if(y_min<y and y<y_max) continue;
					if(abs(y-peak)>Ywidth) continue;
					h_zx[ip]->Fill(z,x);
					h_zxraw[ip]->Fill(z,x);
					h_zy[ip]->Fill(z,y);
				}
			}
			double min_z[9];
			vector<double> bp0;
			vector<double> bp1;
			vector<double> bp2;
			vector<double> by;
			vector<double> bv;
			for(int ip=0;ip<npeaks;++ip){
				double peak = peaks.at(ip);
				cout<<"Peak: "<<peak<<" -> ";
				TString quadt = Form("quad_%d",ip);
				TString lint = Form("lin_%d",ip);
				double p0,p1,p2,v;
				for(int itr=0;itr<3;++itr){
					min_z[ip]=9999;
					if(h_zx[ip]->GetEffectiveEntries()<5)continue;
					lin[ip]->SetParameter(0,peak);
					lin[ip]->SetParLimits(1,-0.1,0.1);
					quad[ip]->SetParameter(2,-7.6e-5);//~1.8GeV;
					//quad[ip]->SetParLimits(0,-50,50);
					quad[ip]->SetParLimits(1,-0.1,-0.01);
					quad[ip]->SetParLimits(2,-3e-4,-3e-5);
					h_zx[ip]->Fit(quadt,"QR0B");
					h_zy[ip]->Fit(lint,"QR0B");
					p0 = quad[ip]->GetParameter(0);
					p1 = quad[ip]->GetParameter(1);
					p2 = quad[ip]->GetParameter(2);
					peak = lin[ip]->GetParameter(0);
					v = lin[ip]->GetParameter(1);
					h_zx[ip]->Reset("ices");	
					h_zy[ip]->Reset("ices");	
					for(auto pos:posCont){
						double x = pos.X(),y=pos.Y(),z=pos.Z();
						if(abs(peak+v*z-y)>dist_cutY[itr] ) continue;
						double val = p0+p1*z+p2*z*z;
						if(abs(val-x)>dist_cut[itr])continue;
						h_zx[ip]->Fill(z,x);
						h_zy[ip]->Fill(z,y);
						min_z[ip]=min(z,min_z[ip]);
					}
				}
				c2->cd(ip+1);
//				h_zxraw[ip]->Draw("colz");
				h_zx[ip]->Draw("colz");
				h_zx[ip]->Fit(quadt,"QRB");
	//			quad[ip]->Draw("same");
				p0 = quad[ip]->GetParameter(0);
				p1 = quad[ip]->GetParameter(1);
				p2 = quad[ip]->GetParameter(2);
				peak = lin[ip]->GetParameter(0);
				v = lin[ip]->GetParameter(1);
				cout<<h_zx[ip]->GetTitle()<<" : Pos: "<<peak<<endl;
				double mom = (1./(2*p2))*(0.299792458)*0.9;
				cout<<Form("p1, ent , mom ,minz = (%f,%f,%f,%f)",p1,h_zx[ip]->GetEntries(),mom,min_z[ip])<<endl;
				if(!(h_zx[ip]->GetEntries()>10  and abs(mom)>500 and min_z[ip]<-150))continue;	
				bp0.push_back(p0);
				bp1.push_back(p1);
				bp2.push_back(p2);
				bv.push_back(v);
				by.push_back(peak);
			}
			c1->cd(1);
			h->Draw("colz");
			h->SetTitle(Form("Event%d",evnum));
			int cnt=0;
			for(int i=0;i<nc2;++i){
				if(cnt<by.size()){	
				quad[i]->SetParameter(0,bp0[cnt]);
				quad[i]->SetParameter(1,bp1[cnt]);
				quad[i]->SetParameter(2,bp2[cnt]);
				quad[i]->Draw("same");
				cnt++;
				}
				quad[i]->SetLineColor(i+1);
				quad[i]->SetLineWidth(3);
			}
			c1->cd(2);
			h2->Draw("colz");
			cnt=0;
			for(int i=0;i<nc2;++i){
				if(cnt<by.size()){
					lin[i]->SetParameter(0,by[cnt]);
					lin[i]->SetParameter(1,bv[cnt]);
					lin[i]->Draw("same");
					cnt++;
				}
				lin[i]->SetLineColor(i+1);
				lin[i]->SetLineWidth(3);
			}
			c1->Modified();
			c1->Update();
			c2->Modified();
			c2->Update();
			c3->Modified();
			c3->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			gTPCManager.ClearHistogram();
			hist->Reset("ice");
			c2->Clear("D");
			for(int ih=0;ih<nc2;++ih){
					h_zx[ih]->Reset("ice");
					h_zxraw[ih]->Reset("ice");
					h_zy[ih]->Reset("ice");
			}
		}
	}
//	h->Draw("colz");
//	h->SetTitle(Form("Run%d",runnum));
}
