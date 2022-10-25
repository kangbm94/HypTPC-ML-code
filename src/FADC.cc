#include "../include/FADC.hh"
#include "TPCManager.cc"
int flayers[10]={31,31,31,31,31,31,31,31,31,31};
int frows[10]={60,61,62,63,64,65,66,67,68,69};
void FADC(){
	int runnum=5754;
	TString dir = base_dir+"MayRun/rootfiles/FADC/";
	TString filename = Form("run0%d_TPCWaveform.root",runnum);
	FADCManager FMan;
	FMan.LoadFile(dir+filename);
	

	TString tpcdir = base_dir+"MayRun/rootfiles/Defocus/";
	TString tpcfilename = Form("run0%d_DstTPCBcOutOld.root",runnum);
	gTPCManager.LoadFile(tpcdir+tpcfilename);
	gTPCManager.LoadBcOut();
	gTPCManager.LoadClusterChain();
	gTPCManager.InitializeHistograms();
	
	TCanvas* c1 = new TCanvas("c2","c2",1500,1200);
	c1->Divide(2,1);
	TPad* pad1 = new TPad("pad1","pad1",0,0,750,600);
	TPad* pad2 = new TPad("pad2","pad2",750,0,1500,600);
	pad2->Divide(3,4);
	int ent = FMan.GetEntries();
	int ent2 = gTPCManager.GetEntries();
	TGraph* BcGr ;
	TF1* Bcf = new TF1("Bcf","pol1",-300,300);	
	if(ent!=ent2) cout<<"Warning! event number does not match!"<<endl;
	auto* h = gTPCManager.GetPadHistogram();

	for(int i=0;i<ent;++i){	
		//TPCDisplay	
		gTPCManager.SetEvent(i);
		int nh = gTPCManager.GetNhits(0);
		int ncl = gTPCManager.GetNhits(1);
		if(gTPCManager.GetBCnt()!=1||ncl>35){
			continue;
		}
		Track Bc = gTPCManager.GetTrack();
		double bcx=Bc.GetPosition(-K18HS).X();
		int cnt =0;
		cout<<"FillingHist"<<endl;
		for(int j=0;j<nh;++j){
			gTPCManager.FillHist(j);
		}
		double bcpx[100]={0},bcpz[100]={0};
		cout<<ncl<<endl;
		for(int j=0;j<ncl;++j){
			TVector3 ClPos = gTPCManager.GetClusterPosition(j);
			double clx=ClPos.X(),clz=ClPos.Z();
			double bcx=Bc.GetPosition(clz).X();
			bcpx[j]=bcx;bcpz[j]=clz;
			int l,r;
			int bcpad = tpc::findPadID(clz,bcx);
			GetLayerRow(bcpad,l,r);
			cout<<Form("BcPos: (%f,%f)",bcpz[j],bcpx[j])<<endl;
		}
		BcGr = new TGraph(ncl,bcpz,bcpx);
		pad1->cd();	
		cout<<"Hist"<<endl;
		h->SetTitle(Form("Run05754,Event%d",i));
		h->Draw("colz");
		h->GetXaxis()->SetRangeUser(220,250);
		h->GetYaxis()->SetRangeUser(20,90);
		gTPCManager.SetPadContent(31,59,2);
		gTPCManager.SetPadContent(31,70,2);
		cout<<"BCGR"<<endl;
		BcGr->Draw("Psame");
		BcGr->SetMarkerStyle(20);
		BcGr->Fit("Bcf");
		c1->cd(1);
		pad1->Draw();
	
		cout<<"BCTrack"<<endl;

		//WaveformDisplay	
		FMan.SetEvent(i);
		FMan.LoadWaveform();
		int wn = FMan.GetWaveNum();
		bool flag = false;
		double par[10];
		TH1D* wh[11];TH1D* rh[11];
		for(int j=0;j<Min(wn,10);++j){
			wh[j] = FMan.GetWaveformHist(flayers[j],frows[j]);
			rh[j] = FMan.GetRawWaveformHist(flayers[j],frows[j]);
			if(wh[j]){
				pad2->cd(j+2);
				FMan.DoFit(flayers[j],frows[j],par);	
				wh[j]->Draw();
				rh[j]->SetLineColor(kGreen);
				rh[j]->Draw("same");
				rh[j]->SetTitle(Form("LR=(%d,%d),ev=%d",flayers[j],frows[j],i));
			}
			else{
				pad2->cd(j+1);
			//	h[j]= new TH1D(Form("dum%d_%d",j,i),Form("dum%d_%d",j,i),10,0,10);
				rh[j]= new TH1D(Form("rdum%d_%d",j,i),Form("rdum%d_%d",j,i),10,0,10);
				rh[j]->Draw("col");
			}
		}
		c1->cd(2);
		pad2->Draw();





		pad1->Modified();
		pad1->Update();
		pad2->Modified();
		pad2->Update();
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
		gTPCManager.ClearHistogram();
		/*
		c1->cd(9);
		auto hb = FMan.GetBaseline();
		auto rhb = FMan.GetRawBaseline();

		if(!hb){
				hb=new TH1D(Form("dum%d",i),Form("dum%d",i),10,0,10);
				rhb=new TH1D(Form("rdum%d",i),Form("rdum%d",i),10,0,10);
		}
		if(!FMan.Check()){
			hb->Draw("col");
			rhb->SetLineColor(kGreen);
			rhb->Draw("same");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
		}
		*/
		//	delete hb;delete rhb;
			for(int j=0;j<8;++j){
//				if(h[j]) delete h[j];
	//			if(rh[j]) delete rh[j];
			}
//			delete h;delete rh;delete hb;delete rhb;
//		c`->Reset("ICES");
			FMan.Clear();
	}
}
