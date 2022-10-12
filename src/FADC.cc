#include "../include/FADC.hh"
#include "TPCManager.cc"
int flayers[8]={31,31,31,31,19,20,21,22};
int frows[8]={63,65,67,69,129,150,129,150};
void FADC(){
	int runnum=5754;
	TString dir = base_dir+"MayRun/rootfiles/FADC/";
	TString filename = Form("run0%d_TPCWaveform.root",runnum);
	FADCManager FMan;
	FMan.LoadFile(dir+filename);
	int ent = FMan.GetEntries();
	for(int i=114;i<ent;++i){	
//	for(int i=0;i<ent;++i){	
		FMan.SetEvent(i);
		FMan.LoadWaveform();
		int wn = FMan.GetWaveNum();
		bool flag = false;
		double par[10];
		TCanvas* c1 = new TCanvas("c2","c2",1500,1200);
		c1->Divide(3,3);
		TH1D* h[9];TH1D* rh[9];
		for(int j=0;j<Min(wn,8);++j){
			h[j] = FMan.GetWaveformHist(flayers[j],frows[j]);
			rh[j] = FMan.GetRawWaveformHist(flayers[j],frows[j]);
			if(h[j]){
				c1->cd(j+1);
				FMan.DoFit(flayers[j],frows[j],par);	
				h[j]->Draw();
				rh[j]->SetLineColor(kGreen);
				rh[j]->Draw("same");
			}
			else{
				c1->cd(j+1);
			//	h[j]= new TH1D(Form("dum%d_%d",j,i),Form("dum%d_%d",j,i),10,0,10);
				rh[j]= new TH1D(Form("rdum%d_%d",j,i),Form("rdum%d_%d",j,i),10,0,10);
				rh[j]->Draw("col");
			}
		}
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
		//	delete hb;delete rhb;
			for(int j=0;j<8;++j){
//				if(h[j]) delete h[j];
	//			if(rh[j]) delete rh[j];
			}
			delete	c1;
//			delete h;delete rh;delete hb;delete rhb;
//		c`->Reset("ICES");
			FMan.Clear();
	}
}
