#include "TPCManager.hh"

const int max_ntrk = 16;
TPCManager T;
void TPC(){
	T.InitializeHistograms();
	cout<<"ViewTPC(TString Filename)"<<endl;
	cout<<"TagTPC(TString Filename)"<<endl;
} void ViewTPC(TString Filename){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	int entries = T.GetNEvent();
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);	
	c1->Divide(3,1);
	for(int i=0;i<entries;++i){
		int ThisEvent = 0;
		T.SetEvent(i);
		if(i%1000==0)cout<<i<<endl;
		ThisEvent=T.WhichEvent();
		if(ThisEvent==Else){
			//			cout<<"Else"<<endl;
		}
		else if(ThisEvent==L2PPi){
			//			cout<<"Lambda -> P Pi"<<endl;
		}
		else if(ThisEvent==L2NPi){
			//			cout<<"Lambda -> N Pi"<<endl;
		} for(int j=0;j<T.GetNpadG4();++j){
			int padID= T.GetPadIDG4(j);
			TVector3 vec = T.GetG4Position(j);
			double x = vec.X();double z = vec.Z();
			double x_t = tpc::getPosition(padID).X();
			double z_t = tpc::getPosition(padID).Z();
			T.FillFlatHist(padID);
			T.FillHist(z_t,x_t);
		}
		if(ThisEvent==L2PPi||ThisEvent==L2NPi){
			TString title = Form("evt=%d",i);
			T.SetTitle(title);
			cout<<"wait"<<endl;
			c1->cd(1);
			T.DrawHist();
			c1->cd(2);
			T.DrawPosHist();
			c1->cd(3);
			T.DrawFlatHist();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
		}
		if(i%1==0){
			T.ClearHistogram();
		}
	}
	cout<<"End"<<endl;
}

void TagTPC(TString Filename){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	int entries = T.GetNEvent();
	TFile* outfile = new TFile("TrainDataTagged_wo_bg.root","recreate");
	TTree* outtree = new TTree("tree","tree");
	int evt_per_tag = 10000;
	const int max_trk=1000;
	short x[max_trk];
	short y[max_trk];
	short z[max_trk];
	int TPCEventTag=0,ntrk=0;
	//	outtree->Branch("TPCEvent",TPCEvent,Form("TPCEvent[%d][%d][%d]/S",nbin,nbin,depth));
	outtree->Branch("TPCEventTag",&TPCEventTag,"TPCEventTag/I");
	outtree->Branch("ntrk",&ntrk,"ntrk/I");
	outtree->Branch("x",x,"x[ntrk]/S");
	outtree->Branch("y",y,"y[ntrk]/S");
	outtree->Branch("z",z,"z[ntrk]/S");
	
	int TagNumber[10] = {0};

	for(int i=0;i<entries;++i){
		int ThisEvent = 0;
		T.SetEvent(i);
		ntrk=T.GetNpadG4();
		if(i%10000==0)cout<<i<<endl;
		ThisEvent=T.WhichEvent();
		for(int j=0;j<max_ntrk;++j){
			x[j]=0;y[j]=0;z[j]=0;
		}
		T.AssignEvent(x,y,z);
		TPCEventTag=ThisEvent;
		if(TagNumber[ThisEvent]<evt_per_tag&&ThisEvent!=Else){
			outtree->Fill();
		}
		TagNumber[ThisEvent]++;
	}
	outfile->Write();
}
