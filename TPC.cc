#include "TPCManager.hh"

const int max_ntrk = 16;
TPCManager T;
void TPC(){
	T.InitializeHistograms();
	cout<<"ViewTPC(TString Filename)"<<endl;
	cout<<"TagTPCEvents(TString Filename, TString OutFileName)"<<endl;
	cout<<"TagTPCTracks(TString Filename)"<<endl;
	cout<<"TestML(TString Filename)"<<endl;
	cout<<"ConvertRealTPC(TString Filename)"<<endl;
	cout<<"CheckRealTPCTraining(TString Filename)"<<endl;
	gStyle->SetPalette(22);
}
void ConvertRealTPC(TString Filename){
	TString dir = "../MayRun/rootfiles/CH2/TPC";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadFile(TPCFile);
	int ent =T.GetNEvent();
	TFile* outfile = new TFile("RealData05641.root","recreate");
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];
	int nhtpc=0,evnum=0;
	TTree* outtree = new TTree("tree","tree");
	outtree->Branch("evnum",&evnum,"evnum/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");

	for(int i=0;i<ent;++i){
		if(i%1000==0)cout<<Form("Processing %d th event...",i)<<endl;
		T.SetEvent(i);	
		evnum=i;
		nhtpc =Min(T.GetNpad(), max_nh);
		T.AssignRealEvent(x,y,z,dedx);
		outtree->Fill();
	}
	outfile->Write();
}
void CheckRealTPCTraining(TString Filename,TString PredFile){
	TFile* datafile = new TFile(Filename,"read");
	TTree* datatree = (TTree*)datafile->Get("tree");

	int nhtpc;
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	int ntrk;
	datatree->SetBranchAddress("nhtpc",&nhtpc);
	datatree->SetBranchAddress("x",x);
	datatree->SetBranchAddress("y",y);
	datatree->SetBranchAddress("z",z);
//	datatree->SetBranchAddress("ntrk",&ntrk);

	TFile* predfile = new TFile(PredFile,"read");
	TTree* predtree = (TTree*)predfile->Get("tree");
	int pred;
	int tag;
	double prb[10];
	int evnum;
	predtree->SetBranchAddress("evnum",&evnum);
	predtree->SetBranchAddress("pred",&pred);
	predtree->SetBranchAddress("prb",prb);
	predtree->SetBranchAddress("tag",&tag);

	TH2D* hist = new TH2D("hist","zx",128,-250,250,128,-300,350);
	TH2D* hist2 = new TH2D("hist2","zy",128,-250,250,128,-300,350);
	TH2D* hist3 = new TH2D("hist3","xy",128,-250,250,128,-300,350);
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(3,1);
	TCanvas* c2 = new TCanvas("c2","c2",1500,800);
	c2->Divide(3,1);
	TH1I* nhhist1 = new TH1I("h1","nh=1",100,0,300);
	TH1I* nhhist2 = new TH1I("h2","nh=2",100,0,300);
	TH1I* nhhist3 = new TH1I("h3","nh=3",100,0,300);
	int ent = predtree->GetEntries();
	for(int i=0;i<ent;++i){
		int evt = i;
		predtree->GetEntry(evt);
		datatree->GetEntry(evnum);
		for(int j=0;j<nhtpc;++j){
//			T.FillHist(z[j],x[j]);
			hist->Fill(z[j],x[j]);
			hist2->Fill(z[j],y[j]);
			hist3->Fill(x[j],y[j]);
		}
		if(nhtpc<50){
			cout<<"Number of Hits: "<<nhtpc<<endl;
			cout<<"Predicted Number of Tracks: "<<pred<<endl;
			c1->cd(1);
			hist->Draw("colz");
			c1->cd(2);
			hist2->Draw("colz");
			c1->cd(3);
			hist3->Draw("colz");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
		}
		if(pred==1){
			nhhist1->Fill(nhtpc);
		}
		if(pred==2){
			nhhist2->Fill(nhtpc);
		}
		if(pred==3){
			nhhist3->Fill(nhtpc);
		}
		T.ClearHistogram();
		hist->Reset("ICES");
		hist2->Reset("ICES");
		hist3->Reset("ICES");
	}
	c2->cd(1);
	nhhist1->Draw();
	c2->cd(2);
	nhhist2->Draw();
	c2->cd(3);
	nhhist3->Draw();
}

void ViewTPC(TString Filename){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	int entries = T.GetNEvent();
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);	
	c1->Divide(3,1);
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];

	for(int i=2996;i<entries;++i){
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
		} 
		T.AssignG4EventD(x,y,z,dedx);
		for(int j=0;j<T.GetNpadG4();++j){
			int padID= T.GetPadIDG4(j);
			TVector3 vec = T.GetG4Position(j);
			T.FillFlatHist(padID);
			T.FillHist(z[j],x[j]);
		}
		//		if(ThisEvent==L2PPi||ThisEvent==L2NPi){
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
		//		}
		if(i%1==0){
			T.ClearHistogram();
		}
	}
	cout<<"End"<<endl;
}

void TagTPCEvents(TString Filename,TString OutFileName){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	int entries = T.GetNEvent();
	TFile* outfile = new TFile(OutFileName,"recreate");
	//		TFile* outfile = new TFile("TaggedTrainDataKBeam.root","recreate");
	//		TFile* outfile = new TFile("TaggedTrainDataKPXi.root","recreate");
	//	TFile* outfile = new TFile("TaggedTrainDataP500.root","recreate");
	//	TFile* outfile = new TFile("TaggedTrainDataP700.root","recreate");
	TTree* outtree = new TTree("tree","tree");
	int evt_per_tag = 1000;
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];
	//	short x[max_nh];
	//	short y[max_nh];
	//	short z[max_nh];
	int TPCEventTag=0,nhtpc=0,ntrk=0;
	//	outtree->Branch("TPCEvent",TPCEvent,Form("TPCEvent[%d][%d][%d]/S",nbin,nbin,depth));
	int evnum;
	outtree->Branch("evnum",&evnum,"evnum/I");
	outtree->Branch("TPCEventTag",&TPCEventTag,"TPCEventTag/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	outtree->Branch("ntrk",&ntrk,"ntrk/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");

	TCanvas* c1 = new TCanvas("c1","c1",1500,700);	
	c1->Divide(2,1);
	int TagNumber[20] = {0};
	cout<<"Processing..."<<endl;
	for(int i=0;i<entries;++i){
		int ThisEvent = 0;
		//		evnum=(i+entries*i/3)%entries;
		evnum=i;
		T.SetEvent(evnum);//For mixing 3 files...
		nhtpc=T.GetNpadG4();
		if(i%1000==0)cout<<i<<endl;
		ThisEvent=T.WhichEvent();
		T.AssignG4EventD(x,y,z,dedx);
		TPCEventTag=ThisEvent;
		ntrk = T.NumberOfTracks(3); 
		if(TagNumber[ntrk]<evt_per_tag){
			outtree->Fill();
		}
		TagNumber[ntrk]++;
		if(ntrk>1000){

			T.FillEvent();
			cout<<"Ntrk = "<<ntrk<<endl;
			TString title = Form("evt=%d",evnum);
			T.SetTitle(title);
			c1->cd(1);
			T.DrawHist();
			c1->cd(2);
			T.DrawPosHist();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
		}
		if(i%1==0){
			T.ClearHistogram();
		}
	}
	outfile->Write();
}
void TestML(TString Filename="TrainData.root",TString TagName="PredictedData.root"){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	TFile* predfile = new TFile(tpcdir+TagName,"read");
	cout<<"FileOpen"<<endl;
	TTree* tree = (TTree*)predfile->Get("tree");
	cout<<"treeGet"<<endl;
	double Background,l2ppi,l2npi,kbeam;
	int evnum;
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("Background",&Background);
	tree->SetBranchAddress("Background",&Background);
	tree->SetBranchAddress("L2PPi",&l2ppi);
	tree->SetBranchAddress("L2NPi",&l2npi);
	tree->SetBranchAddress("KBeam",&kbeam);
	cout<<"BranchSet"<<endl;
	int scale = 1;
	int entries =int( T.GetNEvent()/scale) -2;
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);	
	c1->Divide(2,1);
	for(int i=2990;i<entries;++i){
		int ThisEvent = 0;
		tree->GetEntry(i);
		T.SetEvent(evnum);
		if(i%1000==0)cout<<i<<endl;
		ThisEvent=T.WhichEvent();
		for(int j=0;j<T.GetNpadG4();++j){
			int padID= T.GetPadIDG4(j);
			TVector3 vec = T.GetG4Position(j);
			double x = vec.X();double z = vec.Z();
			double x_t = tpc::getPosition(padID).X();
			double z_t = tpc::getPosition(padID).Z();
			T.FillFlatHist(padID);
			T.FillHist(z_t,x_t);
		}
		TString title = Form("evt=%d",evnum);
		T.SetTitle(title);
		cout<<"wait"<<endl;
		c1->cd(1);
		T.DrawHist();
		c1->cd(2);
		T.DrawPosHist();
		cout<<Form("Probability : (%.3f,%.3f,%.3f,%.3f)",Background,l2ppi,l2npi,kbeam)<<endl;
		cout<<Form("True: %d",ThisEvent)<<endl;
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
		if(i%1==0){
			T.ClearHistogram();
		}
	}
	cout<<"End"<<endl;
}




