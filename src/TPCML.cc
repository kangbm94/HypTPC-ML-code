#include "TPCManager.cc"
#ifndef TPCML_C
#define TPCML_C
TPCManager T;
void ConvertRealTPC(int runnum){
	TString dir = "../MayRun/rootfiles/CH2/TPC";
	TString tpcdir = dir+"/";
	//	TString Filename = "TPCHit0"+to_string(runnum)+".root";
	TString Filename = "TPCCluster0"+to_string(runnum)+".root";
	TString TPCFile = tpcdir+Filename;
	T.LoadClusterFile(TPCFile);
	int ent =T.GetNEvent();
	TString OutFileName = "RealClData0"+to_string(runnum)+".root";
	TFile* outfile = new TFile(OutFileName,"recreate");
	T.WriteTag("Data",TPCFile);
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	int htofnhits;
	int htofhitpat[34];
	double dedx[max_nh];
	int nhtpc=0,evnum=0;
	TTree* outtree = new TTree("tree","tree");
	outtree->Branch("evnum",&evnum,"evnum/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	//	outtree->Branch("htofnhits",&htofnhits,"htofnhits/I");
	//	outtree->Branch("htofhitpat",htofhitpat,"htofhitpat[34]/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");

	for(int i=0;i<ent;++i){
		if(i%1000==0)cout<<Form("Processing %d th event...",i)<<endl;
		T.SetEvent(i);	
		evnum=i;
		nhtpc =Min(T.GetNhits(), max_nh);
		//		htofnhits=T.GetHTOFMT();
		//		T.GetHTOFHitPat(htofhitpat);
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
	TH1I* nhhist1 = new TH1I("h1","nh=1",100,0,300);
	TH1I* nhhist2 = new TH1I("h2","nh=2",100,0,300);
	TH1I* nhhist3 = new TH1I("h3","nh=3",100,0,300);
	TH1I* nhhist4 = new TH1I("h4","nh=4",100,0,300);
	TH1I* nhhist5 = new TH1I("h5","nh=5",100,0,300);
	TH1I* nhhist6 = new TH1I("h6","nh=6",100,0,300);
	int ent = predtree->GetEntries();
	int trks=2;
	bool force=false;
	int nh_cutm[6]={0,60,80,90,125,150};
	int nh_cutM[6]={50,75,90,120,155,200};
	for(int i=000000;i<ent;++i){
		int evt = i;
		if(i%10000==0){
			cout<<i<<endl;
		}
		predtree->GetEntry(evt);
		datatree->GetEntry(evnum);
		ntrk=pred;
		for(int j=0;j<nhtpc;++j){
			//			T.FillHist(z[j],x[j]);
			hist->Fill(z[j],x[j]);
			hist2->Fill(z[j],y[j]);
			hist3->Fill(x[j],y[j]);
		}
		//		if(pred<trks&&pred!=0&&nh_cutm[pred-1]-1<nhtpc&&nhtpc<nh_cutM[pred-1]+1){
		//	}
		if(ntrk!=-1){
			cout<<"Number of Hits: "<<nhtpc<<endl;
			if(ntrk==1){
				cout<<"Xi"<<endl;
			}
			else{
				cout<<"Other"<<endl;
			}
			//			cout<<"Predicted Number of Tracks: "<<pred<<endl;
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
		if(pred==4){
			nhhist4->Fill(nhtpc);
		}
		if(pred==5){
			nhhist5->Fill(nhtpc);
		}
		if(pred==6){
			nhhist6->Fill(nhtpc);
		}
		T.ClearHistogram();
		hist->Reset("ICES");
		hist2->Reset("ICES");
		hist3->Reset("ICES");
	}
}
void TagRealEvents(TString Filename,TString PredFile){
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

	TFile* outfile = new TFile("train_real.root","recreate");
	TTree* outtree = new TTree("tree","tree");
	double dedx[max_nh]={0};
	int TPCEventTago=0,nhtpco=0,ntrko=0,evnumo;
	outtree->Branch("evnum",&evnumo,"evnum/I");
	outtree->Branch("TPCEventTag",&TPCEventTago,"TPCEventTag/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	outtree->Branch("ntrk",&ntrk,"ntrk/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");


	TH2D* hist = new TH2D("hist","zx",128,-250,250,128,-300,350);
	TH2D* hist2 = new TH2D("hist2","zy",128,-250,250,128,-300,350);
	TH2D* hist3 = new TH2D("hist3","xy",128,-250,250,128,-300,350);
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(3,1);
	TCanvas* c2 = new TCanvas("c2","c2",1500,800);
	c2->Divide(3,2);
	TH1I* nhhist1 = new TH1I("h1","nh=1",100,0,300);
	TH1I* nhhist2 = new TH1I("h2","nh=2",100,0,300);
	TH1I* nhhist3 = new TH1I("h3","nh=3",100,0,300);
	TH1I* nhhist4 = new TH1I("h4","nh=4",100,0,300);
	TH1I* nhhist5 = new TH1I("h5","nh=5",100,0,300);
	TH1I* nhhist6 = new TH1I("h6","nh=6",100,0,300);
	int ent = predtree->GetEntries();
	int trks=6;
	bool force=false;
	int nh_cutm[6]={0,60,95,120,160,195};
	int nh_cutM[6]={50,85,120,150,185,220};
	for(int i=000000;i<ent;++i){
		int evt = i;
		if(i%10000==0){
			cout<<i<<endl;
		}
		predtree->GetEntry(evt);
		datatree->GetEntry(evnum);
		if(pred<trks&&pred!=0&&nh_cutm[pred-1]-1<nhtpc&&nhtpc<nh_cutM[pred-1]+1){
			ntrk=pred;
			outtree->Fill();
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
		if(pred==4){
			nhhist4->Fill(nhtpc);
		}
		if(pred==5){
			nhhist5->Fill(nhtpc);
		}
		if(pred==6){
			nhhist6->Fill(nhtpc);
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
	c2->cd(4);
	nhhist4->Draw();
	c2->cd(5);
	nhhist5->Draw();
	c2->cd(6);
	nhhist6->Draw();
	outfile->Write();
}

void ViewTPC(TString Filename){
	TString dir = ".";
	TString tpcdir = dir+"/";
	TString TPCFile = tpcdir+Filename;
	T.LoadG4File(TPCFile);
	int entries = T.GetNEvent();
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];
	int trkid[max_nh],pid[max_nh];
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
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
		} 
		T.AssignG4EventD(trkid,pid,x,y,z,dedx);
		for(int j=0;j<T.GetNhitsG4();++j){
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
	int evt_per_tag = 10000;
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];
	int trkid[max_nh];
	int pid[max_nh];
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
	outtree->Branch("trkid",trkid,"trkid[nhtpc]/I");
	outtree->Branch("pid",pid,"pid[nhtpc]/I");
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
		nhtpc=T.GetNhitsG4();
		if(i%1000==0)cout<<i<<endl;
		ThisEvent=T.WhichEvent();
		T.AssignG4EventD(trkid,pid,x,y,z,dedx);
		TPCEventTag=ThisEvent;
		ntrk = T.NumberOfTracks(3); 
		if(TagNumber[ntrk]<evt_per_tag&&ntrk==0){
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
		for(int j=0;j<T.GetNhitsG4();++j){
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



void TagRealTPCByHand(int runnum , int mode,TString comment = "no comment"){
	//TString Filename = "RealClData0"+to_string(runnum)+".root";
	TString Filename = "RealData0"+to_string(runnum)+".root";
	TFile* datafile = new TFile(Filename,"read");
	TTree* datatree = (TTree*)datafile->Get("tree");
	TString Xidir = "../MayRun/RunSummary/CH2/";
	TString XiFile = Xidir+Form("Xi_Production0%d.root",runnum);
	TFile* file = new TFile(XiFile,"READ");
	TTree* tree = (TTree*)file->Get("tree");
	int XiEv;
	double XiM2;
	tree->SetBranchAddress("XiEv",&XiEv);
	tree->SetBranchAddress("XiM2",&XiM2);

	int evnum;
	int nhtpc;
	int TPCEventTag;
	int ntrk;
	double x[max_nh];
	double y[max_nh];
	double z[max_nh];
	double dedx[max_nh];
	int htofnhits;
	datatree->SetBranchAddress("nhtpc",&nhtpc);
	datatree->SetBranchAddress("x",x);
	datatree->SetBranchAddress("y",y);
	datatree->SetBranchAddress("z",z);
	datatree->SetBranchAddress("htofnhits",&htofnhits);


	TString outfilename; 
	int start = 0;
	if(mode==1){
		outfilename = "XiTaggedReal"+to_string(runnum)+".root";
	}
	else{
		outfilename = "Background"+to_string(runnum)+".root";
	}
	TFile* outfile = new TFile(outfilename,"recreate");
	T.WriteTag("Data",Filename);
	T.WriteTag("Conf","XiTag");
	T.WriteTag("misic",comment);
	TTree* outtree = new TTree("tree","tree");
	outtree->Branch("evnum",&evnum,"evnum/I");
	outtree->Branch("TPCEventTag",&TPCEventTag,"TPCEventTag/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	outtree->Branch("ntrk",&ntrk,"ntrk/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");

	TH2D* hist = new TH2D("hist","zx",128,-250,250,128,-300,350);
	TH2D* hist2 = new TH2D("hist2","zy",128,-250,250,128,-300,350);
	TH2D* hist3 = new TH2D("hist3","xy",128,-250,250,128,-300,350);
	TCanvas* c1 = new TCanvas("c1","c1",1000,800);
	c1->Divide(2,2);
	TH1D* nhhist = new TH1D("h1","M2",100,1,2);
	//	int ent = datatree->GetEntries();
	int ent = tree->GetEntries();
	bool force=false;
	int cnt=0;
	double xim = 1.317;
	double xiw = 0.012;
	double xistarm = 1.539;
	double xistarw = 0.010;
	double nsig = 3.;
	if(mode==1){
		for(int i=start;i<ent;++i){
			tree->GetEntry(i);
			int evt = XiEv; 
			datatree->GetEntry(evt);
			evnum=XiEv;
			bool Xi=false,XiStar=false;
			if(nhtpc>10){
				if( sqrt((XiM2-xim)*(XiM2-xim))<nsig*xiw ){
					Xi=true;
					cout<<"Xi, M2 = "<<XiM2<<endl;
					TPCEventTag = 1;
				}
				if( sqrt((XiM2-xistarm)*(XiM2-xistarm))<nsig*xistarw ){
					XiStar=true;
					cout<<"XiStar, M2 = "<<XiM2<<endl;
					TPCEventTag = 2;
				}
				if(!Xi&&!XiStar){
					cout<<"BG, M2 = "<<XiM2<<endl;
					TPCEventTag = 0;
					if(cnt>80){
						continue;
					}
					cout<<cnt<<endl;
					cnt++;
				}
			}
			else{
					TPCEventTag = 0;
					if(cnt>80){
						continue;
					}
					cout<<cnt<<endl;
					cnt++;
			}
			if(cnt%100==0&&cnt!=0){
				cout<<i<<endl;
				outfile->Write();
			}
			for(int j=0;j<nhtpc;++j){
				hist->Fill(z[j],x[j]);
				hist2->Fill(z[j],y[j]);
				hist3->Fill(x[j],y[j]);
			}
			c1->cd(1);
			hist2->Draw("colz");
			c1->cd(2);
			hist3->Draw("colz");
			c1->cd(3);
			hist->Draw("colz");
			c1->Modified();
			c1->Update();
//			cout<<"Number of Hits: "<<nhtpc<<endl;
//			cout<<"htofnhits : "<<htofnhits<<endl;
			gSystem->ProcessEvents();
			//			cout<<"What's the type of this event? : (0->Else, 1->Xi-like, 2->Xi*-like)"<<endl;
//			cin.ignore();
			if(ntrk<0){
				cout<<"quitting..."<<endl;
				break;
			}
			if(ntrk<11){
				outtree->Fill();
//				cnt++;
			}
			else{
				hist->Reset("ICES");
				hist2->Reset("ICES");
				hist3->Reset("ICES");
				continue;
			}
			nhhist->Fill(XiM2);
			c1->cd(4);
			nhhist->SetTitle(Form("Evt%d",i));
			nhhist->Draw();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			hist->Reset("ICES");
			hist2->Reset("ICES");
			hist3->Reset("ICES");
		}
		outfile->Write();
	}
	else{
		for(int i=start;i<100;++i){
			tree->GetEntry(i);
			int evt = XiEv; 
			datatree->GetEntry(1000*i);
			for(int j=0;j<nhtpc;++j){
				hist->Fill(z[j],x[j]);
				hist2->Fill(z[j],y[j]);
				hist3->Fill(x[j],y[j]);
			}
			c1->cd(1);
			hist2->Draw("colz");
			c1->cd(2);
			hist3->Draw("colz");
			c1->cd(3);
			hist->Draw("colz");
			c1->Modified();
			c1->Update();
			cout<<"Number of Hits: "<<nhtpc<<endl;
			cout<<"htofnhits : "<<htofnhits<<endl;
			gSystem->ProcessEvents();

			TPCEventTag=0;
			outtree->Fill();
			cnt++;
			nhhist->Fill(ntrk);
			c1->cd(4);
			nhhist->SetTitle(Form("Evt%d",i));
			nhhist->Draw();
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			hist->Reset("ICES");
			hist2->Reset("ICES");
			hist3->Reset("ICES");
		}
		outfile->Write();
	}
}
void TagRealTPC(){
	for(int i=5641;i<5667;++i){
		TagRealTPCByHand(i,1);
	}
}
#endif
