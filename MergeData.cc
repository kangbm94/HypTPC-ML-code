const int max_nh=2500;
void MergeData(){
	cout<<"Merge(TString out, TStirng file1, TString file2)"<<endl;
}
void Merge(TString out, TString file1, TString file2){
	TFile* f1 = new TFile(file1,"read");
	TFile* f2 = new TFile(file2,"read");
	TTree* tree1 = (TTree*)f1->Get("tree");
	TTree* tree2 = (TTree*)f2->Get("tree");
	cout<<"FileLoaded"<<endl;
	int nhtpc1,nhtpc2,evnum1,evnum2,tag1,tag2,ev1,ev2,ntrk1,ntrk2;
	double x1[max_nh],y1[max_nh],z1[max_nh],x2[max_nh],y2[max_nh],z2[max_nh];
	double x[max_nh],y[max_nh],z[max_nh];
	tree1->SetBranchAddress("nhtpc",&nhtpc1);
	tree1->SetBranchAddress("ntrk",&ntrk1);
	tree1->SetBranchAddress("x",x1);
	tree1->SetBranchAddress("y",y1);
	tree1->SetBranchAddress("z",z1);
	tree2->SetBranchAddress("nhtpc",&nhtpc2);
	tree2->SetBranchAddress("ntrk",&ntrk2);
	tree2->SetBranchAddress("x",x2);
	tree2->SetBranchAddress("y",y2);
	tree2->SetBranchAddress("z",z2);
	TFile* outfile = new TFile("taggedReal.root","recreate");
	TTree* outtree = new TTree("tree","tree");
	int nhtpc,ntrk,TPCEventTag=0,evnum;
	double dedx[max_nh]={0};
	outtree->Branch("evnum",&evnum,"evnum/I");
	outtree->Branch("TPCEventTag",&TPCEventTag,"TPCEventTag/I");
	outtree->Branch("nhtpc",&nhtpc,"nhtpc/I");
	outtree->Branch("ntrk",&ntrk,"ntrk/I");
	outtree->Branch("x",x,"x[nhtpc]/D");
	outtree->Branch("y",y,"y[nhtpc]/D");
	outtree->Branch("z",z,"z[nhtpc]/D");
	outtree->Branch("dedx",dedx,"dedx[nhtpc]/D");
	for(int i=0;i<10000;++i){
		if(i<1000){
			tree1->GetEntry(i);	
			nhtpc=nhtpc1;
			ntrk=ntrk1;
			for(int j=0;j<nhtpc1;++j){
				x[j]=x1[j];
				y[j]=y1[j];
				z[j]=z1[j];
			}
			outtree->Fill();
		}
		tree2->GetEntry(i);
		nhtpc=nhtpc2;
		ntrk=ntrk2;
		for(int j=0;j<nhtpc2;++j){
			x[j]=x2[j];
			y[j]=y2[j];
			z[j]=z2[j];
		}
		outtree->Fill();
	}
	outfile->Write();
}
