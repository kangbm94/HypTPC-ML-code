void view(){
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(4,2);
//	TFile* file = new TFile("ValidationData.root","read");
	TFile* file = new TFile("TrainedData.root","read");
	TTree* tree = (TTree*)file->Get("tree");
	TH2D* hists[8];
	cout<<"Histo"<<endl;
	for(int i=0;i<8;++i){
		hists[i]=new TH2D(Form("h%d",i),Form("h%d",i),8,0,8,20,0,1);
	}
	int tag,pred;
	tree->SetBranchAddress("tag",&tag);
	tree->SetBranchAddress("pred",&pred);
	cout<<"Fill"<<endl;
	const int ncl=7;
	double ntag[ncl]={0};
	double acc[ncl]={0};
	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		ntag[tag]++;
		if(tag==pred){
			acc[tag]++;
		}
	}
	for(int i=0;i<ncl;++i){
		cout<<Form("%d Tracks : Accuracy %.2f",i,acc[i]/ntag[i])<<endl;
	}
	tree->Draw("pred:tag","","colz");
}
