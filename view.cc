TString ht[8]={"0track","1track","2tracks","3tracks","4tracks","5tracks","6tracks","7orMoreTracks"};
void view(){
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(4,2);
	TFile* file = new TFile("ValidationData.root","read");
//	TFile* file = new TFile("TrainedData.root","read");
//	TFile* file = new TFile("PredictedData.root","read");
	//	TFile* file = new TFile("PredictedData700.root","read");
	TTree* tree = (TTree*)file->Get("tree");
	const int ncl=2;
	TH2D* hists[ncl];
	cout<<"Histo"<<endl;
	for(int i=0;i<ncl;++i){
		hists[i]=new TH2D(Form("h%d",i),Form("h%d",i),8,0,8,20,0,1);
	}
	int tag,pred;
	tree->SetBranchAddress("tag",&tag);
	tree->SetBranchAddress("pred",&pred);
	cout<<"Fill"<<endl;
	double ntag[ncl]={0};
	double acc[ncl]={0};
	int ent = tree->GetEntries();
	TH1I* n_hists[ncl];
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		ntag[tag]++;
		if(tag==pred){
			acc[tag]++;
		}
	}
	for(int i=0;i<ncl;++i){
		n_hists[i]=new TH1I(ht[i],ht[i],8,0,8);
		c1->cd(i+1);
		cout<<Form("%d Tracks : Accuracy %.2f",i,acc[i]/ntag[i])<<endl;
		tree->Draw("pred>>"+ht[i],Form("tag==%d",i));
	}
	TCanvas* c2 = new TCanvas("c2","c2",1500,800);
	tree->Draw("pred:tag","","colz");
}
