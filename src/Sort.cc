#include "../include/Utils.hh"
void Sort(){
	TString filename = base_dir + "TPCInv05643.root";
	TFile* file = new TFile(filename);
	TTree* tree = (TTree*)file->Get("tree");
	TFile* Sorted = new TFile("Sorted.root","recreate");
	TTree* SortedTree = new TTree("tree","tree");
	int runnum,evnum;
	bool FlgLd,FlgXi;
	double MM,InvMXi;
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("FlgLd",&FlgLd);
	tree->SetBranchAddress("FlgXi",&FlgXi);
	tree->SetBranchAddress("MM",&MM);
	tree->SetBranchAddress("InvMXi",&InvMXi);
	SortedTree->Branch("runnum",&runnum);
	SortedTree->Branch("evnum",&evnum);
	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		bool Acpt = true;
		//		bool Acpt = abs(MM-1.315)<0.1 and abs(InvMXi-1.315)<0.1;
		if(FlgLd and FlgXi and Acpt) SortedTree->Fill();
//		SortedTree->Fill();
	}
	cout<<SortedTree->GetEntries()<<endl;
	Sorted->Write();
}
