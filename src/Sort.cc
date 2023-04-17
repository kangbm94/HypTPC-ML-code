#include "../include/Utils.hh"
void Sort(){
	TString filename = base_dir + "TPCInv_cd30.root";
	TFile* file = new TFile(filename);
	TTree* tree = (TTree*)file->Get("tree");
	TFile* Sorted = new TFile("Sorted.root","recreate");
	TTree* SortedTree = new TTree("tree","tree");
	int runnum,evnum;
	int ntK18,ntKurama;
	bool FlgLd,FlgXi;
	double MM,InvMXi,mmpi0;
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("FlgLd",&FlgLd);
	tree->SetBranchAddress("FlgXi",&FlgXi);
	tree->SetBranchAddress("ntK18",&ntK18);
	tree->SetBranchAddress("ntKurama",&ntKurama);
	tree->SetBranchAddress("MM",&MM);
	tree->SetBranchAddress("mmpi0",&mmpi0);
	tree->SetBranchAddress("InvMXi",&InvMXi);
	SortedTree->Branch("runnum",&runnum);
	SortedTree->Branch("evnum",&evnum);
	SortedTree->Branch("MM",&MM);
	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
//		bool Acpt = true;
				bool Acpt = abs(MM-1.535)<0.1 and abs(InvMXi-1.315)<0.1 and mmpi0 > 0;
//				bool Acpt = abs(MM-1.321)<0.1 and abs(InvMXi-1.321)<0.1;
		if(FlgLd and FlgXi and Acpt and runnum == 5644 ){
			//and ntK18==1){
			cout<<Form("Run, Ev = (%d,%d)",runnum,evnum)<<endl;
			SortedTree->Fill();
		}
//		SortedTree->Fill();
	}
	cout<<SortedTree->GetEntries()<<endl;
	Sorted->Write();
}
