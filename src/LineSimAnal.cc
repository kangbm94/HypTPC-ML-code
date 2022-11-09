#include "TPCManager.cc"
void TPCManager::Process(double* vals){
};
double TransverseDistance(double xc,double zc, double x, double z){
	double dist = sqrt((xc-x)*(xc-x)+(zc-z)*(zc-z));
	if(xc<x) dist=-dist;
	return dist;
};
void LineSimAnal(){

	TFile* file = new TFile("Simul05754.root","READ");
	TTree* tree = (TTree*)file->Get("tree");
	vector<double>* z = new vector<double>;
	vector<double>* x = new vector<double>;
	vector<double>* de = new vector<double>;
	vector<double>* cluster = new vector<double>;
	vector<double>* clz = new vector<double>;
	vector<double>* clx = new vector<double>;
	vector<double>* clde = new vector<double>;
	vector<double>* layer = new vector<double>;

	tree->SetBranchAddress("z",&z);
	tree->SetBranchAddress("x",&x);
	tree->SetBranchAddress("de",&de);
	tree->SetBranchAddress("cluster",&cluster);
	tree->SetBranchAddress("clz",&clz);
	tree->SetBranchAddress("clx",&clx);
	tree->SetBranchAddress("clde",&clde);
	tree->SetBranchAddress("layer",&layer);


	TCanvas* c1 = new TCanvas("c1","c1",1400,700);
	c1->Divide(8,4);
	TH2D* hist[32];
	for(int i=0;i<32;++i){
		TString title = Form("layer%d",i);
		hist[i] = new TH2D(title,title,300,-15,15,100,0,1);
	}
	double A;
	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		int ncl = clde->size();
		int nh = de->size();
		for(int j=0;j<ncl;++j){
			double ecl = 0;
			double cl_x = clx->at(j);
			double cl_z = clz->at(j);
			double cl_de = clde->at(j);
			for(int k = 0; k<nh;++k){
				if(j==cluster->at(k)){
					ecl+=de->at(k);	
				}	
			}
			for(int k = 0; k<nh;++k){
				cout<<j<<" -> Cluster : "<<cluster->at(k)<<endl;
				if(ecl!=cl_de) cout<<"Warning! wrong cluster, clde vs desum = "<<cl_de<<" , "<<ecl<<endl;
				if(j==cluster->at(k)){
					A = de->at(k)/ecl;	
					int lay = (int)layer->at(k);
					double h_x = x->at(k);
					double h_z = z->at(k);
					double dist = TransverseDistance(cl_x,cl_z,h_x,h_z);
					hist[lay]->Fill(dist,A);
				}	
			}
		}
	}
	gStyle->SetOptFit(1100);
	for(int i=0;i<32;++i){
		c1->cd(i+1);
		hist[i]->Draw("colz");
		hist[i]->Fit("gaus");
	}
}
