#include "TPCManager.cc"
#include "../include/Math.hh"
double TransverseDistance(double xc,double zc, double x, double z){
	double dist = sqrt((xc-x)*(xc-x)+(zc-z)*(zc-z));
	if(xc<x) dist=-dist;
	return dist;
};



double Gaussianf(double* x, double* par){
	double dx = x[0];
	double peak = par[0];
	double mean = par[1];
	double sig = par[2];
	double val = (dx - mean)*(dx - mean)/sig/sig;//to avoid besselK0 divergence at 0
	double gaus = par[0]*exp(- val );
	return gaus;
}
double Diffusion(double* x, double* par){
	double dx = x[0];
	double peak = par[0];
	double mean = par[1];
	double sig = par[2];
	double val = (dx - mean)*(dx - mean)/sig/sig;//to avoid besselK0 divergence at 0
	double gaus = par[0]*exp(- val / 2);
	double besselK = ROOT::Math::cyl_bessel_k(0,val/2);
	return gaus*besselK;
}
TF1* fDiffusion = new TF1("fDiffusion","Diffusion",-30,30,3);
TF1* fGaus = new TF1("fGaus","Gaussianf",-30,30,3);

void LineSimAnal(){

	TFile* file = new TFile("DiffusionSim2.root");
//	TFile* file = new TFile("Simul.root");
	TTree* tree = (TTree*)file->Get("tree");

	vector<double>* z = new vector<double>;
	vector<double>* x = new vector<double>;
	vector<double>* de = new vector<double>;
	vector<double>* cluster = new vector<double>;
	vector<double>* clsize = new vector<double>;
	vector<double>* clz = new vector<double>;
	vector<double>* clx = new vector<double>;
	vector<double>* clde = new vector<double>;
	vector<double>* layer = new vector<double>;

	tree->SetBranchAddress("z",&z);
	tree->SetBranchAddress("x",&x);
	tree->SetBranchAddress("de",&de);
	tree->SetBranchAddress("cluster",&cluster);
	tree->SetBranchAddress("clsize",&clsize);
	tree->SetBranchAddress("clz",&clz);
	tree->SetBranchAddress("clx",&clx);
	tree->SetBranchAddress("clde",&clde);
	tree->SetBranchAddress("layer",&layer);


	TCanvas* c1 = new TCanvas("c1","c1",1400,700);
	c1->Divide(8,4);
//	TCanvas* c2 = new TCanvas("c2","c2",1400,700);
//	c2->Divide(8,4);
	TCanvas* c3 = new TCanvas("c3","c3",1000,600);
	c3->Divide(2,1);
	TH2D* hist[32];
	for(int i=0;i<32;++i){
		TString title = Form("layer%d",i);
		hist[i] = new TH2D(title,title,300,-15,15,100,0,1);
	}
	double A;
	int ent = tree->GetEntries();
	TGraph* g1 = new TGraph();
	TGraph* g2 = new TGraph();
	for(int i=0;i<ent;++i){
		if(!(i%1000)) cout<<i<< " / "<<ent<<endl;
		tree->GetEntry(i);
		int ncl = clde->size();
		int nh = de->size();
		for(int j=0;j<ncl;++j){
			double ecl = 0;
			double cl_x = clx->at(j);
			double cl_z = clz->at(j);
			double cl_de = clde->at(j);
			int cl_size = 0;
			for(int k = 0; k<nh;++k){
				if(j==cluster->at(k)){
					ecl+=de->at(k);	
					cl_size++;
				}	
			}
			for(int k = 0; k<nh;++k){
				if(j==cluster->at(k)){
					if(ecl!=cl_de){
						cout<<"Warning! wrong cluster,clsize vs cl_size = "<< clsize->at(j)<<" , "<<cl_size<<" clde vs desum = "<<cl_de<<" , "<<ecl<<endl;
					}
					A = de->at(k)/ecl;	
					int lay = (int)layer->at(j);
					double h_x = x->at(k);
					double h_z = z->at(k);
					double dist = TransverseDistance(cl_x,cl_z,h_x,h_z);
					hist[lay]->Fill(dist,A);
				}	
			}
		}
	}
	gStyle->SetOptFit(1100);
	fGaus->SetParLimits(2,0.1,10);
	fDiffusion->SetLineColor(kBlue);
	fDiffusion->SetParLimits(0,0.,100.);
	fDiffusion->SetParLimits(2,0.1,10);
	for(int i=0;i<32;++i){
		c1->cd(i+1);
		gPad->SetLogz();
		cout<<"Layer "<<i<<endl;
		hist[i]->Draw("colz");
		if(hist[i]->GetEffectiveEntries()<100) continue;
		hist[i]->Fit("fGaus","Q0");
		hist[i]->Fit("fDiffusion","Q0");
		fGaus->Draw("same");
		fDiffusion->Draw("same");
		c1->Modified();
		gPad->SetLogz();
		double sig = fGaus->GetParameter(2)/sqrt(30);
		double sig2 = fDiffusion->GetParameter(2)/sqrt(30);
		g1->AddPoint(i+1,sig);
		g2->AddPoint(i+1,sig2);
		cout<<"Gaus: "<<sig<<endl;
		cout<<"Bessel: "<< sig2<<endl;
	}
	c3->cd(1);
	g1->Draw("AP");
	g1->SetMarkerStyle(140);
	g1->SetTitle("Gaussian");
	g1->GetXaxis()->SetTitle("Layer");
	g1->GetYaxis()->SetTitle("DiffusionWidth(Fit) [mm]");
	c3->cd(2);
	g2->Draw("AP");
	g2->SetMarkerStyle(140);
	g2->SetTitle("GausXBesselK0");
	g2->GetXaxis()->SetTitle("Layer");
	g2->GetYaxis()->SetTitle("DiffusionWidth(Fit) [mm]");
}
