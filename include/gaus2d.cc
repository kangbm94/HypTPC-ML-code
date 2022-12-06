#include "Math.hh"
void gaus2d(){

	double sig = 1;
	double x0=0,z0=0;
	TH2D* hist = new TH2D("h2","h2",200,-2.5,2.5,200,-2.5,2.5);
	TH1D* hist2 = new TH1D("h3","h3",200,-3,3);
	for(int i=0;i<1e5;++i){
		double rand = gRandom->Uniform(-5,5);
		auto v = Generate2DGaus(z0,x0+rand,sig);
		double x = v.X(),y=v.Y();
		hist->Fill(x,y);
		hist2->Fill(x);
	}
	auto h2 = hist->ProjectionY("ht",90,110);
	//h2->Draw();
	//	hist2->Draw();
//	hist2->Draw("colz");
	hist2->Fit("gaus");
}
