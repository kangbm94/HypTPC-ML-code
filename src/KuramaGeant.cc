#include "GeantManager.cc"
TString file = "../SimData/KPXI.root";
void KuramaGeant(){
	TH2D* h = new TH2D("scat","scat",100,-1,2,100,0,2.5);
	gGeantManager.LoadFile(file);
	int nev = gGeantManager.GetEntries();
	for(int iev = 0;iev < nev; ++iev){
		gGeantManager.SetEvent(iev);
		int nhtof = gGeantManager.GetNhToF();
		for(int it = 0; it < nhtof; ++it){
//			double q = gGeantManager.GetQToF(it);
			double q = 1;
			h->Fill(q*sqrt(gGeantManager.GetM2(it)),gGeantManager.GetPToF(it));
		}
	}
	h->Draw("colz");
}
