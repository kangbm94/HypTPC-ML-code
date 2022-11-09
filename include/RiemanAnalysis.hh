#include "Geometry.hh"
#include "TPCManager.hh"
class RiemanManager: public TPCManager{
	private:
		RiemanSphere Rieman(-143,0,0,100);
	public:
		RiemanManager(){}
		void ProjectPoints();
};

void RiemanManager::ProjectPoints(){
	int nh = GetNHits(0);
	for(int i=0;i<nh;++i){
		TVector3 v = GetPosition(i);
		double vx = v.x(),vy=v.y(),vz=v.z();
		TVector3 vp(z,x,y);
		RS.ProjectPoint(vp);
	}
	auto smc = RS.LeastSquarePlane();
	smc.Show();
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	auto* h = RS.GetRiemanHist();
//	RS.DrawRiemanHist();
//		h->Draw("col");
	h->SetMarkerStyle(kFullDiamond);
	h->Draw("");
	cout<<h->GetEffectiveEntries()<<endl;
}



