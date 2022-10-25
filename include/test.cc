#include "Geometry.hh"
void test(){
	RiemanSphere RS(2,0,0,3);
	Plane Pl;
	for(int i=0;i<30;++i){
		Pl.GeneratePlanarCircle(4,3,2,0.0,0.5*Pi());
//		Pl.GeneratePlanarCircle(3000,0,0,0,2*Pi());
	}
	cout<<"AddingRandomPoints"<<endl;
	for(int i=0;i<10;++i){
		double x_ = -1+2*gRandom->Rndm();		
		double y_ = -1+2*gRandom->Rndm();		
		TPCPoint pt(x_,y_,0);
		Pl.AddPoint(pt);
	}
	int np = Pl.GetNPoints();
	cout<<"ProjectingPoints"<<endl;
	for(int i=0;i<np;++i){
		RS.ProjectPoint(Pl.GetPoint(i));
		
//		Pl.GetPoint(i).ListElements();
//		RS.GetPoint(i).ListElements();
	}
	cout<<"CalculatingLSPlane"<<endl;
	auto SmallCircle = RS.LeastSquarePlane();
	SmallCircle.Show();
	cout<<"LSPlane Determined"<<endl;
	int nps = SmallCircle.GetNPoints();
//	cout<<"CalculatingD2s: "<<nps<<endl;
	for(int i=0;i<nps;++i){
		cout<<Form("%d th D2 = %f",i,SmallCircle.D2(SmallCircle.GetPoint(i)))<<endl;
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	auto* h = RS.GetRiemanHist();
//	RS.DrawRiemanHist();
//		h->Draw("col");
	h->SetMarkerStyle(kFullDiamond);
	h->Draw("");
	cout<<h->GetEffectiveEntries()<<endl;
		//	h->Draw("");
}
