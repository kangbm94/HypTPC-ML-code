#include "../include/TPCTrajectory.hh"

bool DscendingY(const TPCPoint& lhs, const TPCPoint& rhs){
	return lhs.Y()>rhs.Y();
}

void TrackCandidate::SortY(){
	sort(points.begin(),points.end(),DscendingY);
}
void TrackCandidate::SortZ(){
}
void TrackCandidate::MakeHistogram(int nh){
	histZX		= new TH2D(Form("histZX_%d",nh),"zx",128/2,-250,250,128/2,-250,250);
	histZY		= new TH2D(Form("histZY_%d",nh),"zy",128/2,-250,250,128/2,-300,350);
	histXY		= new TH2D(Form("histXY_%d",nh),"xy",128/2,-250,250,128/2,-300,350);
	histX		= new TH1D(Form("histX_%d",nh),"x",128/2,-250,250);
	histY		= new TH1D(Form("histY_%d",nh),"y",128/2,-300,350);
	histZ		= new TH1D(Form("histZ_%d",nh),"z",128/2,-250,250);	
	int np = NumberOfPoints();
	for(int i=0;i<np;++i){
		TPCPoint a = GetPoint(i);
		histZX->Fill(a.Z(),a.X());histZY->Fill(a.Z(),a.Y());histXY->Fill(a.X(),a.Y());	
		histX->Fill(a.X());histY->Fill(a.Y());histZ->Fill(a.Z());	
	}
}

TH2D* TrackCandidate::Get2DHist(int conf){
	switch(conf){
		case 1:
			return histZY;
		case 2:
			return histXY;
	deafult:
			return histZX;
	}
	return histZX;
}
TH1D* TrackCandidate::Get1DHist(int conf){
	switch(conf){
		case 1:
			return histY;
		case 2:
			return histZ;
	deafult:
			return histX;
	}
	return histX;
}
double TPCTrack::ClosestDistance(TPCPoint a){
	double cd = 1e9;
	int np =NumberOfPoints();
	for(int i=0;i<np;++i){
		cd=Min(cd,GetPoint(i).ClosestDistance(a));
	}
	return cd;
}

void TrackCandidate::ListElements(){
	int np = points.size();
	for(int i=0;i<np;++i){
		points[i].ListElements();
	}
}

void TPCPoint::ListElements(){
		cout<<Form("(%f,%f,%f)",X(),Y(),Z())<<endl;
}

bool TPCTrack::FitBeamCircle(double mom, double* par){
	double rad = 3000*mom/0.9;
	TString minimizer = "Minuit2",Algorithm="Migrad";
//	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(TString minimizer,TString Algorithm);
	
//	ROOT::Math::Functor fcirc(CircleMetric,3);
	ROOT::Fit::Fitter fitter;


	return true;
}

