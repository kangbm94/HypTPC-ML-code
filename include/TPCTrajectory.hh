#include "TPCPoint.hh"
#ifndef TPCTrajectory_h
#define TPCTrajectory_h
static int nh=0;
class TrackCandidate{
	protected:
		vector<TPCPoint> points;
		TH2D* histZX=NULL;
		TH2D* histZY=NULL;
		TH2D* histXY=NULL;
		TH1D* histX=NULL;
		TH1D* histY=NULL;
		TH1D* histZ=NULL;
	public:
		TrackCandidate(){
			nh++;
		}
		TrackCandidate(int np,double* x,double* y,double* z,double* de){
			for(int i=0;i<np;++i){
				TPCPoint a(x[i],y[i],z[i],de[i]);
				points.push_back(a);
			}
			nh++;
		}
		~TrackCandidate(){
			if(histZX!=NULL){
				delete histZX;delete histZY;delete histXY;delete histX;delete histY;delete histZ;
			}
		}
		void AssignPoint(TPCPoint a){
			points.push_back(a);
		}
		void RemovePoint(int a){
			points.erase(points.begin()+a);
		}
		TPCPoint GetPoint(int i){
			return points.at(i);
		}
		void SortY();
		void SortZ();
		TH2D* Get2DHist(int conf);//0(or default)->ZX, 1-> ZY, 2->XY
		TH1D* Get1DHist(int conf);//0(or default)->X, 1-> Y, 2->Z
		void ListElements();
		int NumberOfPoints(){
			return points.size();
		}
		void MakeHistogram(int n);
		bool CircleHough(double mom,double* pars);
};

class TPCTrack:public TrackCandidate{
	private:
		int TrackNo=0;
		int ParticleId=0;
		bool IsBeam = false;
		//Helix Parameters//
		double cx,cy,r,cz,dz;
		//Line Parameters//
		double x0,y0,z0,u,v;
		TGraphErrors* gZX=NULL;
		TGraphErrors* gZY=NULL;
		TGraphErrors* gXY=NULL;
	public:
		TPCTrack(){
			nh++;
		}
		TPCTrack(int np,double* x,double* y,double* z,double* de)
			:TrackCandidate(np,x,y,z,de){
			}
		~TPCTrack():~TrackCandidate(){
			if(gZX!=NULL){
				delete gZX;delete gZY;delete gXY;
			}
		}
		double ClosestDistance(TPCPoint a);
		bool FitBeamCircle(double mom, double* par);
		function<double(double*)> CircleMetric=[this](double* par){
			double cx=par[0],cz=par[1],r=par[3];
			int np = points.size();
			double val=0;
			for(int i=0;i<np;++i){
				val+=CircleD2(points[i],cx,cz,r);
				val+=1;
			}
			return val;
		};
};
#endif
