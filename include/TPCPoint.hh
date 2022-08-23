#include "Math.hh"
#ifndef TPCPoint_h
#define TPCPoint_h
class TPCPoint:public TVector3{
	private:
		double sx,sy,sz;//sigma;
		double dx,dy,dz;//Direction
		double edep;
		bool InTrack = false;
	public:
		TPCPoint(double x_,double y_,double z_,double edep_=0){
			this->SetX(x_);this->SetY(y_);this->SetZ(z_);edep=edep_;
		}
		TPCPoint();
		void SetDirection(double dx_,double dy_, double dz_){
			double norm_=	Norm(dx_,dy_,dz_);
			dx=dx_/norm_;dy=dy_/norm_;dz=dz_/norm_;
		}
		TVector3 GetDirection(){
			return TVector3(dx,dy,dz);
		}
		double GetEdep(){
			return edep;
		}
		double ClosestDistance(TPCPoint a);
		void SetTrackingStatus(bool tf){
			InTrack = tf;
		}
		void ListElements();
};
double TPCPoint::ClosestDistance(TPCPoint a){
	auto SubVec = *this-a;
	auto DirVec = this->GetDirection();
	double V = SubVec.Mag();
	double L = SubVec*DirVec;
	double D = sqrt(V*V-L*L);
	return D;
}


double Euclidean(TPCPoint a, TPCPoint b){
	return (a-b).Mag();
}

double Distance(TPCPoint a, TPCPoint b, int type=0){//0 for Euclidean metric;
	if(type==0) return Euclidean(a,b);
	else return 0;
}
double CircleD2(TPCPoint a,double cx,double cz,double r){
	double dx = cx-a.X();
	double dz = cz-a.Z();
	double dist= r-sqrt(dx*dx-dz*dz);
	return dist*dist;
}


#endif
