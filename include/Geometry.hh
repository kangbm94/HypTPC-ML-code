#include "Utils.hh"
#include "TPCPoint.hh"
#include "Matrix.hh"
#ifndef Geometry_h
#define Geometry_h
class Surface{
	protected:
		vector<TPCPoint> Points;
	public:
		Surface(){};
		TPCPoint GetPoint(int i);
		void AddPoint(TPCPoint pt){Points.push_back(pt);}
		double Position(int inp,int xyz){
//			cout<<"Pos "<<inp<<endl;
			return Points.at(inp)[xyz];};
};
TPCPoint Surface::GetPoint(int i){
	int np = Points.size();
	if(i>np){
		cout<<Form("np = %d! but request is %d!",np,i)<<endl;
		return Points[0];}
	else return Points[i];
}
class Plane:public Surface{
	private:
		TVector3 PlaneVector;//Equation of plane: ax + by + cz - 1 = 0.

	public:
		Plane(){};
		Plane(double* dir){
			PlaneVector=TVector3(dir);
		};
		Plane(TVector3 dir){
			PlaneVector=dir;
		};
		void GeneratePlanarCircle(double rad,double cx, double cy, double st, double len);
		int GetNPoints(){
			return Points.size();
		}
		void Show(){cout<<Form("Vector(%f,%f,%f)",PlaneVector[0],PlaneVector[1],PlaneVector[2])<<endl;}
		double D1(TPCPoint point){return PlaneVector*point-1;}
		double D2(TPCPoint point){return D1(point)*D1(point);};
		bool IsInside(TPCPoint point, double threshold);
		void Rotate(double theta);
};
void Plane::GeneratePlanarCircle(double rad,double cx,double cy,double st,double len){
	TVector2 pos = GenCircleRandom(rad,st,st+len);
	double x = (pos.X()+cx);
	double y = (pos.Y()+cy);
	double z= 0;
	TPCPoint pt(x,y,z);
	Points.push_back(pt);
}

void Plane::Rotate(double theta){
	int np = Points.size();
	vector<TPCPoint> cont;
	double angle = theta/180*Pi();
	for(int i=0;i<np;++i){
		double x =Points[i].X();	
		double y =Points[i].X();
		double z = Points[i].Z();
		double de = Points[i].GetEdep();
		double x_= x*cos(angle)+y*sin(angle);
		double y_= -x*sin(angle)+y*cos(angle);
		TPCPoint pt(x_,y_,z,de);
		cont.push_back(pt);
	}
	//	TPCPoints.clear();
	Points = cont;
}
bool Plane::IsInside(TPCPoint point,double threshold=1){
	if(D2(point)<threshold) return true;
	else return false;
}
#endif

class RiemanSphere:public Surface{
	private:
		TPCPoint Center;
		TPCPoint NorthPole;
		double radius;
		TH3D* hist;
	public:
		RiemanSphere(double x,double y,double z, double r){
			radius = r;
			Center.SetPosition(x,y,z);
			NorthPole.SetPosition(x,y,z+r);
			hist = new TH3D("RiemanSphereH","RiemanSphereH",22,x-1.1*r,x+1.1*r,22,y-1.1*r,y+1.1*r,22,z-1.1*r,z+1.1*r);
			//			hist = new TH2D("RiemanSphere","RiemanSphere",22,x-1.1*r,x+1.1*r,22,y-1.1*r,y+1.1*r);
		}
		RiemanSphere(){radius=1;}
		void ProjectPoint(TPCPoint point);
		TH3D* GetRiemanHist(){
			return hist;
		}
		void DrawRiemanHist(){
			hist->Draw();
		}
		Plane LeastSquarePlane();
};
void RiemanSphere::ProjectPoint(TPCPoint point){
	double dist = Euclidean(NorthPole,point);
	double de = point.GetEdep();
	auto direction = point-NorthPole;
	double t = 2*radius*(point.Z()-NorthPole.Z())/dist/dist;
	auto ProjectedPosition = NorthPole-t*direction;
	TPCPoint ProjectedPoint(ProjectedPosition,de);
	Points.push_back(ProjectedPoint);
	double px=ProjectedPosition.X() ;double py=ProjectedPosition.Y() ;double pz=ProjectedPosition.Z();
	hist->Fill(px,py,pz);

	double chk = Euclidean(Center,ProjectedPoint);
}
Plane RiemanSphere::LeastSquarePlane(){
	 double MatEl[5][5]={0};
	 double Constants[3]={0};
	int np = Points.size();
	cout<<"Number of Points : "<<np<<endl;
	 double planevector[3];
	for(int inp=0;inp<np;++inp){
		for(int ixyz = 0;ixyz<3;++ixyz){
			for(int jxyz = 0;jxyz<3;++jxyz){
				MatEl[ixyz][jxyz]+=(Position(inp,ixyz))*(Position(inp,jxyz));
			}
			Constants[ixyz]+=Position(inp,ixyz);
		}
	}
	Matrix Mat(3,3,MatEl);
	Mat.SetConstants(Constants);
	Mat.GaussianElimination(planevector);
	double pv[3];
	for(int i=0;i<3;++i){
		pv[i]=planevector[i];
	}
	auto Pl= Plane(pv);
	for(int i=0;i<np;++i){
		Pl.AddPoint(Points.at(i));
	}
	return Pl;
}
