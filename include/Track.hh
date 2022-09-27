#include "Utils.hh"
#ifndef Track_hh
#define Track_hh
class Track{
	private:
		double x0,y0,u0,v0;
		double chisqr;
	public: 
		Track(){}
		Track(double x0_, double y0_,double u0_,double v0_){
			x0=x0_,y0=y0_,u0=u0_,v0=v0_;
		}
		TVector2 GetPosition(double z){
			z+=K18HS;
			return TVector2(x0+z*u0,y0+z*v0);
		}
};
#endif
