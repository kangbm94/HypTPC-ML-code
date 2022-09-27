#include "Utils.h"
#ifndef RiemanSphere_h
#define RiemanSphere_h
TVector3 Project2Rieman(double r,double phi){
	double x=r*cos(phi)/(1+r*r);
	double y=r*sin(phi)/(1+r*r);
	double z=r*r/(1+r*r);
	return TVector3(x,y,z);
}
#endif
