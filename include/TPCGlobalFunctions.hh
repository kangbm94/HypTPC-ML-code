#include "Math.hh"
#ifndef TPCGlobalFunctions_h
#define TPCGlobalFunctions_h
bool IsInsideTPC(double x, double y, double z){
	if(sqrt(x*x+z*z)<250) return true;
	else return false;
}
#endif
