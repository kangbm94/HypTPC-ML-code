#include "../include/KinematicLV.hh"
#ifndef KinematicLV_cc
#define KinematicLV_cc
KinematicLV::KinematicLV(double M_, vector<TLorentzVector> D_){
	Mass = M_;
	TLorentzVector V(0,0,0,0);
	for(auto v : D_){
		D.push_back(v);
		V += v;
	}
	this = V;
}
#endif
