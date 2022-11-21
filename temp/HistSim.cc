#include "HistSim.hh"
void HistSim(){
	TH2Poly* hist = tpc::InitializeHistogram();
	TestClass TC(hist);
	TC.FillHist(48);
	TC.Draw();
}
