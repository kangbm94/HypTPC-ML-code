#include "TPCPadHelper.hh"
using namespace tpc;
void test(){
	for(int i=1;i<5500;++i){
		TVector3 pos = getPosition(i);
		double x = pos.X(),z=pos.Z();
		int padID = findPadID(z,x);
		if(i!=padID){
			cout<<Form("Warning! padID!= i.  i,padID = (%d,%d)",i,padID)<<endl;
		}
	}
}
