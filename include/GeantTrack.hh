#ifndef GeantTrack.hh
#define GeantTrack.hh
class GeantTrack{
	private:
		int TrackID,PID;
		//TPC
		vector<TVector3> TPCHitPos;
		vector<TVector3> TPCHitMom;
		vector<double>	TPCHitDeDx;
		vector<double>	TPCHitEdep;
		
		//FToF
		double length;
		double tof;
		TVector3 Mom;
	public:
		GeantTrack(int TrackID){}
		void PushTPCHit(int id, TVector3 pos, TVector3 mom)
	
};
#endif
