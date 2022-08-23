#include "TPCTrajectory.hh"
#ifndef TPCEvent_h
#define TPCEvent_h
class TPCEvent{
	protected:
		vector<TPCTrack> Tracks;
	public:
		void AssignTrack(TPCTrack a){
			Tracks.push_back(a);
		}
		TPCTrack GetTrack(int i){
			return Tracks.at(i);
		}
		int NumberOfTracks(){
			return Tracks.size();
		}
};
class PresortedEvent: public TPCEvent{
	protected:

	public:
		void Presort(TrackCandidate TC,double cd_cut);
};
#endif
