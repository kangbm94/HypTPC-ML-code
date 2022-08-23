#include "../include/TPCEvent.hh"
void PresortedEvent::Presort(TrackCandidate TC,double cd_cut=10){
	TC.SortY();
	TPCTrack a;
	int cnt=0;
	a.AssignPoint(TC.GetPoint(0));
	AssignTrack(a);
	cout<<"1st Track: "<<endl;
	a.ListElements();
	int np = TC.NumberOfPoints();
	cout<<"NP: "<<np<<endl;
	for(int i=1;i<np;++i){
		TPCPoint point = TC.GetPoint(i);
		int nt = Tracks.size();
		bool cd_flag=false;
		for(int j=0;j<nt;++j){
			cnt++;
			cout<<i<<" th point, LOOP :"<<j<<endl;
			double cd = Tracks.at(j).ClosestDistance(point);
			if(cd<cd_cut){
				(Tracks[j]).AssignPoint(point);
				cd_flag=true;
			}
//			Tracks.at(j).ListElements();
		}
		if(!cd_flag){
			cout<<"New Track!" <<endl;
			TPCTrack Track_new;
			Track_new.AssignPoint(point);
			AssignTrack(Track_new);
		}
	}
}

