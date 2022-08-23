#include "src/TPCEvent.cc"
#include "src/TPCTrajectory.cc"
#include "src/TPCManager.cc"
void test(){
}
void TestDistance(){
	TPCPoint a(23.,22.,11.,1.5);
	auto b = a;
	cout<<b.X()<<endl;
	cout<<a.GetEdep()<<endl;
	a.SetDirection(0,0,1);
	cout<<a.GetDirection().X()<<endl;
	TPCPoint perp(3./5,4./5,0,2);
	perp*=7;
	cout<<"perp_edep : "<<perp.GetEdep()<<endl;
	cout<<"perp_mag : "<<perp.Mag()<<endl;
	TPCPoint c=a;
	c+=(10.738*(a.GetDirection())); c+=perp;
	cout<<"Dist : "<<Distance(a,c)<<endl;
	cout<<"CDist : "<<a.ClosestDistance(c)<<endl;
	cout<<"c_edep : "<<c.GetEdep()<<endl;
}
void TestSort(){
	TPCManager TM;
	TString dir = "../MayRun/rootfiles/CH2/TPC";
	TString tpcdir = dir+"/";
	int runnum = 5641;
	TString Filename = "TPCCluster0"+to_string(runnum)+".root";
	TString TPCFile = tpcdir+Filename;
	TM.LoadClusterFile(TPCFile);
	cout<<"FileLoaded"<<endl;
	double x[max_nh],y[max_nh],z[max_nh],de[max_nh]; 
	TCanvas* c2 = new TCanvas("c2","c2",1600,700);	
	c2->Divide(4,4);
	TCanvas* c1 = new TCanvas("c1","c1",1600,700);	
	c1->Divide(2,1);
	for(int i=0;i<10;++i){
		TM.SetEvent(i);
		cout<<"EventSet"<<endl;
		TPCTrack TC(TM.AssignRealEvent(x,y,z,de),x,y,z,de);
		vector<TPCTrack> Tracks;
		c1->cd(1);
		TH2D* hzx = TC.Get2DHist(2);
		hzx->Draw();
		c1->cd(2);
		TH1D* hy = TC.Get1DHist(1);
		hy->Draw();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
	}
}
void TestPresort(){
	TPCManager TM;
	TString dir = "../MayRun/rootfiles/CH2/TPC";
	TString tpcdir = dir+"/";
	int runnum = 5641;
	TString Filename = "TPCCluster0"+to_string(runnum)+".root";
	TString TPCFile = tpcdir+Filename;
	TM.LoadClusterFile(TPCFile);
	cout<<"FileLoaded"<<endl;
	double x[max_nh],y[max_nh],z[max_nh],de[max_nh]; 
	TCanvas* c2 = new TCanvas("c2","c2",1600,700);	
	c2->Divide(4,4);
		TM.SetEvent(2);
		cout<<"EventSet"<<endl;
		TPCTrack TC(TM.AssignRealEvent(x,y,z,de),x,y,z,de);
	PresortedEvent PE;
	PE.Presort(TC,30);
	cout<<"NumberOfTracks: "<<PE.NumberOfTracks()<<endl;

	for(int i=0;i<PE.NumberOfTracks()&&i<16;++i){
		TPCTrack TR = PE.GetTrack(i);
		TR.MakeHistogram(i);
		c2->cd(i+1);
		TR.Get2DHist(2)->Draw("colz");
	}
}
