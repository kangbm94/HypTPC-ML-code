#include "TPCManager.cc"
static int nev = 0;
void TPCPadSimulation(){
	cout<<"TPCLineSimulation()"<<endl;
	cout<<"TPCRKSimulation()"<<endl;
}
void TPCLineSimulation(){
	gTPCManager.InitializeHistograms();
	auto h = gTPCManager.GetPadHistogram();
	int nt = 100;
	double sig0 = 0.204;
	double sigd = 0.57*sqrt(2);
	//*sqrt(30);
	int nevt = 500;//Number of Electrons
	TString dir = "../../MayRun/rootfiles/Defocus/";
	int runnum = 5754;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);

	gTPCManager.LoadTPCBcOut(dir+filename);
	int ent = gTPCManager.GetEntries();

	//	cout<<"hi"<<endl;
	gTPCManager.CreateFile(Form("Simul%d_5.root",runnum));
	gTPCManager.OutBranch("x0",0,1);
	gTPCManager.OutBranch("u0",1,1);
	gTPCManager.OutBranch("nev",0,0);
	gTPCManager.OutBranch("nhits",1,0);
	gTPCManager.OutBranch("ncls",2,0);

	gTPCManager.OutBranch("z",0,2);
	gTPCManager.OutBranch("x",1,2);
	gTPCManager.OutBranch("de",2,2);

	gTPCManager.OutBranch("PadID",3,2);
	gTPCManager.OutBranch("Bcx",4,2);
	gTPCManager.OutBranch("cluster",5,2);

	gTPCManager.OutBranch("clsize",6,2);
	gTPCManager.OutBranch("clz",7,2);
	gTPCManager.OutBranch("clx",8,2);
	gTPCManager.OutBranch("clde",9,2);
	gTPCManager.OutBranch("clBcx",10,2);
	gTPCManager.OutBranch("layer",11,2);
	int cnt = 0;
	ent/=5;
	TCanvas*c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(2,1);
	for(int evt = 4*ent ; evt< 5*ent; evt ++){
		gTPCManager.SetEvent(evt);
		int bcnt = gTPCManager.GetBCnt();
		int nhtpc = gTPCManager.GetNhits(0);
		if((evt+1)%(ent/100)==0) {cnt++;cout<<cnt<<"\% Processing"<<endl;}
		if(bcnt!=1) continue;
		if(nhtpc<9) continue;
		Track Trk = gTPCManager.GetTrack(0);
		double x0 = Trk.GetPosition(-K18HS).X();
		double u0 = (Trk.GetPosition(0).X()-x0)/K18HS;
		for(int i=0;i<nevt;++i){
			double z = -250+ (double)(500)/nevt *(i+1);
			double y = Trk.GetPosition(z).Y()/10 + 30;//mm -> cm, drift length 30cm +- y
			double sig = sigd*sqrt(y);
			double x = Trk.GetPosition(z).X(); 
			auto pos = Generate2DGaus(z,x,sig);
			//			cout<<"Position: "<<pos.X()<<" , "<<pos.Y()<<endl;
			gTPCManager.FillHist(pos.X(),pos.Y());
		}
		double val[10];val[0]=x0;val[1]=u0;
#if 0
		h->Draw("colz");
		TF1*fLine = new TF1("fLine","[0]+[1]*(x+1318.9)",-250,250);
		fLine->SetParameter(0,x0);
		fLine->SetParameter(1,u0);
		fLine->Draw("same");
		fLine->SetLineWidth(3);
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
#endif	
		gTPCManager.Process(val);
		gTPCManager.ClearHits();
		gTPCManager.ClearHistogram();
	}
	gTPCManager.WriteFile();
	auto* chain = gTPCManager.GetOutChain();
	c1->cd(1);
	chain->Draw("clsize");
	c1->cd(2);
	chain->Draw("clde");
}
void TPCRKSimulation(){
	gTPCManager.InitializeHistograms();
	auto h = gTPCManager.GetPadHistogram();
	int nt = 100;
	double sig0 = 0.204;
	double sigd = 0.18*sqrt(2);
	//*sqrt(30);
	int nevt = 500;//Number of Electrons
	TString dir = "../../MayRun/rootfiles/Defocus/";
	int runnum = 5754;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	gTPCManager.LoadTPCBcOut(dir+filename);
	int ent = gTPCManager.GetEntries();
	gTPCManager.CreateFile(Form("RKSimul%d_5.root",runnum));
	gTPCManager.OutBranch("x0",0,1);
	gTPCManager.OutBranch("u0",1,1);
	gTPCManager.OutBranch("nev",0,0);
	gTPCManager.OutBranch("nhits",1,0);
	gTPCManager.OutBranch("ncls",2,0);

	gTPCManager.OutBranch("z",0,2);
	gTPCManager.OutBranch("x",1,2);
	gTPCManager.OutBranch("de",2,2);

	gTPCManager.OutBranch("PadID",3,2);
	gTPCManager.OutBranch("Bcx",4,2);
	gTPCManager.OutBranch("cluster",5,2);

	gTPCManager.OutBranch("clsize",6,2);
	gTPCManager.OutBranch("clz",7,2);
	gTPCManager.OutBranch("clx",8,2);
	gTPCManager.OutBranch("clde",9,2);
	gTPCManager.OutBranch("clBcx",10,2);
	gTPCManager.OutBranch("layer",11,2);
	int cnt = 0;
	ent/=5;
	TCanvas*c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(2,1);
	for(int evt = 4*ent ; evt< 5*ent; evt ++){
		gTPCManager.SetEvent(evt);
		int bcnt = gTPCManager.GetBCnt();
		int nhtpc = gTPCManager.GetNhits(0);
		if((evt+1)%(ent/100)==0) {cnt++;cout<<cnt<<"\% Processing"<<endl;}
		if(bcnt!=1) continue;
		if(nhtpc<9) continue;
		Track Trk = gTPCManager.GetTrack(0);
		double x0 = Trk.GetPosition(-K18HS).X();
		double u0 = (Trk.GetPosition(0).X()-x0)/K18HS;
		for(int i=0;i<nevt;++i){
		
			double z = -250+ (double)(500)/nevt *(i+1);
			double y = Trk.GetPosition(z).Y()/10 + 30;//mm -> cm, drift length 30cm +- y
			double sig = sigd*sqrt(y);
			double x = Trk.GetPosition(z).X(); 
			auto pos = Generate2DGaus(z,x,sig);
			//			cout<<"Position: "<<pos.X()<<" , "<<pos.Y()<<endl;
			gTPCManager.FillHist(pos.X(),pos.Y());
		
		}
		double val[10];val[0]=x0;val[1]=u0;
		gTPCManager.Process(val);
		gTPCManager.ClearHits();
		gTPCManager.ClearHistogram();
	}



}


void TPCManager::Process(double* vals){
	ClearBuffer();
	AssignHits();
	MakeUpClusters(0);
	double x0=vals[0],u0=vals[1];
	OutBufD[0]=x0;OutBufD[1]=u0;
#if 0
#endif
	int nh = GetNumberOfMHits();
	int nclh = GetNumberOfMCls();
	OutBufI[0]=nev;
	OutBufI[1]=nh;
	OutBufI[2]=nclh;
	for(int i=0;i<nclh;++i){
		auto hcl = GetMCl(i);
		double x = hcl.GetPosition().X();
		double z = hcl.GetPosition().Z();
		double de = hcl.GetDe();
		double bcx = x0+u0*(z+K18HS);
		int clsize = hcl.GetSize();
		double cls = clsize;
		double layer = hcl.GetLayer();
		OutBufV[6]->push_back(cls);
		OutBufV[7]->push_back(z);
		OutBufV[8]->push_back(x);
		OutBufV[9]->push_back(de);
		OutBufV[10]->push_back(bcx);
		OutBufV[11]->push_back(layer);
		for(int j=0;j<cls;++j){
			auto hit = hcl.GetHit(j);
			double xh = hit.GetPosition().X();
			double zh = hit.GetPosition().Z();
			double deh = hit.GetDe();
			double bcxh = x0+u0*(zh+K18HS);
			int padID= tpc::findPadID(zh,xh);
			double pid=padID;
			int cluster = hit.GetCluster();
			double cl =cluster;
			OutBufV[0]->push_back(zh);
			OutBufV[1]->push_back(xh);
			OutBufV[2]->push_back(deh);
			OutBufV[3]->push_back(pid);
			OutBufV[4]->push_back(bcxh);
			OutBufV[5]->push_back(cl);

		}
	}


	OutChain->Fill();
	nev++;
}
