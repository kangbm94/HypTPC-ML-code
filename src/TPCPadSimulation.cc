#include "TPCManager.cc"
static int nev = 0;
double Gaussianf(double* x, double* par){
	double dx = x[0];
	double peak = par[0];
	double mean = par[1];
	double sig = par[2];
	double val = (dx - mean)*(dx - mean)/sig/sig;//to avoid besselK0 divergence at 0
	double gaus = par[0]*exp(- 0.5*val );
	return gaus;
}
double Diffusion(double* x, double* par){
	double dx = x[0];
	double peak = par[0];
	double mean = par[1];
	double sig = par[2];
	double val = (dx - mean)*(dx - mean)/sig/sig;//to avoid besselK0 divergence at 0
	double gaus = par[0]*exp(- val / 2);
	double besselK = ROOT::Math::cyl_bessel_k(0,val/2);
	return gaus*besselK;
}
TF1* fDiffusion = new TF1("fDiffusion","Diffusion",-30,30,3);
TF1* fGaus = new TF1("fGaus","Gaussianf",-30,30,3);
void TPCPadSimulation(){
	cout<<"TPCLineSimulation()"<<endl;
	cout<<"TPCRKSimulation()"<<endl;
}
void TPCLineSimulation(){
	gTPCManager.InitializeHistograms();
	auto h = gTPCManager.GetPadHistogram();
	int nt = 100;
	double sig0 = 0.204;
	//*sqrt(30);
	int nevt = 500;//Number of Electrons
	TString dir = "../../MayRun/rootfiles/Defocus/";
	int runnum = 5754;
//	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);

//	gTPCManager.LoadTPCBcOut(dir+filename);
//	int ent = gTPCManager.GetEntries();

	//	cout<<"hi"<<endl;
//	gTPCManager.CreateFile(Form("Simul%d_5.root",runnum));
	gTPCManager.CreateFile("DiffusionSim2.root");
	gTPCManager.OutBranch("x0",0,1);
	gTPCManager.OutBranch("u0",1,1);
	gTPCManager.OutBranch("sig",2,1);
	gTPCManager.OutBranch("nel",3,1);
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
	TCanvas*c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(3,1);
	int maxn = 1000;
	double sigd = 0.57;
	TH2D* hist = new TH2D("hist","hist",100,-300,300,501,-10,10);
	for(int evt = 0 ; evt< maxn; evt ++){
		if((evt+1)%(maxn/100)==0) {cnt++;cout<<cnt<<"\% Processing"<<endl;}
		double x0 = 0;
		double u0 = (0.);
//		int nel = gRandom->Uniform(500,1500);
		int nel = 3000;
		for(int i=0;i<nel;++i){
			double z = -250+ (double)(500)/nel *(i+1);
			double y =  + 30;//mm -> cm, drift length 30cm +- y
			double sig = sigd*sqrt(y);
			double x = x0 +u0 * z; 
			auto pos = Generate2DGaus(z,x,sig);
//			TVector2 pos = TVector2(x+gRandom->Gaus(0,sig),z);
//			cout<<"Position: "<<pos.X()<<" , "<<pos.Y()<<endl;
			gTPCManager.FillHist(pos.X(),pos.Y());
			hist->Fill(pos.X(),pos.Y());
		}
		double val[10];val[0]=x0;val[1]=u0;val[2]=sigd;val[3]=nel;
		gTPCManager.AssignHits();
		gTPCManager.MakeUpClusters(-0.1);
#if 0
		h->Draw("colz");
		TF1*fLine = new TF1("fLine","[0]+[1]*(x+1318.9)",-250,250);
		fLine->SetParameter(0,x0); fLine->SetParameter(1,u0);
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
	gTPCManager.GetPadHistogram()->Draw("colz");
	c1->cd(2);
	hist->Draw("colz");
	c1->cd(3);
	auto H = hist->ProjectionY();
	H->Draw();
	fDiffusion->SetParLimits(2,0.3,10);
	H->Fit("fDiffusion");
	c1->Modified();
	fGaus->SetParLimits(2,0.3,10);
	H->Fit("fGaus");
	fGaus->SetLineColor(kBlue);
	cout<<"DiffusionWidth: "<<fDiffusion->GetParameter(2)/sqrt(30)<<endl;
	cout<<"GausDiffusionWidth: "<<fGaus->GetParameter(2)/sqrt(30)<<endl;
	//	chain->Draw("clsize");
}


void TPCManager::Process(double* vals){
	ClearBuffer();
	AssignHits();
	MakeUpClusters(0);
	double x0=vals[0],u0=vals[1],sig = vals[2],nel=vals[3];
	OutBufD[0]=x0;OutBufD[1]=u0;
	OutBufD[2]=x0;OutBufD[3]=nel;
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
			double xh = hit.GetPosition().X(); double zh = hit.GetPosition().Z();
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

