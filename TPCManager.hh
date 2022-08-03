#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
TF1* f_gaus = new TF1("f_gaus","gaus");
TF1* f_landau = new TF1("f_landau","landau");
TF1* slew_func= new TF1("slew_func",SlewFunc,0,5,3);
const int max_nh=2500;
enum{
	Else=0,
	L2PPi=1,
	L2NPi=2,
	KBeam=3
};
const int nbin=250;const int depth=1;
const double tpc_size=250;
short ToPixel(double x){
	x+=250;
	short x_pix = int(x* (double)nbin/tpc_size/2);
	return x_pix;
}
int ToShort(double y){
	y+=350;
	y*=10;
	return short(y);
}

static const int nhtpcmax = 300;
class TPCManager:public FileManager{
	private:
		TFile* hist_file;
		TH2Poly* PadHist;
		TH2I* FlatHist;
		TH2D* PosHist;
		TH3D* SpaceHist;
		TGraph2D* SpaceGraph;
		TGraph2D* SpaceGraphBase;
		vector<int> *padTpc;
		int iPadtpc[nhtpcmax];
		double xtpc[nhtpcmax];
		double ytpc[nhtpcmax];
		double ztpc[nhtpcmax];
		double dedxtpc[nhtpcmax];
		int idtpc[nhtpcmax];
		int ititpc[nhtpcmax];
		int ntrk[nhtpcmax];
		vector<double>* dlTpc;
		vector<double>* deTpc;
		int htofnhits;
		int htofhitpat[34];
		int nhittpc; double htofua[34];
		int gp = 0;
		int gpb = 0;
	public:
		TPCManager(){};
		virtual void LoadFile(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadChain("tpc");
		}
		virtual void LoadG4File(TString FileName){
			DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadG4Chain("TPC_g");
		}
		void LoadChain(TString ChainName);
		void LoadG4Chain(TString ChainName);
		int GetNEvent(){
			return DataChain->GetEntries();
		};
		void SetEvent(int evt){
			DataChain->GetEntry(evt);
		};
		void ClearHistogram(){
			PadHist->Reset("ICE");
			PosHist->Reset("ICE");
			FlatHist->Reset("ICE");
			gp=0;
		}
		int GetNpad(){
			return Min(padTpc->size(),max_nh);
		};
		int GetPadID(int i){
			return padTpc->at(i);
		};
		int GetNpadG4(){
			return Min(nhittpc,max_nh);
		}
		int GetPadIDG4(int i){
			return iPadtpc[i];
		}
		int Getidtpc(int i){
			return idtpc[i];
		}
		int Getititpc(int i){
			return ititpc[i];
		}
		int Getntrk(int i){
			return ntrk[i];
		}
		int GetTrackNum(){
			return ntrk[nhittpc-1];
		}
		double GetDL(int i){
			return dlTpc->at(i);
		}
		double GetDE(int i){
			return deTpc->at(i);
		}
		double Getdedxtpc(int i){
			return dedxtpc[i];
		}
		TVector3 GetPosition(int padID){
			TVector3 pos =  tpc::getPosition(padID);
			return pos;
		}
		TVector3 GetG4Position(int i){
			return TVector3(xtpc[i],ytpc[i],ztpc[i]);
		}
		int GetHTOFMT(){
			return htofnhits;
		}
		void GetHTOFHitPat(int* hp){
			for(int i=0;i<34;++i){
				hp[i]=htofhitpat[i];
			}
		}
		void SetTitle(TString title){
			PadHist->SetTitle(title);
			FlatHist->SetTitle(title);
		}
		void InitializeHistograms();
		void FillHist(double x, double z){
			PadHist->Fill(x,z);
			PosHist->Fill(x,z);
		};
		void FillHist(int padID);
		void FillFlatHist(int padID);
		void Fill3DHist(TVector3 vect){
			SpaceHist->Fill(vect.Z(),vect.X(),vect.Y());	
		}
		void FillBaseGraph(TVector3 vect){
			SpaceGraphBase->SetPoint(gpb,vect.Z(),vect.X(),vect.Y());	
			gpb++;
		}
		void Fill3DGraph(TVector3 vect){
			SpaceGraph->SetPoint(gp,vect.Z(),vect.X(),vect.Y());	
			gp++;
		}
		void MakeTPCPad();
		void DrawHist(){
			PadHist->Draw("colz");
		}
		void DrawFlatHist(){
			FlatHist->Draw("colz");
		}
		void DrawPosHist(){
			PosHist->Draw("colz");
		}
		void Draw3DHist(){
			SpaceHist->Draw("colz");
		}
		void Draw3DGraph(){
			SpaceGraph->Draw("APSAME");
		}
		TVector3 GetRTheta(int padID);
		TVector2 GetLayerRow(int padID);
		int WhichEvent();
		void AssignG4Event(short * x,short* y,short* z,double* dedx);
		void AssignG4EventD(double * x,double* y,double* z,double* dedx);
		void AssignRealEvent(double * x,double* y,double* z,double* dedx);
		void FillEvent();
		int NumberOfTracks(int min_points=6);
};
void TPCManager::LoadChain(TString ChainName ){
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("padTpc",&padTpc);
	DataChain->SetBranchAddress("dlTpc",&dlTpc);
	DataChain->SetBranchAddress("deTpc",&deTpc);
	DataChain->SetBranchAddress("htofnhits",&htofnhits);
	DataChain->SetBranchAddress("htofhitpat",htofhitpat);
	//	DataChain->SetBranchAddress("htofua",htofua);
};
void TPCManager::LoadG4Chain(TString ChainName ){
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("nhittpc",&nhittpc);
	DataChain->SetBranchAddress("iPadtpc",iPadtpc);
	DataChain->SetBranchAddress("xtpc",xtpc);
	DataChain->SetBranchAddress("ytpc",ytpc);
	DataChain->SetBranchAddress("ztpc",ztpc);
	DataChain->SetBranchAddress("idtpc",idtpc);
	DataChain->SetBranchAddress("ititpc",ititpc);
	DataChain->SetBranchAddress("ntrk",ntrk);
	DataChain->SetBranchAddress("dedxtpc",dedxtpc);
}
void TPCManager::InitializeHistograms(){
	PadHist = tpc::InitializeHistograms();
	FlatHist = new TH2I("PadRTheta","PadRTheta",32,0,32,240,0,240);
	PosHist = new TH2D("PosHisto","PosHisto",128,-250,250,128,-250,250);
	/*
		 SpaceHist = new TH3D("TPCTrack","TPCTrack",260,-260,260,260,-260,260,300,-300,300);
		 SpaceGraph = new TGraph2D();
		 SpaceGraph->GetXaxis()->SetRangeUser(-300,300);
		 SpaceGraphBase->GetXaxis()->SetTitle("Z");
		 SpaceGraphBase->GetYaxis()->SetRangeUser(-300,300);
		 SpaceGraphBase->GetYaxis()->SetTitle("Y");
		 SpaceGraphBase->GetZaxis()->SetRangeUser(-300,300);
		 SpaceGraphBase->GetXaxis()->SetTitle("X");
		 SpaceGraphBase->SetMarkerColor(kRed);
		 SpaceGraphBase->SetMarkerSize(3);
		 MakeTPCPad();
		 SpaceGraphBase->SetMarkerStyle(kCircle);
		 SpaceGraphBase->SetMarkerColor(kBlue);
		 SpaceGraphBase->SetMarkerSize(3);
		 */
}
TVector3 TPCManager::GetRTheta(int padID){
	TVector3 pos = GetPosition(padID);
	return pos;
}
TVector2 TPCManager::GetLayerRow(int padID){
	int layer = tpc::getLayerID(padID);
	int row = tpc::getRowID(padID);
	TVector2 idvec(layer,row);
	return idvec;
}
void TPCManager::FillFlatHist(int padID){
	TVector2 lr = GetLayerRow(padID);
	int l = lr.X();
	int r = lr.Y();
	FlatHist->Fill(l,r);
}
void TPCManager::FillHist(int padID){
	TVector3 hitv = GetPosition(padID);
	double x = hitv.X();
	double z = hitv.Z();
	PadHist->Fill(z,x);
};
void TPCManager::MakeTPCPad(){
	for(int i=0;i<max_padid;++i){
		TVector3 vect = GetPosition(i+1);
		vect.SetY(-300);
		FillBaseGraph(vect);
	}
	SpaceGraph->SetMarkerStyle(kCircle);
	SpaceGraph->SetMarkerColor(kBlue);
	SpaceGraph->SetMarkerSize(3);
}
int TPCManager::WhichEvent(){
	const int max_ntrk=20;
	int particle[max_ntrk]={0};
	int ThisEvent=0,npi=0,nk=0,np=0;
	int nh = GetNpadG4();
	for(int j=0;j<nh;++j){
		particle[Getntrk(j)]=Getidtpc(j);
	}
	for(int j=0;j<max_ntrk;++j){
//		cout<<((particle[j]))<<endl;
		switch(abs(particle[j])){
			case PionID:
				npi++;
				break;
			case KaonID:
				nk++;
				break;
			case ProtonID:
				np++;
				break;
			case 3312:
				//				cout<<"Xi detected"<<endl;
				break;
			case 0:
				break;
			default:
				//				cout<<Form("Warning! pid: %d",particle[j])<<endl;
				break;
		}
	}
	if(npi==2&&nk==2&&np==1){
		ThisEvent=L2PPi;
	}
	else if(npi==1&&nk==2&&np==0){
		ThisEvent=L2NPi;
	}
	else if(npi==0&&nk==1&&np==0){
		ThisEvent=KBeam;
	}
	else{
		ThisEvent=Else;
	}
	//	cout<<ThisEvent<<endl;
	return ThisEvent;
}
int TPCManager::NumberOfTracks(int min_points=6){
	const int max_ntrk=20;
	int counter[max_ntrk]={0};
	int nh = GetNpadG4();
	int TrackCount=0;
	for(int j=0;j<nh;++j){
		counter[Getntrk(j)]+=1;
	}
	for(int j=0;j<max_ntrk;++j){
		if(counter[j]>(min_points-1)) TrackCount++;
	}
	return TrackCount;
}
void TPCManager::AssignG4Event( short* x,short* y,short* z,double* dedx){
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;
	}
//	cout<<"Initialized"<<endl;
	for(int j=0;j<GetNpadG4();++j){
		TVector3 vec = GetG4Position(j);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		short x_pix=ToPixel(x_);short y_short = ToShort(y_);short z_pix=ToPixel(z_);
		x[j]=x_pix;y[j]=y_short;z[j]=z_pix;dedx[j]=Getdedxtpc(j);
	}
}
void TPCManager::AssignG4EventD( double* x,double* y,double* z,double* dedx){
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;
	}
//	cout<<"Initialized"<<endl;
	for(int j=0;j<GetNpadG4();++j){
		TVector3 vec = GetG4Position(j);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		x[j]=x_;y[j]=y_;z[j]=z_;dedx[j]=Getdedxtpc(j);
	}
}
void TPCManager::AssignRealEvent( double* x,double* y,double* z,double* dedx){
	int nitr = GetNpad();
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;dedx[j]=0;
	}
	for(int j=0;j<nitr;++j){
		TVector3 vec = GetPosition(GetPadID(j));
		double dl= GetDL(j);
		vec.SetY(dl);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		x[j]=x_;y[j]=y_;z[j]=z_;dedx[j]=GetDE(j);
	}
}
void TPCManager::FillEvent(){
	int nitr = GetNpadG4();
	double x[max_nh],y[max_nh],z[max_nh],dedx[max_nh];
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;
	}
	AssignG4EventD(x,y,z,dedx);
	for(int j=0;j<nitr;++j){
		FillHist(z[j],x[j]);
	}
}




