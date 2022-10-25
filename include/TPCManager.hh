#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
#include "TPCGlobalFunctions.hh"
#include "Track.hh"
#ifndef TPCManager_h
#define TPCManager_h
const int max_nh=2500;
double Target_pos=-143,Target_x=34,Target_z=24;
double frame_width = 12.5;
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
bool IsInside(double z,double x){
	if(sqrt(z*z+x*x)<250){
		return true;
	}
	else{
		return false;
	}
}
bool IsTarget(double z, double x){
	if(abs(z-Target_pos)<Target_z/2&&abs(x)<Target_x/2){
		return true;
	}
	else{
		return false;
	}
}
bool IsActiveArea(double z, double x){
	double zr = z-Target_pos;
	if(abs(abs(z)-abs(x))<12.){
		return false;
	}
	else{
		return !IsTarget(z,x) and IsInside(z,x);
	}
}
static const int nhtpcmax = 300;
class TPCManager:public FileManager{
	protected:
		TFile* hist_file;
		TH2Poly* PadHist=nullptr;
		TH2I* FlatHist=nullptr;
		TH2D* PosHist=nullptr;
		vector<int> *padTpc;
		int iPadtpc[nhtpcmax];
		double xtpc[nhtpcmax];
		double ytpc[nhtpcmax];
		double ztpc[nhtpcmax];
		double dedxtpc[nhtpcmax];
		int idtpc[nhtpcmax];
		int ititpc[nhtpcmax];
		int ntrk[nhtpcmax];
		vector<double>* dlTpc = new vector<double>;
		vector<double>* deTpc = new vector<double>;
		vector<double>* clxTpc = new vector<double>;
		vector<double>* clyTpc = new vector<double>;
		vector<double>* clzTpc = new vector<double>;
		vector<int>* clsize = new vector<int>;
		int htofnhits;
		int htofhitpat[34];
		int nhittpc; double htofua[34];
		int gp = 0;
		int gpb = 0;
		bool cluster = false;
		int ntBcOut=0;
		vector<double>* x0BcOut = new vector<double>;
		vector<double>* y0BcOut = new vector<double>;
		vector<double>* u0BcOut = new vector<double>;
		vector<double>* v0BcOut = new vector<double>;
	public:
		TPCManager(){};

		//Data I/O
		virtual void LoadFile(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadChain("tpc");
		}
		virtual void LoadClusterFile(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadClusterChain("tpc");
		}
		virtual void LoadG4File(TString FileName){
			DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadG4Chain("TPC_g");
		}
		void LoadChain(TString ChainName);
		void LoadClusterChain(TString ChainName);
		void LoadG4Chain(TString ChainName);
		void LoadBcOut();
		int GetNEvent(){
			return DataChain->GetEntries();
		};
		void SetEvent(int evt){
			DataChain->GetEntry(evt);
		};


		//Histogram Methods//
		void InitializeHistograms();
		void SetTitle(TString title){
			PadHist->SetTitle(title);
			FlatHist->SetTitle(title);
		}
		void FillHist(double x, double z){
			PadHist->Fill(x,z);
			PosHist->Fill(x,z);
		};
		void FillHist(int itr);
		void FillFlatHist(int padID);
		void SetPadContent(int padID,double cont);
		void SetPadContent(int layer,int row,double cont);
		void DrawHist(){
			PadHist->Draw("colz");
		}
		void DrawFlatHist(){
			FlatHist->Draw("colz");
		}
		void DrawPosHist(){
			PosHist->Draw("colz");
		}
		void ClearHistogram(){
			PadHist->Reset("ICE");
			PosHist->Reset("ICE");
			FlatHist->Reset("ICE");
			gp=0;
		}
		TH2Poly* GetPadHistogram(){
			return PadHist;
		}



		int GetNhits(int conf){
			if(!conf)	return Min(padTpc->size(),max_nh);
			else 					return Min(clsize->size(),max_nh);
		};
		int GetPadID(int i){
			return padTpc->at(i);
		};
		int GetNhitsG4(){
			return Min(nhittpc,max_nh);
		}
		int GetPadIDG4(int i){
			return iPadtpc[i];
		}
		int GetClSize(int i){
			return clsize->at(i);
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
		TVector3 GetPosition(int itr){
			TVector3 pos;
			pos =  tpc::getPosition(GetPadID(itr));
			pos.SetY(GetDL(itr));
			return pos;
		}
		TVector3 GetClusterPosition(int itr){
			TVector3 pos;
			pos=TVector3(clxTpc->at(itr),clyTpc->at(itr),clzTpc->at(itr)); 
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
		TVector3 GetRTheta(int padID);
		TVector2 GetLayerRow(int padID);
		int GetBCnt(){
			return ntBcOut;
		}
		Track GetTrack(int it=0){
			return Track(x0BcOut->at(it),y0BcOut->at(it),u0BcOut->at(it),v0BcOut->at(it));
		}
		


		int WhichEvent();
		void AssignG4Event(short * x,short* y,short* z,double* dedx);
		void AssignG4EventD(int* trkid,int* pid, double * x,double* y,double* z,double* dedx);
		int AssignRealEvent(double * x,double* y,double* z,double* dedx);
		void FillEvent();
		int NumberOfTracks(int min_points=6);
};





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
#endif
