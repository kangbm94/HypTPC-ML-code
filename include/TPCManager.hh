#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
#include "TPCGlobalFunctions.hh"
#include "Track.hh"
#include "TPCCluster.hh"
#ifndef TPCManager_h
#define TPCManager_h

class TPCManager:public FileManager{
	protected:
		TFile* hist_file;
		TH2Poly* PadHist=nullptr;
		TH2I* FlatHist=nullptr;
		TH2D* PosHist=nullptr;
		TGeoVolume *TPC3D;
		TPolyMarker3D *tpcHit3d;
		TCanvas* TPCCanv;
		vector<TPCHit> m_Hits;
		vector<TPCCluster> m_Clusters;

		vector<int> *padTpc = new vector<int>;
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
		void LoadTPCBcOut(TString FileName);
		void LoadChain(TString ChainName);
		void LoadClusterChain(TString ChainName);
		void LoadTPCBcOutChain(TString ChainName);
		void LoadG4Chain(TString ChainName);
		void LoadBcOut();
		int GetNEvent(){
			return DataChain->GetEntries();
		};
		int GetEntries(){return GetNEvent();}
		void SetEvent(int evt){
				x0BcOut->clear();
				y0BcOut->clear();
				u0BcOut->clear();
				v0BcOut->clear();
				clxTpc->clear();
				clyTpc->clear();
				clzTpc->clear();
//			padTpc->clear();
//			dlTpc->clear();
//			deTpc->clear();
			DataChain->GetEntry(evt);
		};

		//Histogram Methods//
		void InitializeHistograms();
		void InitializeTPC();
		void SetTitle(TString title){
			PadHist->SetTitle(title);
			FlatHist->SetTitle(title);
		}
		void FillHist(double z, double x);
		void LoadTPC3D();
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
			PadHist->Reset("");
		}
		void ClearTPC(){
//			TPC3D->Reset("");
		}
		TH2Poly* GetPadHistogram(){
			return PadHist;
		}

		void AssignHits();
		bool MakeUpClusters(double Vth);
		int GetNumberOfMHits(){return m_Hits.size();}
		int GetNumberOfMCls(){return m_Clusters.size();}
		void ClearHits(){
			m_Hits.clear();
			m_Clusters.clear();
		}
		TPCHit GetMHit(int i){return m_Hits[i];}
		TPCCluster GetMCl(int i){return m_Clusters[i];}



		virtual void Process(double* vals);


		void DrawTPC(){
	//		auto* dir = gDirectory()->cd();
			TPCCanv = new TCanvas("c1","c1",1200,600);
			TPCCanv->cd();
			TView3D *view = (TView3D*) TView::CreateView(1);
			TPC3D->Draw("");
//			dir->cd();
		}
		int GetNhits(int clusters){
			if(!clusters)	return Min(nhittpc,max_nh);
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
	PadHist = tpc::InitializeHistogram();
	//	FlatHist = new TH2I("PadRTheta","PadRTheta",32,0,32,240,0,240);
	//	PosHist = new TH2D("PosHisto","PosHisto",128,-250,250,128,-250,250);
}
void TPCManager::InitializeTPC(){
	TPC3D=TPCGeometry();
}
#endif
