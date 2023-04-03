#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "TPCGlobalFunctions.hh"
#include "Kinematics.hh"
#include "GeantTrack.hh"
#ifndef GeantManager_h
#define GeantManager_h
class GeantManager:public FileManager{
	protected:
		static const int MaxHits = 500;
		int evnum;
		//TPC
		int nttpc = -1;
		int nhittpc;
		int ntrk[MaxHits];
		int ititpc[MaxHits];
		int idtpc[MaxHits];
		double xtpc[MaxHits];
		double ytpc[MaxHits];
		double ztpc[MaxHits];
		double x0tpc[MaxHits];
		double y0tpc[MaxHits];
		double z0tpc[MaxHits];
		double pxtpc[MaxHits];
		double pytpc[MaxHits];
		double pztpc[MaxHits];
		double dedxtpc[MaxHits];
		double edeptpc[MaxHits];
		int ipadtpc[MaxHits];
		int laytpc[MaxHits];
		int rowtpc[MaxHits];
		//TPC
		//HTOF
		int nhHtof;
		int tidHtof[MaxHits];
		int pidHtof[MaxHits];
		int qHtof[MaxHits];
		double xHtof[MaxHits];
		double yHtof[MaxHits];
		double zHtof[MaxHits];
		double pxHtof[MaxHits];
		double pyHtof[MaxHits];
		double pzHtof[MaxHits];
		double deHtof[MaxHits];
		double tHtof[MaxHits];
		double lengthHtof[MaxHits];
		//HTOF
		//FTOF
		int nhFtof;
		int tidFtof[MaxHits];
		int pidFtof[MaxHits];
		int qFtof[MaxHits];
		double xFtof[MaxHits];
		double yFtof[MaxHits];
		double zFtof[MaxHits];
		double pxFtof[MaxHits];
		double pyFtof[MaxHits];
		double pzFtof[MaxHits];
		double deFtof[MaxHits];
		double tFtof[MaxHits];
		double lengthFtof[MaxHits];
		//FTOF

	protected:
		vector<GeantTrack> Tracks;
		void Initialize();
	public:
		GeantManager(){};
		virtual void LoadFile(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadChain("TPC_g");
			Initialize();
		}
		void SetEvent(int ievt){
			DataChain->GetEntry(ievt);
		}

		int GetNtTPC();


		int GetNhToF(){
			return nhFtof;
		}
		double GetPToF(int i ){
			return sqrt(pxFtof[i]*pxFtof[i]+pyFtof[i]*pyFtof[i]+pzFtof[i]*pzFtof[i])*1./1000.;
		}
		double GetQToF(int i){
			return qFtof[i];
		}
		double GetM2(int itof){
			return MassSquare(GetPToF(itof),lengthFtof[itof],tFtof[itof]);
		}


		void ConstructTracks();
		int GetNtracks();
		GeantTrack GetTrack(int it);
};


#endif
