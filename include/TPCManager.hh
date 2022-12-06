#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
#include "TPCGlobalFunctions.hh"
#include "Track.hh"
#include "TPCCluster.hh"
#ifndef TPCManager_h
#define TPCManager_h
TVector3 VertexPointHelix(const Double_t par1[5], const Double_t par2[5],
                 Double_t& dist, Double_t& t1, Double_t& t2)
{
  //helix function 1
  //x = [0] + [3]*cos(t);
  //y = [1] + [3]*sin(t);
  //z = [2] + [3]*[4]*t;

  //helix function 2
  //x = [5] + [8]*cos(t);
  //y = [6] + [8]*sin(t);
  //z = [7] + [8]*[9]*t;

  TF2 fvert_helix("fvert_helix",
                         "pow(([0]+[3]*cos(x))-([5]+[8]*cos(y)),2)+pow(([1]+[3]*sin(x))-([6]+[8]*sin(y)),2)+pow(([2]+[3]*[4]*x)-([7]+[8]*[9]*y),2)",
                         -5.,5.,-5.,5.);

  fvert_helix.SetParameter(0, par1[0]);
  fvert_helix.SetParameter(1, par1[1]);
  fvert_helix.SetParameter(2, par1[2]);
  fvert_helix.SetParameter(3, par1[3]);
  fvert_helix.SetParameter(4, par1[4]);
  fvert_helix.SetParameter(5, par2[0]);
  fvert_helix.SetParameter(6, par2[1]);
  fvert_helix.SetParameter(7, par2[2]);
  fvert_helix.SetParameter(8, par2[3]);
  fvert_helix.SetParameter(9, par2[4]);

  Double_t close_zin, close_zout;
  fvert_helix.GetMinimumXY(close_zin, close_zout);
  t1 = close_zin;
  t2 = close_zout;
  dist = TMath::Sqrt(fvert_helix.GetMinimum());

  Double_t xin = par1[0]+par1[3]*cos(close_zin);
  Double_t xout = par2[0]+par2[3]*cos(close_zout);
  Double_t yin =  par1[1]+par1[3]*sin(close_zin);
  Double_t yout = par2[1]+par2[3]*sin(close_zout);
  Double_t zin = par1[2]+par1[3]*par1[4]*close_zin;
  Double_t zout =  par2[2]+par2[3]*par2[4]*close_zout;

  // Double_t vx = (par1[0]+par1[3]*cos(close_zin) + par2[0]+par2[3]*cos(close_zout))/2.;
  // Double_t vy = (par1[1]+par1[3]*sin(close_zin) + par2[1]+par2[3]*sin(close_zout))/2.;
  // Double_t vz = (par1[2]+par1[3]*par1[4]*close_zin + par2[2]+par2[3]*par2[4]*close_zout)/2.;
  Double_t vx = (xin+xout)/2.;
  Double_t vy = (yin+yout)/2.;
  Double_t vz = (zin+zout)/2.;

  Double_t dist2 = sqrt(pow(xin-xout,2)
			+pow(yin-yout,2)
			+pow(zin-zout,2));
  // std::cout<<"dist ="<<dist<<", dist2="<<dist2<<std::endl;
  // std::cout<<"close_zin="<<close_zin<<", close_zout="<<close_zout<<std::endl;
  dist = dist2;
  Double_t vertx = -1.*vx;
  Double_t verty = vz;
  Double_t vertz = vy + ZTarget;
  return TVector3(vertx, verty, vertz);
}
class TPCManager:public FileManager{
	protected:
		TFile* hist_file;
		TH2Poly* PadHist=nullptr;
		TH2D* ZYHist=nullptr;
		TH2I* FlatHist=nullptr;
		TH2D* PosHist=nullptr;
		TGeoVolume *TPC3D;
		TPolyMarker3D *tpcHit3d;
		TCanvas* TPCCanv;
		vector<TPCHit> m_Hits;
		vector<TPCCluster> m_Clusters;
		TEllipse* HelixTrack[20] ;
		vector<vector<TLine*>> HelixTrackZY ;
		vector<int> *padTpc = new vector<int>;
		int iPadtpc[nhtpcmax];
		double xtpc[nhtpcmax];
		double ytpc[nhtpcmax];
		double ztpc[nhtpcmax];
		double dedxtpc[nhtpcmax];
		int idtpc[nhtpcmax];
		int ititpc[nhtpcmax];
		int ntrk[nhtpcmax];
		int ntTpc;
		vector<double>* dlTpc = new vector<double>;
		vector<double>* deTpc = new vector<double>;
		vector<double>* clxTpc = new vector<double>;
		vector<double>* clyTpc = new vector<double>;
		vector<double>* clzTpc = new vector<double>;
		vector<double>* helcxTpc = new vector<double>;
		vector<double>* helcyTpc = new vector<double>;
		vector<double>* helz0Tpc = new vector<double>;
		vector<double>* helrTpc = new vector<double>;
		vector<double>* heldzTpc = new vector<double>;
vector<double>* vtx = new vector<double>;
vector<double>* vty = new vector<double>;
vector<double>* vtz = new vector<double>;
		vector<int>* clsize = new vector<int>;
		int evnum,runnum;
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
		int GetEvnum(){return evnum;}
		int GetRunnum(){return runnum;}
		void SetEvent(int evt){
				x0BcOut->clear();
				y0BcOut->clear();
				u0BcOut->clear();
				v0BcOut->clear();
				clxTpc->clear();
				clyTpc->clear();
				clzTpc->clear();
				padTpc->clear();
				helcxTpc->clear();
				helcyTpc->clear();
				helz0Tpc->clear();
				helrTpc->clear();
				heldzTpc->clear();
				evnum=-1;
//			dlTpc->clear();
//			deTpc->clear();
			DataChain->GetEntry(evt);
			for(int i=0;i<20;++i){
//				if(HelixTrack[i]) delete HelixTrack[i];
			}
		};
		void InitializeHelix(){
			cout<<"InitializeHelix: "<<endl;
			for(auto h : HelixTrackZY){
				for(auto ht:h){
//					if(ht) delete ht;
				}
			}
			for(auto h : HelixTrack){
	//			if(h) delete h;
			}

			HelixTrackZY.clear();		
			HelixTrackZY.resize(ntTpc);		
			for(int ih = 0; ih< ntTpc;++ih){
				
				double cx = helcxTpc->at(ih);
				double cy = helcyTpc->at(ih);
				double z0 = helz0Tpc->at(ih);
				double r = helrTpc->at(ih);
				double dz = heldzTpc->at(ih);
				cout<<Form("Params = (%f,%f,%f,%f,%f)",cx,cy,r,z0,dz)<<endl;
//				TString title = Form("Helix%d",ih);
				HelixTrack[ih] = new TEllipse(cy+tpc::ZTarget,-cx,r,r);
				HelixTrack[ih]-> SetLineColor(kRed);
				HelixTrack[ih]-> SetFillStyle(0);
				HelixTrack[ih]-> SetLineColor(ih+1);
				HelixTrack[ih]-> SetLineWidth(2);
				for(int ip=0;ip<300;ip++){
					double t1 = 2*3.14*(ip/150.-0.5),t2=2*3.14*((ip+1)/150.-0.5);
					double y1 = r*dz*t1+z0,y2 = r*dz*t2+z0;
					double z1 = r*sin(t1)+cy+tpc::ZTarget,z2 = r*sin(t2)+cy+tpc::ZTarget;
					/*
					double y1 = -300+ip*600./100;
					double y2 = -300+(ip+1)*600./1000;
					double t1 = (y1-z0)/r*dz;
					double t2 = (y2-z0)/r*dz;
					double z1 = r*cos(t1)+cy+tpc::ZTarget;
					double z2 = r*cos(t2)+cy+tpc::ZTarget;
					*/
					HelixTrackZY[ih].push_back(new TLine(z1,y1,z2,y2));
					HelixTrackZY[ih].at(ip)->SetLineColor(ih+1);
					HelixTrackZY[ih].at(ip)->SetLineWidth(2);
				}
			}
		}
		void DrawHelix(){
			for(int ih = 0; ih< ntTpc;++ih){
				HelixTrack[ih]->Draw("psame");
			}
		}
		void DrawHelixZY(){
			for(int ih = 0; ih< ntTpc;++ih){
				for(int ip=0;ip<300;ip++){
					HelixTrackZY[ih].at(ip)->Draw("same");
				}
			}
		}
		void SearchVertex();
		//Histogram Methods//
		void InitializeHistograms();
		void InitializeTPC();
		void SetTitle(TString title){
			PadHist->SetTitle(title);
			ZYHist->SetTitle(title+"ZY");
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
			ZYHist->Reset("");
		}
		void ClearTPC(){
//			TPC3D->Reset("");
		}
		TH2Poly* GetPadHistogram(){
			return PadHist;
		}
		TH2D* GetZYHistogram(){
			return ZYHist;
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



		virtual void Process(double* vals){};


		void DrawTPC(){
	//		auto* dir = gDirectory()->cd();
			TPCCanv = new TCanvas("c1","c1",1200,600);
			TPCCanv->cd();
			TView3D *view = (TView3D*) TView::CreateView(1);
			TPC3D->Draw("");
//			dir->cd();
		}
		int GetNhits(int clusters){
			if(!clusters)	return Min(padTpc->size(),max_nh);//Min(nhittpc,max_nh);
			else 					return Min(clsize->size(),max_nh);
		};
		int GetPadID(int i){
			cout<<i<<" : GetPadID"<<endl;
			if(i<padTpc->size()){
				return padTpc->at(i);
			}
			else return 0;
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
			if(!cluster){
				pos =  tpc::getPosition(GetPadID(itr));
				pos.SetY(GetDL(itr));
			}
			else{
				pos = TVector3(clxTpc->at(itr),clyTpc->at(itr),clzTpc->at(itr)); 
			}
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
	ZYHist = new TH2D("ZYHist","ZYHist",40,-250,250,50,-300,300);
}
void TPCManager::InitializeTPC(){
	TPC3D=TPCGeometry();
}
void TPCManager::SearchVertex(){
		for(int nt1 = 0; nt1<ntTpc;++nt1){
			double hcx = helcxTpc->at(nt1);
			double hcy = helcyTpc->at(nt1);
			double hz0 = helz0Tpc->at(nt1);
			double hr = helrTpc->at(nt1);
			double hdz = heldzTpc->at(nt1);
			double par1[5] = {hcx,hcy,hz0,hr,hdz};
			TVector3 vert(vtx->at(nt1),vty->at(nt1),vtz->at(nt1));
			for(int nt2 = 0; nt2<ntTpc;++nt2){
				if(nt2 <= nt1) continue;
				double hcx2 = helcxTpc->at(nt2);
				double hcy2 = helcyTpc->at(nt2);
				double hz02 = helz0Tpc->at(nt2);
				double hr2 = helrTpc->at(nt2);
				double hdz2 = heldzTpc->at(nt2);
				double par2[5] = {hcx2,hcy2,hz02,hr2,hdz2};
				double cd,t1,t2;
				auto vert = VertexPointHelix(par1,par2,cd,t1,t2);
//				hist->Fill(cd);
				//			TVector3 mom_vtx = tp->CalcHelixMom(par1, vert.y());
			}
		}
}
#endif
