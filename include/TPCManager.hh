#include "Utils.hh"
#include "FileManager.hh"
#include "TPCPadHelper.hh"
#include "PhysicalConstants.hh"
#include "TPCGlobalFunctions.hh"
//#include "Track.hh"
#include "TPCCluster.hh"
#include "ReconTools.hh"
#ifndef TPCManager_h
#define TPCManager_h
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
vector<double>* chisqr = new vector<double>;
vector<double>* helix_cx = new vector<double>;
vector<double>* helix_cy = new vector<double>;
vector<double>* helix_z0 = new vector<double>;
vector<double>* helix_r = new vector<double>;
vector<double>* helix_dz = new vector<double>;
vector<double>* vtx = new vector<double>;
vector<double>* vty = new vector<double>;
vector<double>* vtz = new vector<double>;
 vector<int>* isBeam = new vector<int>;
vector<int>* clsize = new vector<int>;
vector<int>* pid = new vector<int>;
vector<int>* charge = new vector<int>;
		int evnum,runnum;
		int htofnhits;
		int htofhitpat[34];
		int nhittpc; double htofua[34];
		int gp = 0;
		int gpb = 0;
		bool cluster = false;
		int ntBcOut=0;
		int ntk=0;
		vector<double>* x0BcOut = new vector<double>;
		vector<double>* y0BcOut = new vector<double>;
		vector<double>* u0BcOut = new vector<double>;
		vector<double>* v0BcOut = new vector<double>;
		bool ldflg,xiflg;
		Recon Ld,Xi;
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
				helix_cx->clear();
				helix_cy->clear();
				helix_z0->clear();
				helix_r->clear();
				helix_dz->clear();
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
				
				double cx = helix_cx->at(ih);
				double cy = helix_cy->at(ih);
				double z0 = helix_z0->at(ih);
				double r = helix_r->at(ih);
				double dz = helix_dz->at(ih);
	//			cout<<Form("Params = (%f,%f,%f,%f,%f)",cx,cy,r,z0,dz)<<endl;
//				TString title = Form("Helix%d",ih);
//				if(r>4000) continue;
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
		void ReconEvent();
		void DrawHelix(){
			for(int ih = 0; ih< ntTpc;++ih){
				HelixTrack[ih]->Draw("psame");
			}
		}
		void DrawVertex(){
			TVector3 LV,XV;
			bool ldflg = Ld.Exist(),xiflg = Xi.Exist();
			double z1=0,z2=0,x1=0,x2=0,z3=0,x3=0;
			TEllipse* ldvert;TEllipse* xivert;
			if(ldflg) {
				LV = Ld.Vertex();
				z1 = LV.Z();
				x1 = LV.X();
				auto Dir = Ld.Momentum();
				Dir = Dir* (1./Dir.Mag());
				z2 = LV.Z()-Dir.Z()*100;
				x2 = LV.X()-Dir.X()*100;
				cout<<Form("Lambda Vertex (%f,%f,%f) Mass: %f",LV.X(),LV.Y(),LV.Z(),Ld.Mass())<<endl;
			}
			if(xiflg){
				XV = Xi.Vertex();
				z3 = XV.Z();
				x3 = XV.X();
				cout<<Form("Xi Vertex (%f,%f,%f) Mass: %f",XV.X(),XV.Y(),XV.Z(),Xi.Mass())<<endl;
			cout<<"PropDist : "<<(LV-XV).Mag()<<endl;
			}
			auto ld = new TLine(z1,x1,z2,x2);	
			ld->SetLineWidth(3);ld->SetLineColor(kMagenta);
			ld->Draw("psame");
			ldvert = new TEllipse(z1,x1,3,0);
			ldvert->SetLineColor(kMagenta);
			ldvert->Draw("psame");
			if(xiflg){
				auto ld2 = new TLine(z1,x1,z3,x3);	
			ld2->SetLineWidth(3);ld2->SetLineColor(kCyan);
				ld2->Draw("psame");
				xivert = new TEllipse(z3,x3,3,0);
				xivert->SetLineColor(kCyan);
				xivert->Draw("psame");
			}
		}
		void DrawVertexZY(){
			TVector3 LV,XV;
			TEllipse* ldvert;TEllipse* xivert;
			bool ldflg = Ld.Exist(),xiflg = Xi.Exist();
			double z1=0,z2=0,y1=0,y2=0,z3=0,y3=0;
			if(ldflg) {
				LV = Ld.Vertex();
				z1 = LV.Z();
				y1 = LV.Y();
				auto Dir = Ld.Momentum();
				Dir = Dir* (1./Dir.Mag());
				z2 = LV.Z()-Dir.Z()*100;
				y2 = LV.Y()-Dir.Y()*100;
			}
			if(xiflg){
				XV = Xi.Vertex();
				z3 = XV.Z();
				y3 = XV.Y();
			}
			auto ld = new TLine(z1,y1,z2,y2);	
			ld->SetLineWidth(3);ld->SetLineColor(kMagenta);
			ld->Draw("psame");
			ldvert = new TEllipse(z1,y1,3,0);
			ldvert->SetLineColor(kMagenta);
			ldvert->Draw("psame");
			if(xiflg){
			auto ld2 = new TLine(z1,y1,z3,y3);	
			ld2->SetLineWidth(3);ld2->SetLineColor(kCyan);
			ld2->Draw("psame");
			xivert = new TEllipse(z3,y3,3,0);
			xivert->SetLineColor(kCyan);
			xivert->Draw("psame");
			}
		}
		void DrawHelixZY(){
			for(int ih = 0; ih< ntTpc;++ih){
				for(int ip=0;ip<300;ip++){
					HelixTrackZY[ih].at(ip)->Draw("same");
				}
			}
		}
		bool LambdaEvent(){return Ld.Exist();}
		bool XiEvent(){return Xi.Exist();}
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
#if 0
		Track GetTrack(int it=0){
			return Track(x0BcOut->at(it),y0BcOut->at(it),u0BcOut->at(it),v0BcOut->at(it));
		}
#endif





		int WhichEvent();
		void AssignG4Event(short * x,short* y,short* z,double* dedx);
		void AssignG4EventD(int* trkid,int* pid, double * x,double* y,double* z,double* dedx);
		int AssignRealEvent(double * x,double* y,double* z,double* dedx);
		void FillEvent();
//		int NumberOfTracks(int min_points=6);
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
			double hcx = helix_cx->at(nt1);
			double hcy = helix_cy->at(nt1);
			double hz0 = helix_z0->at(nt1);
			double hr = helix_r->at(nt1);
			double hdz = helix_dz->at(nt1);
			double par1[5] = {hcx,hcy,hz0,hr,hdz};
			TVector3 vert(vtx->at(nt1),vty->at(nt1),vtz->at(nt1));
			for(int nt2 = 0; nt2<ntTpc;++nt2){
				if(nt2 <= nt1) continue;
				double hcx2 = helix_cx->at(nt2);
				double hcy2 = helix_cy->at(nt2);
				double hz02 = helix_z0->at(nt2);
				double hr2 = helix_r->at(nt2);
				double hdz2 = helix_dz->at(nt2);
				double par2[5] = {hcx2,hcy2,hz02,hr2,hdz2};
				double cd,t1,t2;
				auto vert = VertexPointHelix(par1,par2,cd,t1,t2);
//				hist->Fill(cd);
				//			TVector3 mom_vtx = tp->CalcHelixMom(par1, vert.y());
			}
		}
}
#endif
