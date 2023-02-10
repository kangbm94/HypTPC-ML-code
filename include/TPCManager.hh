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
double ExplicitHelix(double x,double y,double cx,double cy,double z0,double dz){
	double x0 = x-cx;
	double y0 = y-cy;
	double r = sqrt(x0*x0+y0*y0);
	return 0;
}
Double_t HypTPCdEdx(Double_t Z, Double_t *x, Double_t *p){
  //x : poq
  //p[0] : converting constant p[1] : density effect correction p[2] : mass
  Double_t me  = 0.5109989461;
  Double_t rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672); //[g cm-3]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t ZoverA = 17.2/37.6; //[mol g-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
  Double_t I2 = 0.9*188.0 + 0.1*41.7; I2 = I2*I2; //Mean excitaion energy [eV]
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t mom = 1000.*x[0]*Z; //MeV
  Double_t beta2 = mom*mom/(mom*mom+p[2]*p[2]);
  Double_t gamma2 = 1./(1.-beta2);
  Double_t Wmax = 2*me*beta2*gamma2/((me+p[2])*(me+p[2])+2*me*p[2]*(TMath::Sqrt(gamma2)-1));
  Double_t dedx = p[0]*constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2)-beta2-p[1]);
  return dedx;
}

Double_t HypTPCBethe(Double_t *x, Double_t *p){ return HypTPCdEdx(1, x, p); }
Int_t HypTPCdEdxPID_temp(Double_t dedx, Double_t poq){
  Double_t bethe_par[2] = {7195.92, -10.5616};
  Double_t limit = 0.6; //GeV/c
  Double_t mpi = 139.57039;
  Double_t mk  = 493.677;
  Double_t mp  = 938.2720813;
  Double_t md  = 1875.612762;
  TF1 *f_pim = new TF1("f_pim", HypTPCBethe, -3., 0., 3);
  TF1 *f_km = new TF1("f_km", HypTPCBethe, -3., 0., 3);
  TF1 *f_pip = new TF1("f_pip", HypTPCBethe, 0., 3., 3);
  TF1 *f_kp = new TF1("f_kp", HypTPCBethe, 0., 3., 3);
  TF1 *f_p = new TF1("f_p", HypTPCBethe, 0., 3., 3);
  TF1 *f_d = new TF1("f_d", HypTPCBethe, 0., 3., 3);

  f_pim -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_km -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_pip -> SetParameters(bethe_par[0], bethe_par[1], mpi);
  f_kp -> SetParameters(bethe_par[0], bethe_par[1], mk);
  f_p -> SetParameters(bethe_par[0], bethe_par[1], mp);
  f_d -> SetParameters(bethe_par[0], bethe_par[1], md);

  Int_t pid[3] = {0};
  if(poq >= limit){
    pid[0]=1; pid[1]=1; pid[2]=1;
  }
  else if(limit > poq && poq >= 0.){
    Double_t dedx_d = f_d -> Eval(poq); Double_t dedx_p = f_p -> Eval(poq);
    Double_t dedx_kp = f_kp -> Eval(poq); Double_t dedx_pip = f_pip -> Eval(poq);
    if(dedx_d > dedx && dedx >= dedx_kp) pid[2]=1;
    if(dedx_p > dedx){
      pid[0]=1; pid[1]=1;
    }
  }
  else if(0.> poq && poq >= -limit){
    pid[0]=1; pid[1]=1;
  }
  else{
    pid[0]=1; pid[1]=1;
  }

  delete f_pim;
  delete f_km;
  delete f_pip;
  delete f_kp;
  delete f_p;
  delete f_d;

  Int_t output = pid[0] + pid[1]*2 + pid[2]*4;
  return output;
}

class TPCManager:public FileManager{
	protected:
		TFile* hist_file;
		TH2Poly* PadHist=nullptr;
		TH2Poly* PadHistf=nullptr;
		TH2Poly* PadHistb=nullptr;
		TH2D* ZYHist=nullptr;
		TH2D* ZYHistf=nullptr;
		TH2D* ZYHistb=nullptr;
		TH2I* FlatHist=nullptr;
		TH2D* PosHist=nullptr;
		TH1D* YHist=nullptr;
		TGeoVolume *TPC3D;
		TPolyMarker3D *tpcHit3d;
		TCanvas* TPCCanv;
		vector<TPCHit> m_Hits;
		vector<TPCCluster> m_Clusters;
		TEllipse* HelixTrack[20] ;
		TEllipse* AccidentalTrack[20] ;
		vector<int> *trigflag = new vector<int>;
		vector<vector<TLine*>> HelixTrackZY ;
		vector<vector<TLine*>> AccidentalTrackZY ;
		vector<int> *padTpc = new vector<int>;
		int iPadtpc[nhtpcmax];
		double xtpc[nhtpcmax];
		double ytpc[nhtpcmax];
		double ztpc[nhtpcmax];
		double dedxtpc[nhtpcmax];
		int idtpc[nhtpcmax];
		int ititpc[nhtpcmax];
		int ntrk[nhtpcmax];
		int ntAcc;
		int ntTpc;
		int nclTpc;
		int nclfTpc;
		int nclbTpc;
		vector<double>* dlTpc = new vector<double>;
		vector<double>* deTpc = new vector<double>;
		vector<double>* cldeTpc = new vector<double>;
		vector<double>* clxTpc = new vector<double>;
		vector<double>* clyTpc = new vector<double>;
		vector<double>* clzTpc = new vector<double>;
		vector<double>* resolution_x = new vector<double>;
		vector<double>* resolution_y = new vector<double>;
		vector<double>* resolution_z = new vector<double>;
		vector<double>* resolution = new vector<double>;
		vector<double>* chisqr = new vector<double>;
		vector<double>* accidental_cx = new vector<double>;
		vector<double>* accidental_cy = new vector<double>;
		vector<double>* accidental_z0 = new vector<double>;
		vector<double>* accidental_r = new vector<double>;
		vector<double>* accidental_dz = new vector<double>;
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
		vector<int>* hough_flag = new vector<int>;
		vector<double>* hough_dist = new vector<double>;
		vector<double>* dEdx = new vector<double>;
		vector<double>* mom0 = new vector<double>;
		vector<int>* pid = new vector<int>;
		vector<int>* charge = new vector<int>;
		vector<double>* beam_y = new vector<double>;
		vector<double>* beam_p0 = new vector<double>;
		vector<double>* beam_p1 = new vector<double>;
		vector<double>* beam_p2 = new vector<double>;
		vector<double>* beam_v = new vector<double>;
		vector<vector<double>>* helix_t = new vector<vector<double>>;
		vector<vector<double>>* track_cluster_layer = new vector<vector<double>>;
		vector<vector<double>>* track_cluster_x_center = new vector<vector<double>>;
		vector<vector<double>>* track_cluster_y_center = new vector<vector<double>>;
		vector<vector<double>>* track_cluster_z_center = new vector<vector<double>>;
		int evnum,runnum;
		int htofnhits;
		int htofhitpat[34];
		int nhittpc; double htofua[34];
		int gp = 0;
		int gpb = 0;
		bool cluster = false;
		int ntBcOut=0;
		int ntk=0;
		double npts=100;
		double anpts=5;
		vector<double>* x0BcOut = new vector<double>;
		vector<double>* y0BcOut = new vector<double>;
		vector<double>* u0BcOut = new vector<double>;
		vector<double>* v0BcOut = new vector<double>;
		bool ldflg,xiflg;
		Recon Ld,Xi;
		TF1* f_bethe = new TF1("f_betaP",HypTPCBethe,0.1,3,3);

	Double_t bethe_pars[2] = {7195.92, -10.5616};
	Double_t mprt  = 938.2720813;
		
		TH2D* ZYHistsAcc[20];
		TH2D* CirHistsAcc[20];

	public:
		TPCManager(){};

		//Data I/O
		virtual void LoadFile(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadChain("tpc");
		}
		virtual void LoadClusterFile(TString FileName,TString chainname = "tpc"){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadClusterChain(chainname);
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
			evnum=-1;
			DataChain->GetEntry(evt);
		};
		int GetNTracks(){
			return ntTpc;
		}
		int GetNTracksAcc(){
			return ntAcc;
		}
		bool TagTrig(int i ){
			if(trigflag->at(i)>0){
				return true;
			}
			else{
				return false;
			}
		}
		void InitializeHelix();
		void InitializeAccidental();
		void ReconEvent();
		void DrawHelix();
		void DrawHelix(int it);
		void DrawHelixZY(); 
		void DrawHelixZY(int it); 
		void DrawAccidental();
		void DrawAccidental(int it);
		void DrawAccidentalZY(); 
		void DrawAccidentalZY(int it); 
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
			cout<<"DrawingVertexZY..."<<endl;
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
				xivert->SetLineColor(kCyan); xivert->Draw("psame"); 
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

		void FillAccHists();
		
		void FillHist(double z, double x);
		void LoadTPC3D();
		void FillHist(int itr);
		void FillAntiProtonHist();
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
		void ClearHistogram();
		void ClearTPC(){
			//			TPC3D->Reset("");
		}
		TH2Poly* GetPadHistogram(){
			return PadHist;
		}
		TH2Poly* GetPadHistogramf(){
			return PadHistf;
		}
		TH2Poly* GetPadHistogramb(){
			return PadHistb;
		}
		TH2D* GetZYHistogram(){
			return ZYHist;
		}
		TH2D* GetZYHistogramf(){
			return ZYHistf;
		}
		TH2D* GetZYHistogramb(){
			return ZYHistb;
		}
		TH1D* GetYHistogram(){
			return YHist;
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

		vector<double>* GetAccidentalDist();
		vector<double>* GetHoughDist(){
			return hough_dist;
		}
		vector<int>* GetHoughFlag(){
			return hough_flag;
		}


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
			else 					return Min(clxTpc->size(),max_nh);
		};
		int GetNhits(){
			return nclTpc;
		}
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
		int GetClDe(int i){
			return cldeTpc->at(i);
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
		void FillgHitPos(double peak, double width){
			int nh = GetNhits(1);
			for(int i=0;i<nh;++i){
				TVector3 hitv = GetPosition(i);
				double x = hitv.X();
				double y = hitv.Y();
				double z = hitv.Z();
				TVector3 res(1,1,1);
				if(abs(y-peak)<width) {gHitPos.push_back(hitv);gRes.push_back(res);}
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
	vector<double>GetBeamY(){
		return *beam_y;
	}
	vector<double>GetBeamP0(){
		return *beam_p0;
	}
	vector<double>GetBeamP1(){
		return *beam_p1;
	}
	vector<double>GetBeamP2(){
		return *beam_p2;
	}
	vector<double>GetBeamV(){
		return *beam_v;
	}
	vector<double>GetAccidentalR(){
		return * accidental_r;
	}
	vector<double>GetAccidentalZ0(){
		return * accidental_z0;
	}
	vector<double>GetAccidentalDZ(){
		return * accidental_dz;
	}

	TH2D* GetZYHistAcc(int i);
	TH2D* GetCirHistAcc(int i);



		int WhichEvent();
		void AssignG4Event(short * x,short* y,short* z,double* dedx);
		void AssignG4EventD(int* trkid,int* pid, double * x,double* y,double* z,double* dedx);
		int AssignRealEvent(double * x,double* y,double* z,double* dedx);
		void FillEvent();
		//		int NumberOfTracks(int min_points=6);
		void SetBetheProton(){
			f_bethe ->SetParameters(bethe_pars[0],bethe_pars[1],mprt);
	}
		bool IsAntiProton( int it){
			if(charge->at(it)<0) return false;
			double dedx_p = f_bethe->Eval(mom0->at(it));
			if(dEdx->at(it)>dedx_p and mom0->at(it)<0.6){
				cout<<Form("mom : %f, dedx_p : %f",mom0->at(it),dedx_p)<<endl;
				return true;
			}
			else return false;
		}
};





void TPCManager::InitializeHistograms(){
	PadHist = tpc::InitializeHistogram();
	PadHistf = tpc::InitializeHistogram();
	PadHistb = tpc::InitializeHistogram();
	//	FlatHist = new TH2I("PadRTheta","PadRTheta",32,0,32,240,0,240);
	//	PosHist = new TH2D("PosHisto","PosHisto",128,-250,250,128,-250,250);
	ZYHist = new TH2D("ZYHist","ZYHist",40,-250,250,140,-350,350);
	YHist = new TH1D("YHist","YHist",140,-350,350);
	for(int i = 0; i<20;++i){
			TString titleZY = Form("HistZY_%d",i);
			TString titleCir = Form("HistCir_%d",i);
			ZYHistsAcc[i] = new TH2D(titleZY,titleZY,40,-250,250,50,-350,350);	
			CirHistsAcc[i] = new TH2D(titleCir,titleCir,50,-250,250,50,-250,250);	
	}
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
