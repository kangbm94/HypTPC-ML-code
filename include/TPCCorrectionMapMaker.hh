#include "TPCManager.hh"


class TPCCorrectionMapMaker: public TPCManager{
	protected:
		static const int max_hit = 1000;
		static const int nbin_x=17;
		static const int nbin_y=17;
		static const int nbin_z=21;
		double x1=-40,x2=40,y1=-40,y2=40,z1=-250,z2=250;


		int hbin_x=20,hbin_y=20,hbin_z=20;
		double cx=(x2-x1)/(nbin_x-1);
		double cy=(y2-y1)/(nbin_y-1);
		double cz=(z2-z1)/(nbin_z-1);
		TH1D* hist_x[nbin_x][nbin_y][nbin_z];
		TH1D* hist_y[nbin_x][nbin_y][nbin_z];
		TH1D* hist_z[nbin_x][nbin_y][nbin_z];
		double peakx[nbin_x][nbin_y][nbin_z];	
		double peaky[nbin_x][nbin_y][nbin_z];	
		double peakz[nbin_x][nbin_y][nbin_z];	
		vector<double>* xCorVec = new vector<double>;
		vector<double>* yCorVec = new vector<double>;
		vector<double>* zCorVec = new vector<double>;
		int ntBcOut=0;
		
	public:
		TPCCorrectionMapMaker(){
			for(int ix=0;ix<nbin_x;++ix){
				for(int iy=0;iy<nbin_y;++iy){
					for(int iz=0;iz<nbin_z;++iz){
						TString title = Form("hist%d_%d_%d_",ix,iy,iz);
						hist_x[ix][iy][iz]= new TH1D(title+"x",title+"x",hbin_x,-cx/2,cx/2);
						hist_y[ix][iy][iz]= new TH1D(title+"y",title+"y",hbin_y,-cy/2,cy/2);
						hist_z[ix][iy][iz]= new TH1D(title+"z",title+"z",hbin_z,-cz/2,cz/2);
					}
				}
			}
			cout<<"Histogram initialized"<<endl;
		}
		void WriteHist(){
			for(int ix=0;ix<nbin_x;++ix){
				for(int iy=0;iy<nbin_y;++iy){
					for(int iz=0;iz<nbin_z;++iz){
						hist_x[ix][iy][iz]->Write();
						hist_y[ix][iy][iz]->Write();
						hist_z[ix][iy][iz]->Write();
					}
				}
			}
		}
		void LoadHist(TString FileName){
			TFile* file = new TFile(FileName,"read");
			for(int ix=0;ix<nbin_x;++ix){
				for(int iy=0;iy<nbin_y;++iy){
					for(int iz=0;iz<nbin_z;++iz){
						TString title = Form("hist%d_%d_%d_",ix,iy,iz);
						delete hist_x[ix][iy][iz];delete hist_y[ix][iy][iz];delete hist_z[ix][iy][iz];
						hist_x[ix][iy][iz]=(TH1D*)file->Get(title+"x");
						hist_y[ix][iy][iz]=(TH1D*)file->Get(title+"y");
						hist_z[ix][iy][iz]=(TH1D*)file->Get(title+"z");
						peakx[ix][iy][iz]=GetPeakPosition(hist_x[ix][iy][iz]);
						peaky[ix][iy][iz]=GetPeakPosition(hist_y[ix][iy][iz]);
//						peakz[ix][iy][iz]=GetPeakPosition(hist_z[ix][iy][iz]);
						peakz[ix][iy][iz]=0;
							//GetPeakPosition(hist_z[ix][iy][iz]);
							//hist_x[ix][iy][iz]->GetBinCenter(hist_x[ix][iy][iz]->GetMaximum());
//						peaky[ix][iy][iz]=hist_y[ix][iy][iz]->GetBinCenter(hist_y[ix][iy][iz]->GetMaximum());
//						peakz[ix][iy][iz]=hist_z[ix][iy][iz]->GetBinCenter(hist_z[ix][iy][iz]->GetMaximum());
					}
				}
			}
		}
		void WriteParam(){
			for(int ix=0;ix<nbin_x;++ix){
				for(int iy=0;iy<nbin_y;++iy){
					for(int iz=0;iz<nbin_z;++iz){
						//vector<double>par= {BinPosX(ix),BinPosY(iy),BinPosZ(iz),peakx[ix][iy][iz],peaky[ix][iy][iz],peakz[ix][iy][iz]};
						vector<double>par= {BinPosX(ix),BinPosY(iy),BinPosZ(iz),peakx[ix][iy][iz],peaky[ix][iy][iz],peakz[ix][iy][iz]};
//						cout<<peakz[ix][iy][iz]<<endl;
						vector<int>dum={};
						WriteParameter(dum,par);
					}
				}
			}
		}
		virtual void LoadTPCBcOut(TString FileName){ DataFile = new TFile(FileName,"READ");
			cout<<FileName<<" Opened"<<endl;
			LoadTPCBcOutChain("tpc");
		}
		void MakeCorParameterFile(TString FileName){
			MakeParameterFile(FileName);
			vector<int> nb = {nbin_x,nbin_y,nbin_z};
			vector<double> par = {x1,y1,z1,cx,cy,cz};
			WriteParameter(nb,par);
		}
		void LoadTPCBcOutChain(TString ChainName);
		int Bin(int nbin,double r1,double r2,double r){
			double dr = (r2-r1)/(nbin-1);
			r-=(r1-dr/2);
			return r/dr;
		}
		int BinX(double x){
			return Bin(nbin_x,x1,x2,x);
		}
		int BinY(double y){
			return Bin(nbin_y,y1,y2,y);
		}
		int BinZ(double z){
			return Bin(nbin_z,z1,z2,z);
		}
		double BinPos(int nbin, double r1,double r2, int bin){
			double dr = r2-r1;
			return r1+bin*dr/(nbin-1);
		}
		double BinPosX(double xbin){
			return BinPos(nbin_x,x1,x2,xbin);
		}
		double BinPosY(double ybin){
			return BinPos(nbin_y,y1,y2,ybin);
		}
		double BinPosZ(double zbin){
			return BinPos(nbin_z,z1,z2,zbin);
		}
		bool Accept(double x,double y, double z){
			int xb=BinX(x),yb=BinY(y),zb=BinZ(z);
			bool xf = xb<nbin_x&&xb>-1,yf = yb<nbin_y&&yb>-1,zf = zb<nbin_z&&zb>-1;
			return xf&&yf&&zf;
		}
		void FillHistograms(TVector3 Pos, TVector3 Cor);
		void AssignHit(TVector3& Pos, TVector3& Cor,int nh);
		void AssignHits(TVector3* Poss, TVector3* Cors);
		void Process();
};
void TPCCorrectionMapMaker::LoadTPCBcOutChain(TString ChainName ){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	cout<<"Chain Loaded: "<<DataChain->GetEntries()<<endl;
	DataChain->SetBranchAddress("cluster_hitpos_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_hitpos_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_hitpos_z",&clzTpc);
	DataChain->SetBranchAddress("xCorVec",&xCorVec);
	DataChain->SetBranchAddress("yCorVec",&yCorVec);
	DataChain->SetBranchAddress("zCorVec",&zCorVec);
	DataChain->SetBranchAddress("cluster_de",&deTpc);
	DataChain->SetBranchAddress("cluster_size",&clsize);
	DataChain->SetBranchAddress("ntBcOut",&ntBcOut);
}
void TPCCorrectionMapMaker::FillHistograms(TVector3 Pos, TVector3 Cor){
	double x=Pos.X(),y=Pos.Y(),z=Pos.Z();
	double Corx=Cor.X(),Cory=Cor.Y(),Corz=Cor.Z();
	if(Accept(x,y,z)&&ntBcOut==1){
		hist_x[BinX(x)][BinY(y)][BinZ(z)]->Fill(Corx);
		hist_y[BinX(x)][BinY(y)][BinZ(z)]->Fill(Cory);
		hist_z[BinX(x)][BinY(y)][BinZ(z)]->Fill(Corz);
	}
}
