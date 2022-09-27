#include "TPCManager.hh"
static const int nbin_x=19;
static const int nbin_y=7;
static const int snake_nbin_x = 100;
static const int snake_nbin_z = 40;
static const int nbin_z=snake_nbin_z;
class TPCEnergyCorrector: public TPCManager{
	protected:
		vector<TH1D*> DeHist(max_padid);
		TFile* OutFile=NULL;
		double de1=0,de2=500;
		TF1* flandau = new TF1("flandau","landau",de1,de2);

	public:
		TPCEnergyCorrector(){
		}
		double BinPosZ(double zbin){
			return BinPos(nbin_z,z1+(z2-z1)/(nbin_z)/2,z2-(z2-z1)/(nbin_z)/2,zbin);
		}
		void MakeOutFile(TString FileName){
			OutFile= new TFile(FileName,"recreate");}
		void WriteDeHist();
		void LoadDeHist(TString FileName);
		void WriteParam();
		TH2D* GetSnakeHist(int ix,int iy,bool isX);
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
		void ScanSnake(TH2D* hist,vector<double>& pars, vector<double>& ents); 


		void ParamFilter(vector<double>&pars,vector<int>ent);
		void SetParameter(vector<double>px,vector<double>py,int ix,int iy);
		bool Accept(double x,double y, double z){
			int xb=BinX(x),yb=BinY(y),zb=BinZ(z);
			bool xf = xb<nbin_x&&xb>-1,yf = yb<nbin_y&&yb>-1,zf = zb<nbin_z&&zb>-1;
			return xf&&yf&&zf;
		}
		void FillHists(TVector3 Pos, TVector3 Cor);
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
void TPCCorrectionMapMaker::FillHists(TVector3 Pos, TVector3 Cor){
	double x=Pos.X(),y=Pos.Y(),z=Pos.Z();
	double Corx=Cor.X(),Cory=Cor.Y(),Corz=Cor.Z();
	if(Accept(x,y,z)&&ntBcOut==1){
		hist_x[BinX(x)][BinY(y)][BinZ(z)]->Fill(Corx);
		hist_y[BinX(x)][BinY(y)][BinZ(z)]->Fill(Cory);
		hist_z[BinX(x)][BinY(y)][BinZ(z)]->Fill(Corz);
	}
}
void TPCCorrectionMapMaker::WriteHist(){
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
//	OutFile->Write();
void TPCCorrectionMapMaker::MakeSnakeHist(){
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			cout<<Form("Drawing Hist(%d,%d)... please wait",ix,iy)<<endl;
			TString Title = Form("Snake%d_%d_",ix,iy);
			TString DrawCommandX = "xCorVec:cluster_hitpos_z>>"+Title+"x";
			TString DrawCommandY = "yCorVec:cluster_hitpos_z>>"+Title+"y";
			TCut Cut = Form("ntBcOut==1&&nhTpc>15&&abs(cluster_hitpos_x-%f)<%f&&abs(cluster_hitpos_y-%f)<%f",BinPosX(ix)+0*width_x/2,width_x/2,BinPosY(iy)+0*width_x/2,width_y/2);
			SnakeHistX[ix][iy] = new TH2D(Title+"x",Title+"x",snake_nbin_z,-250,250,snake_nbin_x,-5,5);
			SnakeHistY[ix][iy] = new TH2D(Title+"y",Title+"y",snake_nbin_z,-250,250,snake_nbin_x,-5,5);
			DataChain->Draw(DrawCommandX,Cut,"col0");
			DataChain->Draw(DrawCommandY,Cut,"col0");
			cout<<SnakeHistX[ix][iy]->GetEntries()<<endl;
		}
	}
}
void TPCCorrectionMapMaker::LoadCorrectionHist(TString FileName){
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
				peakz[ix][iy][iz]=0;
			}
		}
	}
}
void TPCCorrectionMapMaker::SetParameter(vector<double>px,vector<double>py,int ix,int iy){
	for(int i=0;i<snake_nbin_z;++i){
		peakx[ix][iy][i]=px[i];
		peaky[ix][iy][i]=py[i];
	}
}
void TPCCorrectionMapMaker::WriteSnakeHist(){
	OutFile->cd();
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			TString Title = Form("Snake%d_%d_",ix,iy);
			cout<<Title+"x"<<endl;
			SnakeHistX[ix][iy]->Write();
			cout<<Title+"y"<<endl;
			SnakeHistY[ix][iy]->Write();
			cout<<Title+"Done"<<endl;
		}
	}
}
void TPCCorrectionMapMaker::LoadSnakeHist(TString FileName){
	if(OutFile!=NULL){
		delete OutFile;
	}
	OutFile = new TFile(FileName,"read");
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			TString Title = Form("Snake%d_%d_",ix,iy);
			cout<<Title+"x"<<endl;
			SnakeHistX[ix][iy]=(TH2D*)OutFile->Get(Title+"x");
			SnakeHistY[ix][iy]=(TH2D*)OutFile->Get(Title+"y");
		}
	}
}
TH2D* TPCCorrectionMapMaker::GetSnakeHist(int ix, int iy,bool isX){
	TString Title = Form("Snake%d_%d_",ix,iy);
	if(isX){
		Title+="x";
	}
	else{
		Title+="y";
	}
	cout<<Title<<endl;
	return (TH2D*)OutFile->Get(Title);
}

void TPCCorrectionMapMaker::WriteParam(){
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
void TPCCorrectionMapMaker::ParamFilter(vector<double>&pars,vector<int>ent){

}
