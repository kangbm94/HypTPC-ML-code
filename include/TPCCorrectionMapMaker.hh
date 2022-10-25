#include "TPCManager.hh"
static const int nbin_x=20;
static const int nbin_y=5;
static const int snake_nbin_x = 100;
static const int snake_nbin_z = 40;
static const int nbin_z=snake_nbin_z;
double FillEmptyParam(vector<double>x,int i){
	const int size = x.size();
	double val=0;
	if(size<2){
		val= x[0];
		cout<<"Vector Size 1!"<<endl;
	}
	if(i>=size){
		cout<<"ExeedingVectorSize!"<<endl;
		return val;
	}
	int bi = i-1,fi=i+1;
	int bd=0,fd=0;
	double bp=0,fp=0;
	if(isnan(x[i])){
		while(bi>-1){
			bd++;
			if(!isnan(x[bi])){
				bp=x[bi];
				break;
			}
			bi--;
		}
		if(bi==-1) bd=0;
		while(fi<size){
			fd++;
			if(!isnan(x[fi])){
				fp=x[fi];
				break;
			}
			fi++;
		}
		if(fi==size) fd=0;
		val=double(bp*fd+fp*bd)/(fd+bd);	
	}
	else{
		val=x[i];
	}
	cout<<val<<endl;
	return val;
}
double Wmean(vector<double>x,vector<double>e,int i,int & ent){
	const int size = x.size();
	double max_dif = 1.;
	double val = 0;
	double decay = 1/2;
	if(size<2){
		val= x[0];
		cout<<"Vector Size 1!"<<endl;
	}
	if(i>=size){
		cout<<"ExeedingVectorSize!"<<endl;
	}
	if(i==0){
		if(abs(x[i]-x[i+1])>max_dif){
			val=(e[i]*x[i]+e[i+1]*x[i+1])/(e[i]+e[i+1]);
		}
		else{
			val = x[i];
		}
	}
	else if(i<size-1){
		bool back = abs(x[i]-x[i-1])>max_dif;
		bool fwrd = abs(x[i]-x[i+1])>max_dif;
		cout<<x[i]<<" and "<<x[i+1]<<endl;
		if(fwrd){
			cout<<"Fwrd!"<<endl;
		}
		if(fwrd and back){
			cout<<"Diff!"<<endl;
			val=(e[i-1]*x[i-1]+e[i]*x[i]+e[i+1]*x[i+1])/(e[i-1]+e[i]+e[i+1]);
			cout<<Form("(%f*%f+%f*%f+%f*%f)/(%f+%f+%f)=%f",e[i-1],x[i-1],e[i],x[i],e[i+1],x[i+1],e[i-1],e[i],e[i+1],val)<<endl;
		}
		else if(fwrd){
			val=(e[i-1]*x[i-1]+e[i]*x[i])/(e[i-1]+e[i]);
		}
		else if(back){
			val=(e[i]*x[i]+e[i+1]*x[i+1])/(e[i]+e[i+1]);
		}
		else{
			val=x[i];
		}
	}
	else if(i==size-1){
		if(abs(x[i]-x[i-1])>max_dif){
			val=(e[i]*x[i]+e[i-1]*x[i-1])/(e[i]+e[i-1]);
		}
		else{
			val = x[i];
		}
	}
	if(e[i]==1){
		val*=decay;
	}
	return val;
}

void ExtrapolateEmptyVector(vector<double>bf,vector<double>&af){
	af.clear();
	double cut=-500;
	int size = bf.size();
	int ShortBin=0;
	af.push_back(bf[0]);
	for(int i=1;i<size;++i){
		bool ForceStop=false;
		if(bf[i]>cut){ af.push_back(bf[i]);continue;}
		int j=i;
		while(bf[i]<cut&&!ForceStop){
			if(i>=size){
				i=size;
				af[i-1]=0;
				ForceStop=true;
				break;
			}
			i++;
		}
		//		if(ForceStop) continue;
		double dif = bf[i+1]-bf[j-1];
		int interval = i+2-j;
		for(int k=j;k<i+1;++k){
			af.push_back(bf[j-1]+(k-j)*dif/interval);	
		}
	}
}
class TPCCorrectionMapMaker: public TPCManager{
	protected:
		static const int max_hit = 1000;
		double center_x=0,center_y=5;
		double width_x=5,width_y=10;
		double x1=center_x-width_x*( nbin_x/2);
		double x2=center_x+width_x*( nbin_x/2);
		double y1=center_y-width_y*( nbin_y/2);
		double y2=center_y+width_y*( nbin_y/2);
		TH2D* SnakeHistX[nbin_x][nbin_y];
		TH2D* SnakeHistY[nbin_x][nbin_y];
		TFile* OutFile=NULL;


		double z1=-250,z2=250;


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
		int BinX(double x){
			return Bin(nbin_x,x1,x2,x);
		}
		int BinY(double y){
			return Bin(nbin_y,y1,y2,y);
		}
		int BinZ(double z){
			return Bin(nbin_z,z1,z2,z);
		}
		double BinPosX(double xbin){
			return BinPos(nbin_x,x1,x2,xbin);
		}
		double BinPosY(double ybin){
			return BinPos(nbin_y,y1,y2,ybin);
		}


	public:
		TPCCorrectionMapMaker(){
			for(int ix=0;ix<nbin_x;++ix){
				for(int iy=0;iy<nbin_y;++iy){
					for(int iz=0;iz<nbin_z;++iz){
						TString title = Form("hist%d_%d_%d_",ix,iy,iz);
						hist_x[ix][iy][iz]= new TH1D(title+"x",title+"x",hbin_x,-cx,cx);
						hist_y[ix][iy][iz]= new TH1D(title+"y",title+"y",hbin_y,-cy,cy);
						hist_z[ix][iy][iz]= new TH1D(title+"z",title+"z",hbin_z,-cz,cz);
					}
				}
			}
			cout<<"Hist initialized"<<endl;
		}
		double BinPosZ(double zbin){
			return BinPos(nbin_z,z1+(z2-z1)/(nbin_z)/2,z2-(z2-z1)/(nbin_z)/2,zbin);
		}
		void MakeSnakeFile(TString FileName);
		void WriteHist();
		void WriteSnakeHist();
		void LoadCorrectionHist(TString FileName);
		void LoadSnakeHist(TString FileName);
		void WriteParam();
		void MakeSnakeHist();
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
		TH2D* ScanSnake(int ix,int iy,vector<double>& pars, vector<double>& ents); 


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
		void RestoreEmptyParams();

};
void TPCCorrectionMapMaker::LoadTPCBcOutChain(TString ChainName ){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	cout<<"Chain Loaded: "<<DataChain->GetEntries()<<endl;
	DataChain->SetBranchAddress("cluster_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_z",&clzTpc);
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
			TString Key = Form("Snake%d_%d_",ix,iy);
			TString DrawCommandX = "xCorVec:cluster_hitpos_z>>"+Key+"x";
			
			TString Title = Form("X(%f,%f)Y(%f,%f)",BinPosX(ix)-width_x/2,BinPosX(ix)+width_x/2,BinPosY(iy)-width_y/2,BinPosY(iy)+width_y/2);

			TCut Cut = Form("ntBcOut==1&&nhTpc>15&&abs(cluster_hitpos_x-%f)<%f&&abs(cluster_hitpos_y-%f)<%f",BinPosX(ix)+0*width_x/2,width_x/2,BinPosY(iy),width_y/2);

			TCut Inside = Form("sqrt(cluster_hitpos_x*cluster_hitpos_x+cluster_hitpos_z*cluster_hitpos_z)<250");
			TCut In_frame = Form("abs(abs(cluster_hitpos_x)-abs(cluster_hitpos_z))<%f",frame_width);
			TCut In_Target = Form("abs(cluster_hitpos_z-%f)<%f&&abs(cluster_hitpos_x)<%f",Target_pos,Target_z,Target_x);
			//			TCut Cut = Form("ntBcOut==1&&nhTpc>15&&abs(cluster_x-%f)<%f&&abs(cluster_y-%f)<%f",BinPosX(ix)+0*width_x/2,width_x/2,BinPosY(iy)+0*width_x/2,width_y/2);
			SnakeHistX[ix][iy] = new TH2D(Key+"x",Title+"x",snake_nbin_z,-250,250,snake_nbin_x,-5,5);
			SnakeHistY[ix][iy] = new TH2D(Key+"y",Title+"y",snake_nbin_z,-250,250,snake_nbin_x,-5,5);
			DataChain->Draw(DrawCommandX,Cut and not (In_frame or In_Target) and Inside,"col0");
			//			DataChain->Draw(DrawCommandY,Cut,"col0");
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
	//	cout<<Title<<endl;
	return (TH2D*)SnakeHistX[ix][iy];
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
