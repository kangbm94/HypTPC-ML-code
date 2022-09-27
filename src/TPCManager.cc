#include "../include/TPCManager.hh"
#ifndef TPCManager_C
#define TPCManager_C
const int max_ntrk = 16;
TPCManager T;
TPCManager gTPCManager;
void TPCManager::LoadChain(TString ChainName ){
	cluster = false;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("raw_padid",&padTpc);
	DataChain->SetBranchAddress("raw_hitpos_y",&dlTpc);
	DataChain->SetBranchAddress("raw_de",&deTpc);
};
void TPCManager::LoadClusterChain(TString ChainName="tpc" ){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);


	DataChain->SetBranchAddress("cluster_hitpos_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_hitpos_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_hitpos_z",&clzTpc);
//	DataChain->SetBranchAddress("cluster_de",&deTpc);
	DataChain->SetBranchAddress("cluster_size",&clsize);
}
void TPCManager::LoadBcOut(){
	DataChain->SetBranchAddress("ntBcOut",&ntBcOut);
	DataChain->SetBranchAddress("x0BcOut",&x0BcOut);
	DataChain->SetBranchAddress("y0BcOut",&y0BcOut);
	DataChain->SetBranchAddress("u0BcOut",&u0BcOut);
	DataChain->SetBranchAddress("v0BcOut",&v0BcOut);
}
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
void TPCManager::FillHist(int itr){
	TVector3 hitv = GetPosition(itr);
	double x = hitv.X();
	double z = hitv.Z();
	PadHist->Fill(z,x);
};
void TPCManager::SetPadContent(int padID,double cont){
	PadHist->SetBinContent(padID,cont);
}





int TPCManager::WhichEvent(){
	const int max_ntrk=20;
	int particle[max_ntrk]={0};
	int ThisEvent=0,npi=0,nk=0,np=0;
	int nh = GetNhitsG4();
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
	int nh = GetNhitsG4();
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
	for(int j=0;j<GetNhitsG4();++j){
		TVector3 vec = GetG4Position(j);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		short x_pix=ToPixel(x_);short y_short = ToShort(y_);short z_pix=ToPixel(z_);
		x[j]=x_pix;y[j]=y_short;z[j]=z_pix;dedx[j]=Getdedxtpc(j);
	}
}
void TPCManager::AssignG4EventD(int* trkid,int* pid, double* x,double* y,double* z,double* dedx){
	for(int j=0;j<max_nh;++j){
		trkid[j]=-999;x[j]=0;y[j]=0;z[j]=0;
	}
	//	cout<<"Initialized"<<endl;
	for(int j=0;j<GetNhitsG4();++j){
		TVector3 vec = GetG4Position(j);
		int trid = Getntrk(j);
		int partid = Getidtpc(j);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		trkid[j]=trid;pid[j]=partid;x[j]=x_;y[j]=y_;z[j]=z_;dedx[j]=Getdedxtpc(j);
	}
}
int TPCManager::AssignRealEvent( double* x,double* y,double* z,double* dedx){
	int nitr = GetNhits(cluster);
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;dedx[j]=0;
	}
	for(int j=0;j<nitr;++j){
		TVector3 vec = GetPosition(j);
		double x_ = vec.X();double y_=vec.Y();double z_ = vec.Z();
		x[j]=x_;y[j]=y_;z[j]=z_;dedx[j]=GetDE(j);
	}
	return nitr;
}
void TPCManager::FillEvent(){
	int nitr = GetNhitsG4();
	double x[max_nh],y[max_nh],z[max_nh],dedx[max_nh];
	int trkid[max_nh],pid[max_nh];
	for(int j=0;j<max_nh;++j){
		x[j]=0;y[j]=0;z[j]=0;
	}
	AssignG4EventD(trkid,pid,x,y,z,dedx);
	for(int j=0;j<nitr;++j){
		FillHist(z[j],x[j]);
	}
}

#endif
