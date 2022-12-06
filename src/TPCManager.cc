#include "../include/TPCManager.hh"
#ifndef TPCManager_C
#define TPCManager_C
const int max_ntrk = 16;
TPCManager gTPCManager;
void TPCManager::LoadChain(TString ChainName ){
	cluster = false;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("evnum",&evnum);
	DataChain->SetBranchAddress("nhTpc",&nhittpc);
	DataChain->SetBranchAddress("padTpc",&padTpc);
	DataChain->SetBranchAddress("dlTpc",&dlTpc);
//	DataChain->SetBranchAddress("cdeTpc",&deTpc);
};
void TPCManager::LoadClusterChain(TString ChainName="tpc" ){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("runnum",&runnum);
	DataChain->SetBranchAddress("evnum",&evnum);

	DataChain->SetBranchAddress("cluster_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_z",&clzTpc);
	//	DataChain->SetBranchAddress("cluster_de",&deTpc);
	DataChain->SetBranchAddress("cluster_size",&clsize);
	DataChain->SetBranchAddress("padTpc",&padTpc);
	DataChain->SetBranchAddress("ntTpc",&ntTpc);
	DataChain->SetBranchAddress("helix_cx",&helcxTpc);	
	DataChain->SetBranchAddress("helix_cy",&helcyTpc);	
	DataChain->SetBranchAddress("helix_z0",&helz0Tpc);	
	DataChain->SetBranchAddress("helix_r",&helrTpc);	
	DataChain->SetBranchAddress("helix_dz",&heldzTpc);	
	DataChain->SetBranchAddress("vtx",&vtx);
	DataChain->SetBranchAddress("vty",&vty);
	DataChain->SetBranchAddress("vtz",&vtz);
}
void TPCManager::LoadTPCBcOut(TString filename){
	LoadFile(filename);	
	LoadTPCBcOutChain("tpc");
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
void TPCManager::LoadBcOut(){
	DataChain->SetBranchAddress("ntBcOut",&ntBcOut);
	DataChain->SetBranchAddress("x0BcOut",&x0BcOut);
	DataChain->SetBranchAddress("y0BcOut",&y0BcOut);
	DataChain->SetBranchAddress("u0BcOut",&u0BcOut);
	DataChain->SetBranchAddress("v0BcOut",&v0BcOut);
}
void TPCManager::LoadTPCBcOutChain(TString ChainName ="tpc"){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	cout<<"Chain Loaded: "<<DataChain->GetEntries()<<endl;
	DataChain->SetBranchAddress("ntBcOut",&ntBcOut);
	DataChain->SetBranchAddress("x0BcOut",&x0BcOut);
	DataChain->SetBranchAddress("y0BcOut",&y0BcOut);
	DataChain->SetBranchAddress("u0BcOut",&u0BcOut);
	DataChain->SetBranchAddress("v0BcOut",&v0BcOut);
	DataChain->SetBranchAddress("nhTpc",&nhittpc);
	DataChain->SetBranchAddress("cluster_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_z",&clzTpc);
	DataChain->SetBranchAddress("cluster_de",&deTpc);
	DataChain->SetBranchAddress("cluster_size",&clsize);


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
void TPCManager::FillHist(double z, double x){
	int padId = tpc::findPadID(z,x);
	if(!tpc::Dead(padId)) PadHist->Fill(z,x);
}
void TPCManager::FillHist(int itr){
	TVector3 hitv = GetPosition(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	PadHist->Fill(z,x);
	ZYHist->Fill(z,y);
};
void TPCManager::SetPadContent(int padID,double cont){
	PadHist->SetBinContent(padID,cont);
}
void TPCManager::SetPadContent(int layer,int row,double cont){
	int padID=tpc::GetPadId(layer,row);
	PadHist->SetBinContent(padID,cont);
}
void TPCManager::LoadTPC3D(){
	TPCCanv->cd();
	if(tpcHit3d) delete tpcHit3d;
	int nh;
	if(!cluster){
		nh= GetNhits(0);
	}
	else{
		nh=GetNhits(1);
	}
	tpcHit3d= new TPolyMarker3D(nh,8);    
	for(int i=0;i<nh;++i){
		TVector3 pos = GetPosition(i);
		double x=pos.x(),y=pos.y(),z=pos.z();
		tpcHit3d->SetPoint(i,z,x,y);
	}
//  TView3D *view = (TView3D*) TView::CreateView(1);
	tpcHit3d->Draw("SAME");
	TPCCanv->Modified();
	TPCCanv->Update();
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
void TPCManager::AssignHits(){
	for(int ipad = 0;ipad<max_pad;++ipad){
		int de = PadHist->GetBinContent(ipad+1);
		if(de)m_Hits.push_back(TPCHit(ipad+1,de));
	}
#if 0
	for(auto hit : m_Hits)hit.Show();
#endif
}

bool TPCManager::MakeUpClusters(double Vth=3){
	
	int nh = m_Hits.size();
  if(nh==0) return false;
	vector<int> joined(nh,0);
	int clnum = 0;
	for(int i=0;i<nh;++i){
		if(joined[i])continue;
		vector<TPCHit> Cand;
		auto hit = m_Hits[i];
		if(hit.GetDe()<=Vth)continue;
		int layer = hit.GetLayer();
		hit.SetCluster(clnum);
		Cand.push_back(hit);
		joined[i]++;
		for(int j=0;j<nh;++j){
			auto thit = m_Hits[j];
			if(i==j or joined[j] or layer != thit.GetLayer()) continue;
      int rowID = (int)thit.GetRow();
      for( auto c_hit: Cand){
        int c_rowID = (int)c_hit.GetRow();
        if(tpc::IsClusterable(layer, rowID, c_rowID,2)){
          thit.SetCluster(clnum);
					Cand.push_back(thit);
          joined[j]++;
          break;
        }
      }
		}
//		for(auto cls : Cand) cls.Show();
		TPCCluster cl(Cand);
		m_Clusters.push_back(cl);
		clnum++;
	}
//	cout<<m_Clusters.size()<<endl;
#if 0
	for(auto cl : m_Clusters)	cl.Show();
#endif
	return true;
}

#endif
