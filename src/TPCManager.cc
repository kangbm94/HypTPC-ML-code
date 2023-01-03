#include "../include/TPCManager.hh"
#include "ReconTools.cc"
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
	DataChain->SetBranchAddress("dlTpc",&dlTpc); //	DataChain->SetBranchAddress("cdeTpc",&deTpc);
};
void TPCManager::LoadClusterChain(TString ChainName="tpc" ){
	cluster = true;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("runnum",&runnum);
	DataChain->SetBranchAddress("evnum",&evnum);

	DataChain->SetBranchAddress("nclTpc",&nclTpc);
	DataChain->SetBranchAddress("cluster_de",&cldeTpc);
	DataChain->SetBranchAddress("cluster_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_z",&clzTpc);
	DataChain->SetBranchAddress("nclTpcf",&nclfTpc);
	DataChain->SetBranchAddress("cluster_xf",&clxfTpc);
	DataChain->SetBranchAddress("cluster_yf",&clyfTpc);
	DataChain->SetBranchAddress("cluster_zf",&clzfTpc);
	DataChain->SetBranchAddress("nclTpcb",&nclbTpc);
	DataChain->SetBranchAddress("cluster_xb",&clxbTpc);
	DataChain->SetBranchAddress("cluster_yb",&clybTpc);
	DataChain->SetBranchAddress("cluster_zb",&clzbTpc);
	//	DataChain->SetBranchAddress("cluster_de",&deTpc);
	DataChain->SetBranchAddress("cluster_size",&clsize);
	DataChain->SetBranchAddress("hough_flag",&hough_flag);
	DataChain->SetBranchAddress("padTpc",&padTpc);
	DataChain->SetBranchAddress("ntTpc",&ntTpc);
	DataChain->SetBranchAddress("helix_cx",&helix_cx);	
	DataChain->SetBranchAddress("helix_cy",&helix_cy);	
	DataChain->SetBranchAddress("helix_z0",&helix_z0);	
	DataChain->SetBranchAddress("helix_r",&helix_r);	
	DataChain->SetBranchAddress("helix_dz",&helix_dz);	
	DataChain->SetBranchAddress("chisqr",&chisqr);	
	DataChain->SetBranchAddress("isBeam",&isBeam);
	DataChain->SetBranchAddress("pid",&pid);
	DataChain->SetBranchAddress("charge",&charge);
	DataChain->SetBranchAddress("mom0",&mom0);
	DataChain->SetBranchAddress("dEdx",&dEdx);
	DataChain->SetBranchAddress("vtx",&vtx);
	DataChain->SetBranchAddress("vty",&vty);
	DataChain->SetBranchAddress("vtz",&vtz);
	DataChain->SetBranchAddress("beam_y",&beam_y);
	DataChain->SetBranchAddress("beam_p0",&beam_p0);
	DataChain->SetBranchAddress("beam_p1",&beam_p1);
	DataChain->SetBranchAddress("beam_p2",&beam_p2);
	DataChain->SetBranchAddress("beam_v",&beam_v);
	DataChain->SetBranchAddress("helix_t",&helix_t);
//	DataChain->SetBranchAddress("track_cluster_layer",&track_cluster_layer);
	DataChain->SetBranchAddress("hitlayer",&track_cluster_layer);
	DataChain->SetBranchAddress("track_cluster_x_center",&track_cluster_x_center);
	DataChain->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
	DataChain->SetBranchAddress("track_cluster_z_center",&track_cluster_z_center);
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
void TPCManager::FillHistf(double z, double x){
	int padId = tpc::findPadID(z,x);
	if(!tpc::Dead(padId)) PadHistf->Fill(z,x);
}
void TPCManager::FillHistb(double z, double x){
	int padId = tpc::findPadID(z,x);
	if(!tpc::Dead(padId)) PadHistb->Fill(z,x);
}
void TPCManager::FillHist(int itr){
	TVector3 hitv = GetPosition(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	int bx,by;
	int bin = ZYHist->FindBin(z,y);
//	ZYHist->SetBinContent(bin,x);
	if(hough_flag->at(itr)<9999){
		PadHist->Fill(z,x);
		ZYHist->Fill(z,y);
//		YHist->Fill(y);
	}
};
void TPCManager::FillAntiProtonHist(){
	int nt = charge->size();
	double dz,rad;
	for(int i=0;i<nt;++i){
		if(!IsAntiProton(i)) continue;
		auto posvx = track_cluster_x_center->at(i); 
		auto posvy = track_cluster_y_center->at(i); 
		auto posvz = track_cluster_z_center->at(i); 
		int nh = posvx.size();
		int min_layer = 100; int max_layer = -1;
		double min_tcal,max_tcal;
		cout<<nt<<endl;
		cout<<track_cluster_layer->size()<<endl;
		cout<<helix_t->size()<<endl;
		auto layerv = track_cluster_layer->at(i);
		auto tcalv = helix_t->at(i);
		dz = helix_dz->at(i);
		rad = helix_r->at(i);
		for(int j=0;j<nh;++j){
			double x = posvx.at(j);
			double y = posvy.at(j);
			double z = posvz.at(j);
			int layer = layerv.at(j);
			double tcal = tcalv.at(j);
			int padId = tpc::findPadID(z,x);
			if(!tpc::Dead(padId)){
				PadHist->Fill(z,x);
				ZYHist->Fill(z,y);
			}
			if(min_layer>layer){
				min_layer=layer;
				min_tcal = tcal;
			}
			if(max_layer<layer){
				max_layer=layer;
				max_tcal = tcal;
			}
		}
		cout<<Form("Min( layer,t ) = (%d,%f), Max( layer,t ) = (%d,%f), rad = %f, dz=%f",min_layer,min_tcal,max_layer,max_tcal,rad,dz)<<endl;
	}
}
void TPCManager::FillHistf(int itr){
	TVector3 hitv = GetPositionf(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	if(abs(y)>50)PadHistf->Fill(z,x);
	if(abs(y)>50)ZYHistf->Fill(z,y);
};
void TPCManager::FillHistb(int itr){
	TVector3 hitv = GetPositionb(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	PadHistb->Fill(z,x);
	ZYHistb->Fill(z,y);
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
#if 0 
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
#endif
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
void
TPCManager::ReconEvent(){
	vector<Vertex> verts;
	vector<Track> parts;
	bool ldflg = false,xiflg=false;
	Ld.Clear();Xi.Clear();
	double chi_cut = 50;
	for(int nt1 = 0; nt1<ntTpc;++nt1){
		if(chisqr->at(nt1)>chi_cut) continue; 
		if(isBeam->at(nt1)) continue; 
		int nh = helix_cx->size();
		double hcx = helix_cx->at(nt1);
		double hcy = helix_cy->at(nt1);
		double hz0 = helix_z0->at(nt1);
		double hr = helix_r->at(nt1);
		double hdz = helix_dz->at(nt1);
		double par1[5] = {hcx,hcy,hz0,hr,hdz};
		int id1 = pid->at(nt1);
		double q1 = charge->at(nt1);
		parts.push_back(Track(id1,q1,par1,nt1));
	}
	int np = parts.size();
	if(np<1) return;
	for(int nt1=0;nt1<np;++nt1){
		Vertex f(parts[nt1]);
		for(int nt2=nt1+1;nt2<np;++nt2){
			f.AddTrack(parts[nt2]);	
		}
		//if(f.NTrack()>1) verts.push_back(f);
		verts.push_back(f);
	}
	vector<Recon>LdCand;
	int nvt = verts.size();
	for(auto vt: verts){
		vt.SearchLdCombination();
		auto Ldc = vt.GetLd();
		LdCand.push_back(Ldc);
	}
	int nld= LdCand.size();
	double comp = 9999;

	for(auto m:LdCand){
		if( abs(mL-m.Mass())<comp) {comp=abs(mL-m.Mass());Ld=m;}
	}
	comp = 9999;
	VertexLH V(Ld);

	for(auto p : parts){
		V.AddTrack(p);
	}
	ldflg=Ld.Exist();
	if(ldflg)V.SearchXiCombination();	
	Xi = V.GetXi();
	xiflg=Xi.Exist();
}
void TPCManager::DrawHelix(int it){
		HelixTrack[it]->Draw("psame");
}
void TPCManager::DrawHelix(){
	cout<<"DrawingHelix..."<<endl;
	for(int it = 0; it< ntTpc;++it)DrawHelix(it);
}

void TPCManager::DrawHelixZY(int it){
	for(int ip=0;ip<npts;ip++){ 
		HelixTrackZY[it].at(ip)->Draw("same");
	} 
}
void TPCManager::DrawHelixZY(){
	cout<<"DrawingHelixZY..."<<endl;
	for(int it = 0; it< ntTpc;++it)DrawHelixZY(it);
}

#endif
