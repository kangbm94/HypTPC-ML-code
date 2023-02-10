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
	cout<<"LoadingFile..."<<endl;
	DataChain	= (TChain*) DataFile->Get(ChainName);
	DataChain->SetBranchAddress("trigflag",&trigflag);
	DataChain->SetBranchAddress("runnum",&runnum);
	DataChain->SetBranchAddress("evnum",&evnum);

	DataChain->SetBranchAddress("nclTpc",&nclTpc);
	DataChain->SetBranchAddress("cluster_de",&cldeTpc);
	DataChain->SetBranchAddress("cluster_x",&clxTpc);
	DataChain->SetBranchAddress("cluster_y",&clyTpc);
	DataChain->SetBranchAddress("cluster_z",&clzTpc);
	DataChain->SetBranchAddress("resolution_x",&resolution_x);
	DataChain->SetBranchAddress("resolution_y",&resolution_y);
	DataChain->SetBranchAddress("resolution_z",&resolution_z);
	DataChain->SetBranchAddress("resolution",&resolution);
	DataChain->SetBranchAddress("hough_flag",&hough_flag);
	DataChain->SetBranchAddress("hough_dist",&hough_dist);
	
	DataChain->SetBranchAddress("padTpc",&padTpc);

	DataChain->SetBranchAddress("ntAcc",&ntAcc);
	DataChain->SetBranchAddress("accidental_cx",&accidental_cx);	
	DataChain->SetBranchAddress("accidental_cy",&accidental_cy);	
	DataChain->SetBranchAddress("accidental_z0",&accidental_z0);	
	DataChain->SetBranchAddress("accidental_r",&accidental_r);	
	DataChain->SetBranchAddress("accidental_dz",&accidental_dz);	
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

	DataChain->SetBranchAddress("track_cluster_x_center",&track_cluster_x_center);
	DataChain->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
	DataChain->SetBranchAddress("track_cluster_z_center",&track_cluster_z_center);
	/*
		 DataChain->SetBranchAddress("beam_y",&beam_y);
		 DataChain->SetBranchAddress("beam_p0",&beam_p0);
		 DataChain->SetBranchAddress("beam_p1",&beam_p1);
		 DataChain->SetBranchAddress("beam_p2",&beam_p2);
		 DataChain->SetBranchAddress("beam_v",&beam_v);
		 DataChain->SetBranchAddress("helix_t",&helix_t);
		 */
	//	DataChain->SetBranchAddress("track_cluster_layer",&track_cluster_layer);
	DataChain->SetBranchAddress("hitlayer",&track_cluster_layer);
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

void TPCManager::FillAccHists(){
	int nh = GetNhits();
	for(int ih = 0; ih < nh; ++ih){
		int hf = hough_flag->at(ih) - 1000;
		if(hf > -1 and hf < 20){
			double x = clxTpc->at(ih);
			double y = clyTpc->at(ih);
			double z = clzTpc->at(ih);
			ZYHistsAcc[hf]->Fill(z,y);
			CirHistsAcc[hf]->Fill(z,x);
		}
	}
}
TH2D* TPCManager::GetZYHistAcc(int i){
	return ZYHistsAcc[i];
}
TH2D* TPCManager::GetCirHistAcc(int i){
	return CirHistsAcc[i];
}

void TPCManager::ClearHistogram(){
	PadHist->Reset("");
	ZYHist->Reset("");
	YHist->Reset("");
	for(int i = 0; i<20;++i){
		ZYHistsAcc[i] ->Reset();
		CirHistsAcc[i]->Reset();
	}
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
	double res = resolution->at(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	int bx,by;
	int bin = ZYHist->FindBin(z,y);
	int padId = tpc::findPadID(z,x);
	PadHist->SetBinContent(padId,res);
	//	ZYHist->SetBinContent(bin,x);
	//	if(hough_flag->at(itr)<9999){
//	PadHist->Fill(z,x);
	ZYHist->Fill(z,y);
	if(hough_flag->at(itr)<990)
	cout<<Form("Pos(%f,%f,%f)",x,y,z)<<endl;
	//		YHist->Fill(y);
	//	}
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
	for(int it = 0; it< ntTpc;++it)DrawHelix(it);
}

void TPCManager::DrawAccidental(int it){
	AccidentalTrack[it]->Draw("psame");
}
void TPCManager::DrawAccidental(){
	for(int it = 0; it< ntAcc;++it)DrawAccidental(it);
}

void TPCManager::DrawHelixZY(int it){
	for(int ip=0;ip<npts;ip++){ 
		HelixTrackZY[it].at(ip)->Draw("same");
	} 
}
void TPCManager::DrawHelixZY(){
	for(int it = 0; it< ntTpc;++it)DrawHelixZY(it);
}

void TPCManager::DrawAccidentalZY(int it){
	for(int ip=0;ip<anpts;ip++){ 
		AccidentalTrackZY[it].at(ip)->Draw("same");
	} 
}
void TPCManager::DrawAccidentalZY(){
	for(int it = 0; it< ntAcc;++it)DrawAccidentalZY(it);
}
void TPCManager::InitializeAccidental(){
	AccidentalTrackZY.clear();
	AccidentalTrackZY.resize(ntAcc);
	for(int it = 0; it< ntAcc;++it){
		double cx = accidental_cx->at(it);
		double cy = accidental_cy->at(it);
		double z0 = accidental_z0->at(it);
		double r = accidental_r->at(it);
		double dz = accidental_dz->at(it);
		double pars[5]={cx,cy,z0,r,dz};
		double theta1=70;
		double theta2=105;
		AccidentalTrack[it] = new TEllipse(cy+ZTarget,-cx,r,r,theta1,theta2);
		AccidentalTrack[it] ->SetNoEdges();
		AccidentalTrack[it] ->SetFillStyle(0);
		AccidentalTrack[it] ->SetLineColor(it+1);
		AccidentalTrack[it] ->SetLineWidth(2);
		AccidentalTrack[it] ->SetLineStyle(10);
		TVector3 pos(0,z0,0);
		double t = GetTcal(pars,pos);
		double th1 =   acos(-1)*(1-20./180); 
		double th2 =   acos(-1)*(1+15./180); 
		double dt = (th2-th1)/(anpts);
		for(int ip=0;ip<anpts;ip++){
			double t1 = th1+ip*dt,t2 = t1+dt;
			double y1 = r*dz*t1+z0,y2 = r*dz*t2+z0;
			double z1 = r*sin(t1)+cy+ZTarget,z2 = r*sin(t2)+cy+ZTarget;
			AccidentalTrackZY[it].push_back(new TLine(z1,y1,z2,y2));
			AccidentalTrackZY[it].at(ip)->SetLineColor(it+1);
			AccidentalTrackZY[it].at(ip)->SetLineWidth(2);
			AccidentalTrackZY[it].at(ip)->SetLineStyle(10);
		}
	}
}
vector<double>* TPCManager::GetAccidentalDist(){
	int nh = GetNhits(1);
	vector<double>*dists = new vector<double>;
	delete dists;
	dists = new vector<double>;
	double t_ = -9999;
	for(int ih = 0; ih< nh; ++ih){
		double min_dist = 9999;
		for(int iac = 0; iac < ntAcc;++iac){
			TVector3 pos(
				clxTpc->at(ih),
				clyTpc->at(ih),
				clzTpc->at(ih));
			double par[8] = {
				accidental_cx->at(iac),
				accidental_cy->at(iac),
				accidental_z0->at(iac),
				accidental_r->at(iac),
				accidental_dz->at(iac),
				pos.X(),pos.Y(),pos.Z()
			};
			double t = GetTcal(par,pos);
			TVector3 Hpos = HelixPos(par,t);
			double dist = (Hpos-pos).Mag();
			if(dist<min_dist){
				min_dist = dist;
				t_=t;
			}
		}
		if(min_dist<50)cout<<"T = "<<t_<<endl;
		dists->push_back(min_dist);
	}
	return dists;
}
void TPCManager::InitializeHelix(){
	//			cout<<"InitializeHelix: "<<endl;
	//			HelixTrack.clear();		

	HelixTrackZY.clear();		
	HelixTrackZY.resize(ntTpc);		
	for(int it = 0; it< ntTpc;++it){
		double cx = helix_cx->at(it);
		double cy = helix_cy->at(it);
		double z0 = helix_z0->at(it);
		double r = helix_r->at(it);
		double dz = helix_dz->at(it);
		//			cout<<Form("Params = (%f,%f,%f,%f,%f)",cx,cy,r,z0,dz)<<endl;
		//				TString title = Form("Helix%d",it);
		//				if(r>4000) continue;
		double pars[5]={cx,cy,z0,r,dz};
		double t_min = 100;
		double t_max = -100;
		auto xcl = track_cluster_x_center->at(it);
		auto ycl = track_cluster_y_center->at(it);
		auto zcl = track_cluster_z_center->at(it);
		vector<double>tvec;
		for(int ih=0;ih<xcl.size();++ih){
			//					cout<<"track : "<< it<<" hit: "<<ih<<endl;
			double x = xcl.at(ih);
			double y = ycl.at(ih);
			double z = zcl.at(ih);
			TVector3 pos(x,y,z);
			double t = GetTcal(pars,pos);
			if(t<t_min) t_min=t;
			if(t>t_max) t_max=t;
			double hx = cos(t),hy=sin(t);
			double theta = atan2(-hx,hy);
			//					theta=fmod(theta*180./acos(-1)+360.*100,360.);
			theta=fmod(theta*180./acos(-1),360.);
			tvec.push_back(theta);
		}
		double theta_med = TMath::Median(tvec.size(),tvec.data());
		//				cout<<t_min<<","<<t_max<<endl;
		double hx_min = cos(t_min),hy_min=sin(t_min);
		double hx_max = cos(t_max),hy_max=sin(t_max);
		double tc_min = atan2(-hx_min,hy_min),tc_max=atan2(-hx_max,hy_max);
		//double theta1=fmod(tc_min*180./acos(-1)+360.*100.,360.),theta2=fmod(tc_max*180./acos(-1)+360.*100.,360.);
		double theta1=fmod(tc_min*180./acos(-1),360.),theta2=fmod(tc_max*180./acos(-1),360.);
		double th1_temp = theta1-theta_med;// = fmod(theta1-theta_med,360.);
		double th2_temp = theta2-theta_med;// = fmod(theta2-theta_med,360.);
																			 //				cout<<Form("t1,tm,t2 = (%f,%f,%f)",theta1,theta_med,theta2)<<endl; 

		if( sin(th1_temp*acos(-1)/180.)>0  and sin(th2_temp*acos(-1)/180.)<0 ){
			theta2=theta2-360;
		}

		//				cout<<Form("t1,tm,t2 = (%f,%f,%f)",theta1,theta_med,theta2)<<endl; 
		//				theta1=360.-theta1;
		//				theta2=360.-theta2;

		/*
			 if(theta1>theta2){
			 double dum = theta2;
			 theta2=theta1;
			 theta1=dum;
			 }*/
		//			cout<<theta1<<","<<theta2<<endl;
		HelixTrack[it] = new TEllipse(cy+ZTarget,-cx,r,r,theta1,theta2);
		HelixTrack[it]->SetNoEdges();
		//				HelixTrack[it] = new TEllipse(cy+ZTarget,-cx,r,r,0.,360.);
		//				HelixTrack[it] = new TEllipse(cy+ZTarget,-cx,r,r);
		HelixTrack[it]-> SetLineColor(kRed);
		HelixTrack[it]-> SetFillStyle(0);
		HelixTrack[it]-> SetLineColor(it+1);
		HelixTrack[it]-> SetLineWidth(2);
		double dt = (t_max-t_min)/npts;
		for(int ip=0;ip<npts;ip++){
			double t1 = t_min+dt*ip,t2 = t1+dt;
			double y1 = r*dz*t1+z0,y2 = r*dz*t2+z0;
			double z1 = r*sin(t1)+cy+ZTarget,z2 = r*sin(t2)+cy+ZTarget;
			HelixTrackZY[it].push_back(new TLine(z1,y1,z2,y2));
			HelixTrackZY[it].at(ip)->SetLineColor(it+1);
			HelixTrackZY[it].at(ip)->SetLineWidth(2);
		}
	}
	cout<<"HelixInit!"<<endl;
}
#endif
