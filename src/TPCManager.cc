#include "../include/TPCManager.hh"
#include "ReconTools.cc"
#ifndef TPCManager_C
#define TPCManager_C
namespace{
	bool IsAccidental(int hf){
		if(399 < hf and hf < 500) return true;
		else return false;
	}
	bool IsKurama(int hf){
		if(299 < hf and hf < 400) return true;
		else return false;
	}
	bool IsK18(int hf){
		if(199 < hf and hf < 300) return true;
		else return false;
	}

}
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
	DataChain->SetBranchAddress("trigflag",&trigflag);
	DataChain->SetBranchAddress("runnum",&runnum);
	DataChain->SetBranchAddress("evnum",&evnum);

	DataChain->SetBranchAddress("nclTpc",&nclTpc);
	DataChain->SetBranchAddress("cluster_de",&cldeTpc);
	DataChain->SetBranchAddress("cluster_x",&cluster_x);
	DataChain->SetBranchAddress("cluster_y",&cluster_y);
	DataChain->SetBranchAddress("cluster_z",&cluster_z);
	DataChain->SetBranchAddress("hough_flag",&hough_flag);
	DataChain->SetBranchAddress("hough_dist",&hough_dist);
	

	DataChain->SetBranchAddress("ntTpc",&ntTpc);
	DataChain->SetBranchAddress("helix_cx",&helix_cx);	
	DataChain->SetBranchAddress("helix_cy",&helix_cy);	
	DataChain->SetBranchAddress("helix_z0",&helix_z0);	
	DataChain->SetBranchAddress("helix_r",&helix_r);	
	DataChain->SetBranchAddress("helix_dz",&helix_dz);	
	DataChain->SetBranchAddress("helix_flag",&helix_flag);	
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
	DataChain->SetBranchAddress("helix_t",&helix_t);
	DataChain->SetBranchAddress("hitlayer",&track_cluster_layer);

	DataChain->SetBranchAddress("ntKurama",&ntKurama);
	DataChain->SetBranchAddress("qKurama",qKurama);
	DataChain->SetBranchAddress("m2",m2);
	DataChain->SetBranchAddress("pKurama",pKurama);

	DataChain->SetBranchAddress("nKK",&nKK);
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
	DataChain->SetBranchAddress("cluster_x",&cluster_x);
	DataChain->SetBranchAddress("cluster_y",&cluster_y);
	DataChain->SetBranchAddress("cluster_z",&cluster_z);
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
			double x = cluster_x->at(ih);
			double y = cluster_y->at(ih);
			double z = cluster_z->at(ih);
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
//	double res = resolution->at(itr);
	double x = hitv.X();
	double y = hitv.Y();
	double z = hitv.Z();
	int bx,by;
	int bin = ZYHist->FindBin(z,y);
	int padId = tpc::findPadID(z,x);
	PadHist->SetBinContent(padId,1);
	//	ZYHist->SetBinContent(bin,x);
	//	if(hough_flag->at(itr)<9999){
//	PadHist->Fill(z,x);
	if(hough_flag->at(itr) < 400 or hough_flag->at(itr)>499)ZYHist->Fill(z,y);
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
//	if(tpcHit3d) delete tpcHit3d;
	int nh;
	if(!cluster){
		nh= GetNhits(0);
	}
	else{
		nh=GetNhits(1);
	}
	tpcHit3d.clear();
	tpcHit3d.resize(nh);
	for(int i=0;i<nh;++i){
		tpcHit3d.at(i) = new TPolyMarker3D();
		int hf = hough_flag->at(i);
		TVector3 pos = GetPosition(i);
		double x=pos.x(),y=pos.y(),z=pos.z();
		tpcHit3d.at(i)->SetPoint(0,z,x,y);
		tpcHit3d.at(i)->SetMarkerStyle(8);	
		if(399 < hf and hf < 500){
			tpcHit3d.at(i)->SetMarkerColor(hf - 399);	
			tpcHit3d.at(i)->SetMarkerStyle(29);	
		}
		if(hf == 200){
			tpcHit3d.at(i)->SetMarkerColor(46);	
			tpcHit3d.at(i)->SetMarkerStyle(101);	
		}
		if(hf == 300){
			tpcHit3d.at(i)->SetMarkerColor(40);	
			tpcHit3d.at(i)->SetMarkerStyle(103);	
		}
		tpcHit3d.at(i)->Draw("SAME");
	}
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
void TPCManager::DrawHelix(int it){
	HelixTrack[it]->Draw("psame");
}
void TPCManager::DrawHelix(){
	for(int it = 0; it< HelixTrack.size();++it){
		if(helix_flag->at(it)> 399 and helix_flag->at(it)<500){
			DrawHelix(it);
//			cout<<it<<" , "<<helix_flag->at(it)<<endl;
		}
	}
}
void TPCManager::DrawHelix3D(int it){
	HelixTrack3D[it]->Draw("same");
}
void TPCManager::DrawHelix3D(){
	for(int it=0;it<HelixTrack3D.size();++it){
		HelixTrack3D[it]->Draw("same");
	}
}
void TPCManager::DrawAccidental3D(){
	for(int it=0;it<AccidentalTrack3D.size();++it){
		AccidentalTrack3D[it]->Draw("same");
	}
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
	for(int it = 0; it< HelixTrack.size();++it)DrawHelixZY(it);
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
		double cx = helix_cx->at(it);
		double cy = helix_cy->at(it);
		double z0 = helix_z0->at(it);
		double r = helix_r->at(it);
		double dz = helix_dz->at(it);
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
				cluster_x->at(ih),
				cluster_y->at(ih),
				cluster_z->at(ih));
			double par[8] = {
				helix_cx->at(iac),
				helix_cy->at(iac),
				helix_z0->at(iac),
				helix_r->at(iac),
				helix_dz->at(iac),
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
	HelixTrack.clear();		
	HelixTrackZY.clear();		
	AccidentalTrack.clear();		
	AccidentalTrackZY.clear();		
	int TrackNo = 0;
	for(int it = 0; it< ntTpc;++it){
		int hf = helix_flag ->at(it);
		double cx = helix_cx->at(it);
		double cy = helix_cy->at(it);
		double z0 = helix_z0->at(it);
		double r = helix_r->at(it);
		double dz = helix_dz->at(it);
		double pars[5]={cx,cy,z0,r,dz};
		double t_min = 100;
		double t_max = -100;
		auto xcl = track_cluster_x_center->at(it);
		auto ycl = track_cluster_y_center->at(it);
		auto zcl = track_cluster_z_center->at(it);
		vector<double>tvec;
		for(int ih=0;ih<xcl.size();++ih){
			double x = xcl.at(ih);
			double y = ycl.at(ih);
			double z = zcl.at(ih);
			TVector3 pos(x,y,z);
			double t = GetTcal(pars,pos);
			if(t<t_min) t_min=t;
			if(t>t_max) t_max=t;
			double hx = cos(t),hy=sin(t);
			double theta = atan2(-hx,hy);
			theta=fmod(theta*180./acos(-1),360.);
			tvec.push_back(theta);
		}
		double theta_med = TMath::Median(tvec.size(),tvec.data());
		double hx_min = cos(t_min),hy_min=sin(t_min);
		double hx_max = cos(t_max),hy_max=sin(t_max);
		double tc_min = atan2(-hx_min,hy_min),tc_max=atan2(-hx_max,hy_max);
		double theta1=fmod(tc_min*180./acos(-1),360.),theta2=fmod(tc_max*180./acos(-1),360.);
		double th1_temp = theta1-theta_med;// = fmod(theta1-theta_med,360.);
		double th2_temp = theta2-theta_med;// = fmod(theta2-theta_med,360.);

		if( sin(th1_temp*acos(-1)/180.)>0  and sin(th2_temp*acos(-1)/180.)<0 ){
			theta2=theta2-360;
		}

		auto Track= new TEllipse(cy+ZTarget,-cx,r,r,theta1,theta2);
		Track->SetNoEdges();
		Track-> SetLineColor(kRed);
		Track-> SetFillStyle(0);
		Track-> SetLineColor(it+1);
		Track-> SetLineWidth(2);
		vector<TLine*> TrackZY;
		double dt = (t_max-t_min)/npts;
		for(int ip=0;ip<npts;ip++){
			double t1 = t_min+dt*ip,t2 = t1+dt;
			double y1 = r*dz*t1+z0,y2 = r*dz*t2+z0;
			double z1 = r*sin(t1)+cy+ZTarget,z2 = r*sin(t2)+cy+ZTarget;
			TrackZY.push_back(new TLine(z1,y1,z2,y2));
		}
		if(hf> 399 and hf < 500){
			TrackNo++;
			Track-> SetLineColor(1);
			Track-> SetLineWidth(2);
			for(auto t: TrackZY){
				t->SetLineColor(1);
				t->SetLineStyle(kDashed);
			}
			AccidentalTrackZY.push_back(TrackZY);
			AccidentalTrack.push_back(Track);
		}
		else{
			Track-> SetLineColor(TrackNo);
			Track-> SetLineWidth(2);
			for(auto t: TrackZY){
				t->SetLineColor(TrackNo);
			}
			HelixTrackZY.push_back(TrackZY);
			HelixTrack.push_back(Track);
		}
	}
}
void TPCManager::LoadAccidental3D(){
	AccidentalTrack3D.clear();
	for(int it = 0; it< ntTpc;++it){
		int track_flag = helix_flag->at(it);
		if(track_flag < 400 or track_flag > 499) continue;
//		cout<<"flag : "<<track_flag<<endl;
		
		double cx = helix_cx->at(it);
		double cy = helix_cy->at(it);
		double z0 = helix_z0->at(it);
		double r = helix_r->at(it);
		double dz = helix_dz->at(it);
		
		double pars[5]={cx,cy,z0,r,dz};
		double t_min = 100;
		double t_max = -100;
		vector<double>tvec;
		vector<TVector3>tpos;
		for(int ih=0;ih<helix_t->at(it).size();++ih){
			double t = helix_t->at(it).at(ih);
			if(t<t_min) t_min=t;
			if(t>t_max) t_max=t;
			tvec.push_back(t);
		}
		sort(tvec.begin(),tvec.end());
		double dt = (t_max-t_min)/npts;
		auto Track = new TPolyLine3D(npts);
		Track ->SetLineColor(2);
		Track ->SetLineWidth(3);
		Track ->SetLineStyle(kDashed);
		for(int ip=0;ip<npts;ip++){
			double t1 = t_min+dt*ip,t2 = t1+dt;
			double x1 = -(r*cos(t1)+cx); 
			double y1 = r*dz*t1+z0;
			double z1 = r*sin(t1)+cy+ZTarget;
			Track->SetPoint(ip,z1,x1,y1);	
		}
//		cout<<"Accidental Track " << it << " radius "<<r<<" rdz = "<<r*dz<<" z0 = "<<z0<<" T = ("<<t_min<<" , "<<t_max<<" )"<<endl; 
		for(auto t:tvec){
//			cout<<"Accidental Track " << it <<" "<<HelixPos(pars,t).Y()<< " T = "<<t<<endl; 
		}
		AccidentalTrack3D.push_back(Track);
	}//ntTpc
}
void TPCManager::LoadHelix3D(){
	HelixTrack3D.clear();
	HelixTrackID.clear();
	HelixTrackMom.clear();
	const double Const = 0.299792458; // =c/10^9
	const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
	int TrackNo=0;
	for(int it = 0; it< ntTpc;++it){
		int track_flag = helix_flag->at(it);
//		cout<<"flag : "<<track_flag<<endl;
		if(track_flag > 399 and track_flag < 500){
			HelixTrackID.push_back(-1);
			HelixTrackMom.push_back(1800);
			continue;
		}
		else{
			HelixTrackID.push_back(TrackNo);
			TrackNo++;	
		}
		double cx = helix_cx->at(it);
		double cy = helix_cy->at(it);
		double z0 = helix_z0->at(it);
		double r = helix_r->at(it);
		double dz = helix_dz->at(it);
		
		double pt = abs(r)* Const*dMagneticField;
		double pz = pt * dz;
		double mom = sqrt(pt*pt+pz*pz);
		HelixTrackMom.push_back(mom);
		double pars[5]={cx,cy,z0,r,dz};
		double t_min = 100;
		double t_max = -100;
		vector<double>tvec;
		vector<double>tvec2;
		vector<TVector3>tpos;
		for(int ih=0;ih<helix_t->at(it).size();++ih){
			double t = helix_t->at(it).at(ih);
			tvec.push_back(t);
			if(t>0){
				tvec2.push_back(t);
			}
			else{
				tvec2.push_back(t+2*acos(-1));
			}
		}
		sort(tvec.begin(),tvec.end());
		sort(tvec2.begin(),tvec2.end());
		double tgap1=-999;
		double tgap2=-999;
		int nt= tvec.size();
		for(int it=0;it<nt-1;++it){
			double tg = tvec.at(it+1)-tvec.at(it);
			if(tgap1 < tg) tgap1 = tg;
			double tg2 = tvec2.at(it+1)-tvec2.at(it);
			if(tgap2 < tg) tgap2 = tg2;
		}
		if(tgap1 >tgap2 and 0){
			t_min = tvec2.at(0);
			t_max = tvec2.at(nt-1);
		}
		else{
			t_min = tvec.at(0);
			t_max = tvec.at(nt-1);
		}
		double dt = (t_max-t_min)/npts;
//		cout<<theta_med<<endl;
//		t_max = 1.1*acos(-1);
//		t_min = 0.9*acos(-1);
		auto Track = new TPolyLine3D(npts);
		Track ->SetLineColor(TrackNo);
		Track ->SetLineWidth(3);
		if(track_flag == 200){
			Track ->SetLineColor(46);
			Track ->SetLineStyle(kDashDotted);
		}
		if(track_flag == 300){
			Track ->SetLineColor(40);
			Track ->SetLineStyle(kDashDotted);
		}
		for(int ip=0;ip<npts;ip++){
			double t1 = t_min+dt*ip;
			double x1 = -(r*cos(t1)+cx); 
			double y1 = r*dz*t1+z0;
//			if(abs(y1)>250) cout<<Form("y>250! t1 = %f, z0 = %f",t1,z0)<<endl;
			double z1 = r*sin(t1)+cy+ZTarget;
			Track->SetPoint(ip,z1,x1,y1);	
		}
//		cout<<"Track "<<it <<"Momentum "<<mom<<"MeV/c"<<endl;
//		cout<<"Track " << it << " radius "<<r<<" rdz = "<<r*dz<<" z0 = "<<z0<<" T = ("<<t_min<<" , "<<t_max<<" )"<<endl; 
		for(auto t:tvec){
//			cout<<"Track " << it <<" "<<HelixPos(pars,t).Y()<< " T = "<<t<<endl; 
		}
		HelixTrack3D.push_back(Track);
	}//ntTpc
}

void
TPCManager::DoZYHough(){
	cout<<"ZYReset"<<endl;
	hist_ZY->Reset();
	ZYLine.clear();
	cout<<"ZYReset"<<endl;
	int nh = GetNhits(1);
	cout<<"nh = "<<nh<<endl;
	vector<bool>countflag(nh,false);
	for(int it = 0; it < 10; ++it){
		for(int ih = 0; ih < nh; ++ih){
			if(countflag.at(ih)) continue;
			if(hough_flag->at(ih)>0) continue;
			auto pos = GetPosition(ih);
			double x = pos.X(),y=pos.Y(),z=pos.Z();
			for(int ti = 0; ti<ZYtheta_ndiv;++ti){
				double theta = hist_ZY->GetXaxis()->GetBinCenter(ti+1);
				double rho = cos(theta)*z + sin(theta)*y;
				hist_ZY->Fill(theta,rho);
			}//ti

		}//ih
		int maxbin = hist_ZY->GetMaximumBin();
		int mx,my,mz;
		hist_ZY->GetBinXYZ(maxbin,mx,my,mz);
		double mtheta = hist_ZY->GetXaxis()->GetBinCenter(mx);
		double mrho = hist_ZY->GetYaxis()->GetBinCenter(my);
		double y0 = mrho/sin(mtheta);
		double v = -cos(mtheta)/sin(mtheta);
		vector<int> index;
		int nc = 0;
		for(int ih = 0; ih < nh; ++ih){
			if(countflag.at(ih)) continue;
			if(hough_flag->at(ih)>0) continue;
			auto pos = GetPosition(ih);
			double x = pos.X(),y=pos.Y(),z=pos.Z();
			double Ydist = sqrt( (y0 + v*z  - y)/(v*v+1));
			if(Ydist < 20){
				index.push_back(ih);
				nc++;
			}
		}
		if(nc> 8){
			for(auto ih : index){
				countflag.at(ih)=true;
			}
			cout<<Form("Track %d: hits = %d , y = %f +  %f z ",it,nc,y0,v)<<endl;
			cout<<Form("Params: (rho , theta ) = ( %f , %f) ",mrho,mtheta)<<endl;
			ZYLine.push_back(new TLine(-250,y0-250*v,250,y0+250*v));		
		}
	}//it
};
void 
TPCManager::DrawZYHough(){
	int nZY = ZYLine.size();
	for(int i=0;i<nZY;++i){
		ZYLine.at(i)->Draw("same");
		ZYLine.at(i)->SetLineColor(i+1);
		ZYLine.at(i)->SetLineWidth(1);
	}
}


void TPCManager::DrawVertex3D(){
			const double Const = 0.299792458; // =c/10^9
			const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
			TVector3 LV,XV,XVCor;
			bool ldflg = Ld.Exist(),xiflg = Xi.Exist();
			double z1=0,z2=0,z3=0,z4=0;
			double x1=0,x2=0,x3=0,x4=0;
			double y1=0,y2=0,y3=0,y4=0;
			Vertex3d.clear();	
			VertexTrack3D.clear();	
			if(ldflg) {
				LV = Ld.Vertex();
				z1 = LV.Z();
				x1 = LV.X();
				y1 = LV.Y();
				auto Dir = Ld.Momentum();
				auto LVM = Ld.Momentum();
				Dir = Dir* (1./Dir.Mag());
				z2 = z1-Dir.Z()*150;
				x2 = x1-Dir.X()*150;
				y2 = y1-Dir.Y()*150;
				cout<<Form("Lambda Vertex (%f,%f,%f) Mass: %f",LV.X(),LV.Y(),LV.Z(),Ld.Mass())<<endl;
				cout<<Form("Lambda Mom (%f,%f,%f) mag: %f",LVM.X(),LVM.Y(),LVM.Z(),LVM.Mag())<<endl;
				auto LdVert = new TPolyMarker3D();
				LdVert->SetPoint(0,z1,x1,y1);
				LdVert->SetMarkerColor(kMagenta);
				LdVert->SetMarkerStyle(39);
				LdVert->SetMarkerSize(2);
				Vertex3d.push_back(LdVert);
				auto LdTrack = new TPolyLine3D(2);
				LdTrack->SetPoint(0,z1,x1,y1);
				LdTrack->SetPoint(1,z2,x2,y2);
				LdTrack->SetLineColor(kMagenta);
				LdTrack->SetLineWidth(2);
//				VertexTrack3D.push_back(LdTrack);
				auto Pi = Ld.GetDaughter(1);
				auto PiP = Pi.Vect()*1000;
//				cout<<"LdPiID = "<<LdPiID<<endl;
//				cout<<Form("LdPiMom = (%f,%f,%f), mag = %f MeV/c, cd = %f mm",PiP.X(),PiP.Y(),PiP.Z(),HelixTrackMom.at(LdPiID),Ld.GetCD())<<endl;
//				cout<<"LdPiTrackID = "<<HelixTrackID.at(LdPiID)<<endl;
				
				auto Proton = Ld.GetDaughter(0);
				auto PP = Proton.Vect()*1000;
//				cout<<"LdProtonID = "<<LdProtonID<<endl;
//				cout<<Form("LdProtonMom = (%f,%f,%f), mag = %f MeV/c",PP.X(),PP.Y(),PP.Z(),HelixTrackMom.at(LdProtonID))<<endl;
//				cout<<"LdProtonTrackID = "<<HelixTrackID.at(LdProtonID)<<endl;
				HelixTrack3D.at(HelixTrackID.at(LdPiID))->SetLineColor(kMagenta);
				HelixTrack3D.at(HelixTrackID.at(LdPiID))->SetLineWidth(3);
				HelixTrack3D.at(HelixTrackID.at(LdPiID))->SetLineStyle(2);
				HelixTrack3D.at(HelixTrackID.at(LdProtonID))->SetLineColor(kMagenta);
				HelixTrack3D.at(HelixTrackID.at(LdProtonID))->SetLineWidth(3);
				HelixTrack3D.at(HelixTrackID.at(LdProtonID))->SetLineStyle(6);
			}
			if(xiflg){
				XV = Xi.Vertex();
//				XVCor = XiCor.Vertex();
				z3 = XV.Z();
				x3 = XV.X();
				y3 = XV.Y();
				auto XiPi = Xi.GetDaughter(1);
				auto XiPiP = XiPi.Vect()*1000;
				auto Dir = Xi.Momentum();
				auto XVM = Xi.Momentum();
				Dir = Dir* (1./Dir.Mag());
				z4 = z3+Dir.Z()*1000;
				x4 = x3+Dir.X()*1000;
				y4 = y4+Dir.Y()*1000;
				cout<<Form("Xi Vertex (%f,%f,%f) Mass: %f",XV.X(),XV.Y(),XV.Z(),Xi.Mass())<<endl;
				cout<<Form("Xi Mom (%f,%f,%f) mag: %f",XVM.X(),XVM.Y(),XVM.Z(),XVM.Mag())<<endl;
//				cout<<"XiPiID = "<<XiPiID<<endl;
//				cout<<Form("XiPiMom = (%f,%f,%f), mag = %f MeV/c, cd = %f mm",XiPiP.X(),XiPiP.Y(),XiPiP.Z(),HelixTrackMom.at(XiPiID),Xi.GetCD())<<endl;
//				cout<<"XiPiTrackID = "<<HelixTrackID.at(XiPiID)<<endl;
//				cout<<Form("XiCor Vertex (%f,%f,%f) Mass: %f",XVCor.X(),XVCor.Y(),XVCor.Z(),XiCor.Mass())<<endl;
				auto XiVert = new TPolyMarker3D();
				XiVert->SetPoint(0,z3,x3,y3);
				XiVert->SetMarkerColor(kCyan);
				XiVert->SetMarkerStyle(39);
				XiVert->SetMarkerSize(2);
				Vertex3d.push_back(XiVert);
//				cout<<"PropDist : "<<(LV-XV).Mag()<<endl;
				auto XiTrack = new TPolyLine3D(2);
				XiTrack->SetPoint(0,z3,x3,y3);
				XiTrack->SetPoint(1,z4,x4,y4);
				XiTrack->SetLineColor(kCyan);
				XiTrack->SetLineWidth(2);
//				VertexTrack3D.push_back(XiTrack);
				HelixTrack3D.at(HelixTrackID.at(XiPiID))->SetLineColor(kCyan);
				HelixTrack3D.at(HelixTrackID.at(XiPiID))->SetLineWidth(3);
				HelixTrack3D.at(HelixTrackID.at(XiPiID))->SetLineStyle(2);
				if(XiStarSearch and XiStar.Exist()){
					auto XiStarVert = new TPolyMarker3D();
					auto XSV = XiStar.Vertex();
					auto XSM = XiStar.GetIniMom()*1000;
					double x5=XSV.x(),y5=XSV.y(),z5=XSV.z();
					double x6=XSM.x(),y6=XSM.y(),z6=XSM.z();
					x6+=x5;
					y6+=y5;
					z6+=z5;
					XiStarVert->SetPoint(0,z5,x5,y5);
					XiStarVert->SetMarkerStyle(39);
					XiStarVert->SetMarkerSize(2);
					XiStarVert->SetMarkerColor(kGreen);
					Vertex3d.push_back(XiStarVert);
					auto XiStarTrack = new TPolyLine3D(2);
					XiStarTrack->SetPoint(0,z5,x5,y5);
					XiStarTrack->SetPoint(1,z6,x6,y6);
					XiStarTrack->SetLineColor(kGreen);
					XiStarTrack->SetLineWidth(2);
//					VertexTrack3D.push_back(XiStarTrack);
					auto KMTrack = new TPolyLine3D(2);
					auto KMMom = XiStar.GetKMMom();
					double kmx=KMMom.x(),kmy=KMMom.y(),kmz=KMMom.z();
					KMTrack->SetPoint(0,z5,x5,y5);
					KMTrack->SetPoint(1,z5-kmz*1000,x5-kmx*1000,y5-kmy*1000);
					KMTrack->SetLineColor(kBlack);
					KMTrack->SetLineWidth(2);
//					VertexTrack3D.push_back(KMTrack);
					
					auto KPTrack = new TPolyLine3D(2);
					auto KPMom = XiStar.GetKPMom();
					double kpx=KPMom.x(),kpy=KPMom.y(),kpz=KPMom.z();
					KPTrack->SetPoint(0,z5,x5,y5);
					KPTrack->SetPoint(1,z5+kpz*1000,x5+kpx*1000,y5+kpy*1000);
					KPTrack->SetLineColor(kBlack);
					KPTrack->SetLineWidth(2);
//					VertexTrack3D.push_back(KPTrack);
					auto Pi0Vert = new TPolyMarker3D();
					auto Pi0V = Pi0.Vertex();
					double x7=Pi0V.x(),y7=Pi0V.y(),z7=Pi0V.z();
					Pi0Vert->SetPoint(0,z7,x7,y7);
					Pi0Vert->SetMarkerStyle(39);
					Pi0Vert->SetMarkerSize(2);
					Pi0Vert->SetMarkerColor(kBlue);
					Vertex3d.push_back(Pi0Vert);
					cout<<Form("Pi0 Mass = %f GeV",Pi0.Mass())<<endl;
				}
			}
			for(auto v :Vertex3d){
				v->Draw("same");
			}
			for(auto t :VertexTrack3D){
				t->Draw("same");
			}
}



void
TPCManager::ReconEvent(){
	vector<Vertex> verts;
	vector<Track> parts;
	vector<Track> kuramas;
	bool ldflg = false,xiflg=false;
	Ld.Clear();Xi.Clear();
	double chi_cut = 150;
	double cd_cut = 30;
	LdPiID = -1;LdProtonID = -1;XiPiID=-1;
	Track K18Track,KuramaTrack;
	for(int nt1 = 0; nt1<ntTpc;++nt1){
		if(chisqr->at(nt1)>chi_cut) continue; 
		int hf = helix_flag->at(nt1);
		if(IsAccidental(hf)) continue;
		int nh = helix_cx->size();
		double hcx = helix_cx->at(nt1);
		double hcy = helix_cy->at(nt1);
		double hz0 = helix_z0->at(nt1);
		double hr = helix_r->at(nt1);
		double hdz = helix_dz->at(nt1);
		double par1[5] = {hcx,hcy,hz0,hr,hdz};
		int id1 = pid->at(nt1);
		double q1 = charge->at(nt1);
		bool kurama = true;
		if(IsK18(hf)){
			K18Track = Track(id1,q1,par1,nt1);
			cout<<"K18"<<endl;
		}
		else if(IsKurama(hf)){
			if(kurama){
				KuramaTrack = Track(id1,q1,par1,nt1);
				cout<<"Kurama"<<endl;
				kurama = false;
			}
			if(q1 <0){
				cout<<Form("Warning! Kurama charge =%f!",q1)<<endl;
				q1 = 1.;
			}
			auto T = Track(id1,q1,par1,nt1);
			T.SetK();
			kuramas.push_back(Track(id1,q1,par1,nt1));
		}
		else{
			parts.push_back(Track(id1,q1,par1,nt1));
		}
	}
	if(XiStarSearch){
		vector<TVector3> KMX,KMP,KPX,KPP;
		for(int ikm=0;ikm<ntK18;++ikm){
			double xkm = xtgtHS[ikm];		
			double ykm = ytgtHS[ikm];		
			double zkm = ztgtHS[ikm];	
			double ukm = utgtHS[ikm];		
			double vkm = vtgtHS[ikm];
			double pkm = pHS[ikm];
			double pzkm = pkm/sqrt(1+ukm*ukm+vkm*vkm);
			double pxkm = pzkm * ukm;
			double pykm = pzkm * vkm;
//			cout<<Form("param : x,y,z,u,v,p = (%f,%f,%f,%f,%f,%f)",xkm,ykm,zkm,ukm,vkm,pkm)<<endl;
			KMX.push_back(TVector3(xkm,ykm,zkm));
			KMP.push_back(TVector3(pxkm,pykm,pzkm));
		}
		for(int ikp=0;ikp<ntKurama;++ikp){
			double xkp = xtgtKurama[ikp];		
			double ykp = ytgtKurama[ikp];		
			double zkp = ztgtKurama[ikp];	
			double ukp = utgtKurama[ikp];		
			double vkp = vtgtKurama[ikp];	
			double pkp = pKurama[ikp];
			double pzkp = pkp/sqrt(1+ukp*ukp+vkp*vkp);
			double pxkp = pzkp * ukp;
			double pykp = pzkp * vkp;
//			cout<<Form("param : x,y,z,u,v,p = (%f,%f,%f,%f,%f,%f)",xkp,ykp,zkp,ukp,vkp,pkp)<<endl;
			KPX.push_back(TVector3(xkp,ykp,zkp));
			KPP.push_back(TVector3(pxkp,pykp,pzkp));
		}
//		cout<<Form("nkm,nkp = %d, %d KMX.size()",KMX.size(),KPX.size())<<endl;
		XiStar.SetK18Track(K18Track);
		XiStar.SetKuramaTrack(KuramaTrack);
		XiStar.Construct(KMX,KMP,KPX,KPP);
		cout<<"ID"<<endl;
		cout<<K18Track.GetID()<<endl;
		cout<<KuramaTrack.GetID()<<endl;
		cout<<"ID"<<endl;
	}
	int np = parts.size();
	if(np<1) return;
	for(int nt1=0;nt1<np;++nt1){
		Vertex f(parts[nt1]);
		f.SetCdCut(cd_cut);
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

	LdProtonID = Ld.GetID1();
	LdPiID = Ld.GetID2();
	comp = 9999;
	bool KinematicFit = false;
	VertexLH V(Ld,KinematicFit);
	for(auto kurama:kuramas){
		V.AddKuramaTrack(kurama);
	}
	V.SetCdCut(cd_cut);
	for(auto p : parts){
		V.AddTrack(p);
	}
	ldflg=Ld.Exist();
	if(ldflg)V.SearchXiCombination();	
	Xi = V.GetXi();
//	XiCor = V.GetXiCor();
	XiPiID = Xi.GetID2();
	xiflg=Xi.Exist();
	if(XiStarSearch){
		Pi0=Recon(XiStar,Xi,mXiStar,mXi);
	}
}











#endif
