#include "../include/ReconTools.hh"
#ifndef ReconTools_C
#define ReconTools_C
#include "KinFit.cc"
namespace{
	TVector3 Mom(double x,double y,double u,double v,double p){
		double pt = p/sqrt(1+u*u+v*v);
		return TVector3(pt*u,pt*v,pt);
	}
	TVector3 Pos(double x,double y,double z){
		return TVector3(x,y,z);
	}
	bool InTarget(TVector3 pos,double cd){
		if(abs(pos.X())<20 
				and abs(pos.Y())<25
				and abs(pos.Z()+30)<100
				and cd < 100) return true;
		else return false;
	}
}
Recon::Recon(vector<TLorentzVector> D,TVector3 vertex,double clos_dist,int id1,int id2, int charge_){
	charge = charge_;
	LV.SetXYZM(0,0,0,0);
	Daughters = D;
	Vert=vertex;
	close_dist=clos_dist;
	for(auto lv:D) LV+=lv;
	CombID = pow(2,id1)+pow(2,id2);
	exist = true;
	auto V_t = GlobalToTarget(Vert);
	auto mom = LV.Vect();
	auto dir_t = GlobalToTargetMom(mom);
	dir_t = dir_t* (1/dir_t.Y());
	auto u = dir_t.X();
	auto v = dir_t.Z();
	par[0]=V_t.X()-V_t.Z()*u,par[1]=V_t.Z()-V_t.Y()*v,par[2]= u,par[3]=v;
	if(charge){
		GetHelixParameter(Vert,mom,charge,par);
	}
	trid1=id1;trid2=id2;
}

bool Vertex::AddTrack(Track p){
	int np = NTrack();
	int nt = p.GetID();
	if(Counted(p)) return false;
	if(np==1){
		auto par1 = Tracks[0].GetPar();
		auto par2 = p.GetPar();
		double cd,t1,t2;
		auto pos = VertexPointHelix(par1,par2,cd,t1,t2);
		if(cd<cdcut){
			Tracks.push_back(p);
			verts.push_back(pos);
			Vert_id+=pow(2,nt);
			SetVert();
			return true;
		}
	}
	else{
		double t = GetTcal(p.GetPar(),vert);
		auto point = HelixPos(p.GetPar(),t);
		double cd = (vert-point).Mag();
		double ct,t1,t2;
		if(cd<cdcut*1.5){
			vector<TVector3> cand;
			for(auto pt:Tracks){
				cand.push_back(VertexPointHelix(pt.GetPar(),p.GetPar(),ct,t1,t2));
			}
			TVector3 vec(0,0,0);
			for(auto v:cand){vec+=v*(1./cand.size());	
				Tracks.push_back(p);Vert_id+=pow(2,nt);verts.push_back(vec); SetVert();
			}
			return true;
		}
	}
	return false;
}
void Vertex::SearchLdCombination(){
	int np = NTrack();
	vector<Track>PCand;PCand.clear();
	vector<Track>PiCand;PiCand.clear();
	for(auto p:Tracks){
		if(p.IsP())PCand.push_back(p);
		if(p.IsPi() and p.IsNegative())PiCand.push_back(p);
	}
	double mom_cut = 0.9;//1.2 for XiStar, 0.9 for Xi
	double mom_cutMin = 0.2;//0.3 for XiStar, 0.2 for Xi
	for(auto p:PCand){
		for(auto pi:PiCand){
			double cd_,t1_,t2_;
			if(p.GetID() == pi.GetID()) continue;
			auto ppivert = VertexPointHelix(p.GetPar(),pi.GetPar(),cd_,t1_,t2_); 
			if(cd_>cdcut) continue;
			auto p1 = CalcHelixMom(p.GetPar(),ppivert.y());
			auto p2 = CalcHelixMom(pi.GetPar(),ppivert.y());
			auto pLV = TLorentzVector(p1,sqrt(mp*mp+p1.Mag2()));
			auto piLV = TLorentzVector(p2,sqrt(mpi*mpi+p2.Mag2()));
//			vector<TLorentzVector> lv1 = {pLV,piLV};
//			auto ldlv1 = pLV+piLV;
//			if(ldlv1.Vect().Mag()<mom_cut and ldlv1.Vect().Mag()>mom_cutMin)LdCand.push_back(Recon(lv1,ppivert,cd_,p.GetID(),pi.GetID()));
			auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv2 = {pLV,piLVInv};
			auto ldlv2 = pLV+piLVInv;
			if(ldlv2.Vect().Mag()<mom_cut and ldlv2.Vect().Mag()>mom_cutMin)LdCand.push_back(Recon(lv2,ppivert,cd_,p.GetID(),pi.GetID()));
		}
	}
}

bool VertexLH::AddTrack(Track p){
	int np = NTrack();
	int nt = p.GetID();
	if(Recons[0].Counted(p)) return false;
	if(!(Recons[0].Exist())) return false;
	if(np==1 or np!=1){
		auto par1 = p.GetPar();
		auto par2 = Recons[0].GetPar();
		double cd,t1,t2;
		auto vv = Recons[0].Vertex();
		auto pos = VertexPointHelixLinear(par1,par2,cd,t1,t2);
		if(cd<cdcut){
				auto prop = vv-pos;
//				cout<<Form("Cd=%f, PropDist(%f,%f,%f)",cd,prop.X(),prop.Y(),prop.Z())<<endl;
				//	cout<<Form("LdVert(%f,%f,%f)",vv.X(),vv.Y(),vv.Z())<<endl;
		//	cout<<Form("HLVert(%f,%f,%f)",pos.X(),pos.Y(),pos.Z())<<endl;
			Tracks.push_back(p);
			verts.push_back(pos);Vert_id+=pow(2,nt);SetVert();
			return true;
		}
	}
	else{
		double t = GetTcal(p.GetPar(),vert);
		auto point = HelixPos(p.GetPar(),t);
		double cd = (vert-point).Mag();
		double ct,t1,t2;
		if(cd<cdcut){
			vector<TVector3> cand;
			for(auto pt:Tracks){
				cand.push_back(VertexPointHelix(pt.GetPar(),p.GetPar(),ct,t1,t2));
			}
			TVector3 vec(0,0,0);
			for(auto v:cand){
				vec+=v*(1./cand.size());	
				Tracks.push_back(p);Vert_id+=pow(2,nt);verts.push_back(vec); SetVert();
			}
			return true;
		}
	}
	return false;
}
void VertexLH::SearchXiCombination(){
	int np = NTrack();
	LdCand.clear();
	vector<Track>PiCand;PiCand.clear();
	for(auto p:Tracks){
		if(p.IsPi() and p.IsNegative())PiCand.push_back(p);
	}
	double mom_cut = 0.9;//Xi->0.9 XiStar->1.2
	double mom_cutMin = 0.4;//Xi-> 0.4 XiStar->0.8
	for(auto ld:Recons) LdCand.push_back(ld);
	for(auto ld:LdCand){
		for(auto pi:PiCand){
			double cd_,t1_,t2_;
			if(ld.Counted(pi)) continue;
			auto ldpivert = VertexPointHelixLinear(pi.GetPar(),ld.GetPar(),cd_,t1_,t2_); 
			TVector3 ldmomdir = ld.Momentum() * (1./ld.Momentum().Mag());
			TVector3 lddir = ld.Vertex() - ldpivert;
			lddir = lddir * (1./lddir.Mag());
			if(ldmomdir * lddir < 0) continue;
			auto p2 = CalcHelixMom(pi.GetPar(),ldpivert.y());
			std::bitset<8>ldb(ld.GetID());
//			cout<<"cd = "<<cd_<<" PiID,LdId = ("<<pi.GetID()<<" , "<<ldb<<" )"<<endl;
			auto ldLV = ld.GetLV();

			auto piLV = TLorentzVector(p2,sqrt(mpi*mpi+p2.Mag2()));
//			vector<TLorentzVector> lv1 = {ldLV,piLV};
//			auto xilv1 = ldLV+piLV;
//			if(xilv1.Vect().Mag()<mom_cut and xilv1.Vect().Mag()>mom_cutMin)XiCand.push_back(Recon(lv1,ldpivert,cd_,ld.GetID(),pi.GetID(),-1));
			auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv2 = {ldLV,piLVInv};
			auto xilv2 = ldLV+piLVInv;
			if(xilv2.Vect().Mag()<mom_cut and xilv2.Vect().Mag()>mom_cutMin)XiCand.push_back(Recon(lv2,ldpivert,cd_,ld.GetID(),pi.GetID(),-1));
		
			
/*
			auto ldPLV = ld.GetDaughter(0);
			auto ldPP = ldPLV.Vect(); 
			auto ldPPres = ldPP * (0.03);
			
			auto ldpiLV = ld.GetDaughter(1);
			auto ldpiP = ldpiLV.Vect(); 
			auto ldpiPres = ldpiP * (0.03);
			Fitter.AssignLorentzVector(ldPLV,ldpiLV);
			Fitter.SetResolution(ldPPres,ldpiPres);
			
			Fitter.DoKinematicFit();
			auto ldPCor = Fitter.GetFittedLV().at(0);
			auto ldpiCor = Fitter.GetFittedLV().at(1);
			auto ldLVCor = ldPCor+ldpiCor;

			auto xiCorlv1 = ldLVCor+piLV;
			vector<TLorentzVector> lv1Cor = {ldLVCor,piLV};
			if(xiCorlv1.Vect().Mag()<mom_cut and xiCorlv1.Vect().Mag()>mom_cutMin)XiCorCand.push_back(Recon(lv1Cor,ldpivert,cd_,ld.GetID(),pi.GetID()));
			auto xiCorlv2 = ldLVCor+piLVInv;
			vector<TLorentzVector> lv2Cor = {ldLVCor,piLVInv};
			if(xiCorlv2.Vect().Mag()<mom_cut and xiCorlv2.Vect().Mag()>mom_cutMin)XiCorCand.push_back(Recon(lv2Cor,ldpivert,cd_,ld.GetID(),pi.GetID()));
		*/
		}
	}
}

void VertexLH::AddKuramaTrack(Track p){
	KuramaFlag = true;
	KuramaTracks.push_back(p);
}

XiStarRecon::XiStarRecon(vector<TVector3>KMX,vector<TVector3>KMP,vector<TVector3>KPX,vector<TVector3>KPP){
	int ntK18 = KMX.size();
	int ntKurama = KPX.size();
	charge = -1;
	vector<double>close_dists;
	vector<double>MMs;
	vector<TLorentzVector>LVs;
	vector<TVector3> verts;
	for(int ikm = 0;ikm<ntK18;++ikm){
		auto xkm = KMX.at(ikm);
		auto pkm = KMP.at(ikm);
		for(int ikp = 0;ikp<ntKurama;++ikp){
			auto xkp = KPX.at(ikp);
			auto pkp = KPP.at(ikp);
			auto kkvert = VertexPoint(xkm, xkp, pkm, pkp);
			double cd = CloseDist(xkm,xkp,pkm,pkp); 
			if(!InTarget(kkvert,cd))continue;
			TLorentzVector LVKM(pkm,sqrt(mk*mk+pkm*pkm));
			TLorentzVector LVKP(pkp,sqrt(mk*mk+pkp*pkp));
			TLorentzVector LVTarget(0,0,0,mp);
			TLorentzVector LVMM = LVKM+LVTarget-LVKP;
			double MM = LVMM.Mag();
			verts.push_back(kkvert);
			close_dists.push_back(cd);
			LVs.push_back(LVMM);
			MMs.push_back(MM);
		}
	}
	double dif = 1e9;
	int id = -1;
	for(int im=0;im<MMs.size();++im){
		if(dif >abs( MMs.at(im) - mXiStar)){
			dif = abs(MMs.at(im)-mXiStar);
			id = im;
		}
	}
	if(id>0){
		exist = true;
		auto vert_ = verts.at(id);
		double x = vert_.X();	
		double y = vert_.Y();	
		double z = vert_.Z();	
		z=ZTarget;
		Vert= TVector3(x,y,z);
		auto LV = LVs.at(id);	
		auto mom = LV.Vect();
		LV.SetXYZM(mom.X(),mom.Y(),mom.Z(),mXiStar);
		GetHelixParameter(Vert,mom,charge,par);
	}
}
VertexMM::VertexMM(XiStarRecon XiStar,Recon Xi){
	auto XiStarPar =	XiStar.GetPar();
	auto XiPar =	Xi.GetPar();
	double cd,t1,t2;
	vert  = VertexPointHelix(XiStarPar,XiPar,cd,t1,t2);
	auto XiStarp = CalcHelixMom(XiStarPar,vert.y());
	auto Xip = CalcHelixMom(XiPar,vert.y());
	TLorentzVector LVXiStar(XiStarp, sqrt(XiStarp.Mag2()+mXiStar*mXiStar));	
	TLorentzVector LVXi(Xip, sqrt(Xip.Mag2()+mXi*mXi));	

	auto LVPi0 = LVXiStar - LVXi;
	Recon Pi0;
	Pi0.SetLV(LVPi0);
	Pi0.SetExistance(true);
}
#endif
