#include "../include/ReconTools.hh"
#ifndef ReconTools_C
#define ReconTools_C
#include "KinFit.cc"
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
		if(cd<cdcut){
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
		if(p.IsPi())PiCand.push_back(p);
	}
	double mom_cut = 0.9;//1.2 for XiStar, 0.9 for Xi
	double mom_cutMin = 0.2;//0.3 for XiStar, 0.2 for Xi
	for(auto p:PCand){
		for(auto pi:PiCand){
			double cd_,t1_,t2_;
			if(p.GetID() == pi.GetID()) continue;
			auto ppivert = VertexPointHelix(p.GetPar(),pi.GetPar(),cd_,t1_,t2_); 
			auto p1 = CalcHelixMom(p.GetPar(),ppivert.y());
			auto p2 = CalcHelixMom(pi.GetPar(),ppivert.y());
			auto pLV = TLorentzVector(p1,sqrt(mp*mp+p1.Mag2()));
			auto piLV = TLorentzVector(p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv1 = {pLV,piLV};
			auto ldlv1 = pLV+piLV;
			if(ldlv1.Vect().Mag()<mom_cut and ldlv1.Vect().Mag()>mom_cutMin)LdCand.push_back(Recon(lv1,ppivert,cd_,p.GetID(),pi.GetID()));
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
		if(p.IsPi())PiCand.push_back(p);
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
			vector<TLorentzVector> lv1 = {ldLV,piLV};
			auto xilv1 = ldLV+piLV;
			if(xilv1.Vect().Mag()<mom_cut and xilv1.Vect().Mag()>mom_cutMin)XiCand.push_back(Recon(lv1,ldpivert,cd_,ld.GetID(),pi.GetID()));
			auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv2 = {ldLV,piLVInv};
			auto xilv2 = ldLV+piLVInv;
			if(xilv2.Vect().Mag()<mom_cut and xilv2.Vect().Mag()>mom_cutMin)XiCand.push_back(Recon(lv2,ldpivert,cd_,ld.GetID(),pi.GetID()));
		
			
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
		}
	}
}
#endif
