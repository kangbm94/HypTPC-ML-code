#include "../include/ReconTools.hh"
#ifndef ReconTools_C
#define ReconTools_C
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
	double mom_cut = 1;
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
			if(ldlv1.Vect().Mag()<mom_cut)LdCand.push_back(Recon(lv1,ppivert,p.GetID(),pi.GetID()));
			auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv2 = {pLV,piLVInv};
			auto ldlv2 = pLV+piLVInv;
			if(ldlv2.Vect().Mag()<mom_cut)LdCand.push_back(Recon(lv2,ppivert,p.GetID(),pi.GetID()));
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
	double mom_cut = 0.9;
	for(auto ld:Recons) LdCand.push_back(ld);
	for(auto ld:LdCand){
		for(auto pi:PiCand){
			double cd_,t1_,t2_;
			if(ld.Counted(pi)) continue;
			auto ldpivert = VertexPointHelixLinear(pi.GetPar(),ld.GetPar(),cd_,t1_,t2_); 
			auto p2 = CalcHelixMom(pi.GetPar(),ldpivert.y());
			std::bitset<8>ldb(ld.GetID());
			cout<<"cd = "<<cd_<<" PiID,LdId = ("<<pi.GetID()<<" , "<<ldb<<" )"<<endl;
			auto ldLV = ld.GetLV();
			auto piLV = TLorentzVector(p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv1 = {ldLV,piLV};
			auto xilv1 = ldLV+piLV;
			if(xilv1.Vect().Mag()<mom_cut)XiCand.push_back(Recon(lv1,ldpivert,ld.GetID(),pi.GetID()));
			auto piLVInv = TLorentzVector(-p2,sqrt(mpi*mpi+p2.Mag2()));
			vector<TLorentzVector> lv2 = {ldLV,piLVInv};
			auto xilv2 = ldLV+piLVInv;
			if(xilv2.Vect().Mag()<mom_cut)XiCand.push_back(Recon(lv2,ldpivert,ld.GetID(),pi.GetID()));
		}
	}
}
#endif
