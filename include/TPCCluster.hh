#include "TPCPadHelper.hh"
#include "Utils.hh"
#ifndef TPCCluster_h
#define TPCCluster_h
using namespace tpc;
class TPCHit{
	protected:
		int m_padID=-1;
		int m_layer=-1;
		double m_row=-1;
		double m_de=0;
		TVector3 m_position = TVector3(0,0,0);
		int m_cluster = -1;
//		double m_phi;
	public:
		TPCHit(){};
		TPCHit(int padID,double de){
			m_padID=padID;
			m_layer = getLayerID(padID); 
			m_row = getRowID(padID); 
			m_position = getPosition(padID);
			m_de = de;
//			cout<<"Row: "<<m_row<<endl;
//			cout<<Form("P=%d,L=%d,R=%f",m_padID,m_layer,m_row)<<endl;
	//		m_phi = ;
			m_cluster = -1;
		}
		double GetPadID(){
			return m_padID;
		}
		int GetLayer(){
			return m_layer;
		}
		double GetDe(){
			return m_de;
		}

		TVector3 GetPosition(){
			return m_position;
		}
		double GetRow(){
			return m_row;
		}
		void SetCluster(int clnum){
			m_cluster = clnum;
		}
		int GetCluster(){
			return m_cluster;
		}

		virtual void Show(){
				cout<<Form("L%d Position: (%f,%f,%f)",m_layer,m_position.X(),m_position.Y(),m_position.Z())<<endl;
		}
};

class TPCCluster:public TPCHit{
	protected:
		vector<TPCHit> m_hits;

	public:
		TPCCluster(vector<TPCHit> hits){
			m_hits = hits;
			for(auto hit : m_hits){
				m_layer = hit.GetLayer();
				m_de+= hit.GetDe();
				auto pos = hit.GetPosition();	
				m_position+=hit.GetDe()*pos;	
			}
			m_position*=1/m_de;
//			cout<<m_hits.size()<<endl;
	}
		virtual void Show(){
				cout<<Form(" L%d Size: %lu, Position: (%f,%f,%f)",m_layer,m_hits.size(),m_position.X(),m_position.Y(),m_position.Z())<<endl;
		}
		int GetSize(){
			return m_hits.size();
		}
		TPCHit GetHit(int i){
			return m_hits[i];
		}
};

#endif
