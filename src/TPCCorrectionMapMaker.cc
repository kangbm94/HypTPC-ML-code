#include "TPCManager.cc"
#include "../include/TPCCorrectionMapMaker.hh"
TPCCorrectionMapMaker Corrector;
void TPCCorrectionMapMaker::AssignHit(TVector3& Pos, TVector3& Cor,int nh){
	double x=clxTpc->at(nh),y=clyTpc->at(nh),z=clzTpc->at(nh);	
	double corx=xCorVec->at(nh),cory=yCorVec->at(nh),corz=zCorVec->at(nh);	
	Pos=TVector3(x,y,z);Cor=TVector3(corx,cory,corz);
}
void TPCCorrectionMapMaker::AssignHits(TVector3* Pos, TVector3* Cor){
	int nhits = GetNhits();
	for(int i=0;i<max_hit;++i){
		Pos[i]={0,0,0};Cor[i]={0,0,0};
	}
	if(ntBcOut==1){
		for(int i=0;i<nhits;++i){
			AssignHit(Pos[i],Cor[i],i);
		}
	}
}
void TPCCorrectionMapMaker::Process(){
	cout<<"Processing..."<<endl;
	int ent = DataChain->GetEntries();
	TVector3 Pos[max_hit],Cor[max_hit];
	for(int i=0;i<ent;++i){
		SetEvent(i);
		AssignHits(Pos,Cor);
		int nh = GetNhits();
		Indicator(i,ent);
		for(int j=0;j<nh;++j){
			FillHistograms(Pos[j],Cor[j]);
		}
	}
}
void TPCCorrectionMapMaker(){
	cout<<"TPCCorrectionMapMaker(int dum)"<<endl;
	cout<<"WriteCorrectionMap()"<<endl;
}
void TPCCorrectionMapMaker(int a){
	TFile* file = new TFile("CorHist.root","recreate");
	TString dir = "../../MayRun/rootfiles/Defocus/";
	int runnum = 5000;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.Process();
	file->cd();
	Corrector.WriteHist();
	file->Write();	
	file->Close();	
}
void WriteCorrectionMap(){
	Corrector.MakeCorParameterFile("TPCCorrectionMap_KBM");
	Corrector.LoadHist("CorHist.root");
	Corrector.WriteParam();
}
