#include "../include/GeantManager.hh"
#ifndef GeantManager_c
#define GeantManager_c
GeantManager gGeantManager;
void GeantManager::Initialize(){
	DataChain->SetBranchAddress("evnum",&evnum);	

//TPC
	DataChain->SetBranchAddress("nhittpc",&nhittpc);	
	DataChain->SetBranchAddress("ntrk",ntrk);	
	DataChain->SetBranchAddress("ititpc",ititpc);	
	DataChain->SetBranchAddress("idtpc",idtpc);	
	DataChain->SetBranchAddress("xtpc",xtpc);	
	DataChain->SetBranchAddress("ytpc",ytpc);	
	DataChain->SetBranchAddress("ztpc",ztpc);	
	DataChain->SetBranchAddress("x0tpc",x0tpc);	
	DataChain->SetBranchAddress("y0tpc",y0tpc);	
	DataChain->SetBranchAddress("z0tpc",z0tpc);	
	DataChain->SetBranchAddress("pxtpc",pxtpc);	
	DataChain->SetBranchAddress("pytpc",pytpc);	
	DataChain->SetBranchAddress("pztpc",pztpc);	
	DataChain->SetBranchAddress("dedxtpc",dedxtpc);	
	DataChain->SetBranchAddress("edeptpc",edeptpc);	
	DataChain->SetBranchAddress("ipadtpc",ipadtpc);	
	DataChain->SetBranchAddress("laytpc",laytpc);	
	DataChain->SetBranchAddress("rowtpc",rowtpc);	
//TPC
//HTOF
	DataChain->SetBranchAddress("nhHtof",&nhHtof);
	DataChain->SetBranchAddress("tidHtof",tidHtof);
	DataChain->SetBranchAddress("pidHtof",pidHtof);
	DataChain->SetBranchAddress("qHtof",qHtof);
	DataChain->SetBranchAddress("xHtof",xHtof);
	DataChain->SetBranchAddress("yHtof",yHtof);
	DataChain->SetBranchAddress("zHtof",zHtof);
	DataChain->SetBranchAddress("pxHtof",pxHtof);
	DataChain->SetBranchAddress("pyHtof",pyHtof);
	DataChain->SetBranchAddress("pzHtof",pzHtof);
	DataChain->SetBranchAddress("deHtof",deHtof);
	DataChain->SetBranchAddress("tHtof",tHtof);
	DataChain->SetBranchAddress("lengthHtof",lengthHtof);
//HTOF
//FTOF
	DataChain->SetBranchAddress("nhFtof",&nhFtof);
	DataChain->SetBranchAddress("tidFtof",tidFtof);
	DataChain->SetBranchAddress("pidFtof",pidFtof);
	DataChain->SetBranchAddress("qFtof",qFtof);
	DataChain->SetBranchAddress("xFtof",xFtof);
	DataChain->SetBranchAddress("yFtof",yFtof);
	DataChain->SetBranchAddress("zFtof",zFtof);
	DataChain->SetBranchAddress("pxFtof",pxFtof);
	DataChain->SetBranchAddress("pyFtof",pyFtof);
	DataChain->SetBranchAddress("pzFtof",pzFtof);
	DataChain->SetBranchAddress("deFtof",deFtof);
	DataChain->SetBranchAddress("tFtof",tFtof);
	DataChain->SetBranchAddress("lengthFtof",lengthFtof);
//FTOF

};
int GeantManager::GetNtTPC(){
	if(nttpc!= -1) return nttpc;
	for(int ih = 0; ih < nhittpc;++ih){
		nttpc = Max(nttpc,ntrk[ih]);
	}
}
void GeantManager::ConstructTracks(){
	int nt = GetNtTPC();
	Tracks.resize(nt);
	for(int it = 0;it < nt; ++it){
			
	}
}
#endif
