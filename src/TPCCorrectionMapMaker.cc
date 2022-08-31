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
			FillHists(Pos[j],Cor[j]);
		}
	}
}
void TPCCorrectionMapMaker::ScanSnake(TH2D* hist,vector<double> &peaks,vector<double> &ents){
	TF1* f = new TF1("f","gaus",-5,5);
	for(int i=1;i<snake_nbin_z+1;i++){
		TString Key = (TString)hist->GetTitle()+Form("_Proj%d",i);
		TH1D* h = hist->ProjectionY(Key,i,i);
		double peakcenter = GetPeakPosition(h);
		cout<<Form("%d: peakcenter %f",i,peakcenter)<<endl;
		f->SetParLimits(1,-5,5);
		h->Fit("f","QR");
		ents.push_back(h->GetEntries());
		if(ents[i-1]>1000){
			peaks.push_back(f->GetParameter(1));
		}
		else{
			ents[i-1]=0;
			peaks.push_back(0);
		}
	}
}
void TPCCorrectionMapMaker(){
	cout<<"TPCCorrectionMapMaker(int dum)"<<endl;
	cout<<"WriteCorrectionMap()"<<endl;
	cout<<"SnakeCorrector()"<<endl;
	cout<<"WriteSnakeCorrectionMap()"<<endl;
}
void TPCCorrectionMapMaker(int a){
	TFile* file = new TFile("CorHist.root","recreate");
	TString dir = "../../MayRun/rootfiles/Defocus/old";
	int runnum = 5755;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.Process();
	file->cd();
	Corrector.WriteHist();
	file->Write();	
	file->Close();	
}
void SnakeCorrector(){
	TString dir = "../../MayRun/rootfiles/Defocus/old/";
	int runnum = 5000;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.MakeSnakeHist();
	Corrector.MakeOutFile("SnakeHist.root");
	Corrector.WriteSnakeHist();
}
void WriteCorrectionMap(){
	Corrector.MakeCorParameterFile("TPCCorrectionMap_KBM");
	Corrector.LoadCorrectionHist("CorHist.root");
	Corrector.WriteParam();
}
void WriteSnakeCorrectionMap(){
	Corrector.LoadSnakeHist("SnakeHist.root");
	Corrector.MakeCorParameterFile("TPCCorrectionMap_Snake");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
//	TH2D* fitHist = new TH2D("FitHist","FitHist",nbin_z,0,nbin_z,100,-5,5);
	TGraph* gr= new TGraph(nbin_z);
	gr->SetMarkerStyle(3);
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			cout<<Form("(%d,%d)",ix,iy)<<endl;
			vector<double>px,pxxp,pxf,py,pyxp,pyf;
			vector<double>ex,exxp,ey,eyxp;
			auto* hx = Corrector.GetSnakeHist(ix,iy,1);
			Corrector.ScanSnake(hx,px,ex);
			auto* hy = Corrector.GetSnakeHist(ix,iy,0);
			Corrector.ScanSnake(hy,py,ey);
			cout<<Form("Scanned(%d,%d)",ix,iy)<<endl;
			cout<<"px: "<<px.size()<<endl;
			ExtrapolateEmptyVector(px,pxxp);
			cout<<"ex: "<<ex.size()<<endl;
			ExtrapolateEmptyVector(ex,exxp);
			cout<<"py: "<<py.size()<<endl;
			ExtrapolateEmptyVector(py,pyxp);
			cout<<"ey: "<<ey.size()<<endl;
			ExtrapolateEmptyVector(ey,eyxp);
			cout<<Form("Extrapolated(%d,%d)",ix,iy)<<endl;
			cout<<pxxp.size()<<endl;
			cout<<exxp.size()<<endl;
			delete gr;
			gr= new TGraph(snake_nbin_z-1);
			for(int i=0;i<px.size();++i){
//				pxf.push_back(Wmean(pxxp,exxp,i));
//				pyf.push_back(Wmean(pyxp,eyxp,i));
				pxf.push_back(Wmean(px,exxp,i));
				pyf.push_back(Wmean(py,eyxp,i));
				gr->SetPoint(i,Corrector.BinPosZ(i),pxf[i]);
//				gr->SetPoint(i,Corrector.BinPosZ(i),px[i]);
			}
			c1->cd();
			hx->Draw("colz");
			gr->Draw("Same");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"SettingParams"<<endl;
			RestoreEmpty(px);
			RestoreEmpty(py);
//			Corrector.SetParameter(px,py,ix,iy);
//			Corrector.SetParameter(pxxp,pyxp,ix,iy);
			Corrector.SetParameter(pxf,pyf,ix,iy);
		}
	}
	Corrector.WriteParam();
	//	Corrector.ScanSnake("SnakeHist.root");
}
