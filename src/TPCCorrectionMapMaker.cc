#include "TPCManager.cc"
#include "../include/TPCCorrectionMapMaker.hh"
TPCCorrectionMapMaker Corrector;
void TPCCorrectionMapMaker(){
	cout<<"TPCCorrectionMapMaker(int dum)"<<endl;
	cout<<"WriteCorrectionMap()"<<endl;
	cout<<"MakeSnakeHistogram()"<<endl;
	cout<<"WriteSnakeCorrectionMap()"<<endl;
	ROOT::EnableThreadSafety();
}
void TPCCorrectionMapMaker(int a){
	TFile* file = new TFile("CorHist.root","recreate");
	TString dir = "../../MayRun/rootfiles/Defocus/";
//	int runnum = 5755;
	int runnum = 5000;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.Process();
	file->cd();
	Corrector.WriteHist();
	file->Write();	
	file->Close();	
}
void MakeSnakeHistogram(){
	TString dir = "../../MayRun/rootfiles/Defocus/";
	int runnum = 5000;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.MakeSnakeHist();
	Corrector.MakeSnakeFile("SnakeHist_node2.root");
	Corrector.WriteSnakeHist();
}
void WriteCorrectionMap(){
	Corrector.MakeCorParameterFile("TPCCorrectionMap_KBM");
	Corrector.LoadCorrectionHist("CorHist.root");
	Corrector.WriteParam();
}
void WriteSnakeCorrectionMap(){
	Corrector.LoadSnakeHist("SnakeHist_node2.root");
	Corrector.MakeCorParameterFile("TPCCorrectionMap_Snake");
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	TGraphErrors* gr= new TGraphErrors(nbin_z);
	gr->SetMarkerStyle(5);
	gr->SetLineColor(kRed);
	gr->SetLineWidth(2);
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			vector<double>px,pxf,pxm1,pxm2,pxm3,py,pyxp,pyf;
			vector<double>ex,ex1,ex2,ey,eyxp;
			auto* hx = Corrector.ScanSnake(ix,iy,px,ex);
			cout<<Form("Scanned(%d,%d)",ix,iy)<<endl;
			cout<<"px: "<<px.size()<<endl;
			cout<<"ex: "<<ex.size()<<endl;
			delete gr;
			int np = px.size();
			int eval;
			gr= new TGraphErrors(np);
			for(int i=0;i<np;++i){
				pxf.push_back(FillEmptyParam(px,i));
			}
			cout<<"pxm1"<<endl;
			for(int i=0;i<np;++i){
				pxm1.push_back(Wmean(pxf,ex,i,eval));
			}
			cout<<"pxm2"<<endl;
			for(int i=0;i<np;++i){
				pxm2.push_back(Wmean(pxm1,ex,i,eval));
			}
			cout<<"pxm3"<<endl;
			for(int i=0;i<np;++i){
				pxm3.push_back(Wmean(pxm2,ex,i,eval));
				pyf.push_back(0);
				gr->SetPoint(i,Corrector.BinPosZ(i),pxm3[i]);
//				gr->SetPoint(i,Corrector.BinPosZ(i),pxf[i]);
				gr->SetPointError(i,1,1/ex[i]);
			}
			c1->cd();
			hx->Draw("colz");
			gr->Draw("Same");
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			cin.ignore();
			cout<<"SettingParams"<<endl;
//			RestoreEmpty(px);
			Corrector.SetParameter(pxm3,pyf,ix,iy);
		}
	}
//	Corrector.SmoothenParam();
	Corrector.WriteParam();
}




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
TH2D* TPCCorrectionMapMaker::ScanSnake(int ix,int iy,vector<double> &peaks,vector<double> &ents){
	auto* hist = GetSnakeHist(ix,iy,0);
	cout<<"Hist: "<<hist->GetEntries()<<endl;
	int cut=200;
	TF1* f = new TF1("f","gaus",-5,5);
	bool active = true;
	for(int i=1;i<snake_nbin_z+1;i++){
		double pos_z = BinPosZ(i);
		double pos_x = BinPosX(ix);
//		active = IsActiveArea(pos_z,pos_x);
		TString Key = (TString)hist->GetTitle()+Form("_Proj%d",i);
		TH1D* h = hist->ProjectionY(Key,i,i);
		double peakcenter = GetPeakPosition(h);
		cout<<Form("%d: peakcenter %f",i,peakcenter)<<endl;
		f->SetParLimits(1,-5,5);
		h->Fit("f","QR");
		ents.push_back(h->GetEffectiveEntries());
		if(!active){ //			peaks.push_back(-9999);
			cout<<"DeadArea!"<<endl;
			peaks.push_back(nand);
			ents[i-1]=1;
		}
		else if(ents[i-1]>cut ){
			peaks.push_back(f->GetParameter(1));
		}
		else{
			ents[i-1]=1;
			peaks.push_back(nand);
		}
	}
	return hist;
}
void TPCCorrectionMapMaker::MakeSnakeFile(TString FileName){
	OutFile= new TFile(FileName,"recreate");
	WriteTag("X_seg",Form("( nbin =%d, center =%f, width =%f)",nbin_x,center_x,width_x));
	WriteTag("Y_seg",Form("( nbin =%d, center =%f, width =%f)",nbin_y,center_y,width_y));
}
void TPCCorrectionMapMaker::RestoreEmptyParams(){
	double peak_temp[nbin_x][nbin_y][nbin_z];
	for(int ix=0;ix<nbin_x;++ix){			
		for(int iy=0;iy<nbin_y;++iy){			
			for(int iz=0;iz<nbin_z;++iz){
				peak_temp[ix][iy][iz]=nand;
				int bx = ix-1,by=iy-1,bz=iz-1;
				int fx = ix+1,fy=iy+1,fz=iz+1;
//				double bpx,bpy,bpz,fpx,fpy,fpz;
				double bpx,fpx;
				int bzd=0,fzd=0;
				if(isnan(peakx[ix][iy][iz])){
					while(bz>-1){
						bzd++;
						if(!isnan(peakx[ix][iy][bz])){
							bpx+=peakx[ix][iy][bz];
							break;
						}
						bz--;
					}
					if(bz==-1) bzd=0;
					while(fz<nbin_z){
						fzd++;
						if(!isnan(peakx[ix][iy][fz])){
							fpx+=peakx[ix][iy][fz];
							break;
						}
						fz++;
					}
					if(fz==nbin_z) fzd=0;
					peak_temp[ix][iy][iz]=(bpx*fzd+fpx*bzd)/(fzd+bzd);	
				};
			}
		}
	}
}
