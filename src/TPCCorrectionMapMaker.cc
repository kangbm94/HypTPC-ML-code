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
//	runnum = 5754;
	TString filename = Form("run0%d_DSTTPCBcOut.root",runnum);
	Corrector.LoadTPCBcOut(dir+filename);
	Corrector.MakeSnakeHist();
	Corrector.MakeSnakeFile("SnakeHist.root");
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
	TGraphErrors* gr= new TGraphErrors(snake_nbin_z);
	gr->SetMarkerStyle(5);
	gr->SetLineColor(kRed);
	gr->SetLineWidth(2);
	for(int ix=0;ix<nbin_x;++ix){
		for(int iy=0;iy<nbin_y;++iy){
			vector<double>px,pxf,pxm1,pxm2,pxm3,py,pyxp,pyf;
			vector<double>ex,ex1,ex2,ey,eyxp;
			auto* hx = Corrector.ScanSnake(ix,iy,px,ex);
			cout<<Form("Scanned(%d,%d)",ix,iy)<<endl;
			delete gr;
			int np = px.size();
			int eval;
			gr= new TGraphErrors(np);
			for(int i=0;i<np;++i){
				pxf.push_back(FillEmptyParam(px,i));
				pyf.push_back(0);
			}
			for(int i=0;i<np;++i){
				pxm1.push_back(Wmean(pxf,ex,i,eval));
			}
			for(int i=0;i<np;++i){
				pxm2.push_back(Wmean(pxm1,ex,i,eval));
			}
			for(int i=0;i<np;++i){
				pxm3.push_back(Wmean(pxm2,ex,i,eval));
				pyf.push_back(0);
//				gr->SetPoint(i,Corrector.BinPosZ(i),pxm3[i]);
				gr->SetPoint(i,Corrector.BinPosZ(i),pxf[i]);
				gr->SetPointError(i,1,1/sqrt(ex[i]));
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
			Corrector.SetParameter(pxf,pyf,ix,iy);
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
	int nhits = GetNhits(1);
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
		int nh = GetNhits(1);
		Indicator(i,ent);
		for(int j=0;j<nh;++j){
			FillHists(Pos[j],Cor[j]);
		}
	}
}
TF1* f = new TF1("f","gaus",-5,5);
TH2D* TPCCorrectionMapMaker::ScanSnake(int ix,int iy,vector<double> &peaks,vector<double> &ents){
	auto* hist = GetSnakeHist(ix,iy,0);
//	cout<<"Hist: "<<hist->GetEntries()<<endl;
	int cut=300;
	bool active = true;
	for(int i=0;i<nbin_z;i++){
		double pos_z = BinPosZ(i);
		double pos_x = BinPosX(ix);
//		active = IsActiveArea(pos_z,pos_x);
		TString Key = (TString)hist->GetTitle()+Form("_Proj%d",i);
		TH1D* h = hist->ProjectionY(Key,i+1,i+1);
//		double peakcenter = GetPeakPosition(h);
		ents.push_back(h->GetEffectiveEntries());
		double peakcenter = h->GetBinCenter(h->GetMaximumBin());
		if(abs(peakcenter)>7.1&&abs(pos_x)>150)peakcenter=0;
		//	f->SetParLimits(1,peakcenter-1,peakcenter+1);
		f->SetParameter(1,peakcenter);
		f->SetParLimits(1,peakcenter-2,peakcenter+2);
		f->SetRange(peakcenter-5,peakcenter+5);
		h->Fit("f","QR0");
		cout<<i<<"-> ent : "<<ents[i]<<", peak : "<<peakcenter<<endl; 	
		if(!active){ //			peaks.push_back(-9999);
			cout<<"DeadArea!"<<endl;
			peaks.push_back(nand);
			ents[i]=1;
		}
		else if(ents[i]>cut ){
			peaks.push_back(f->GetParameter(1));
		}
		else{
			ents[i]=1;
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
