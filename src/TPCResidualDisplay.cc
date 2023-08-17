
#include "TPCManager.cc"
int runnum = 5721;
TString dir = base_dir;
int layer,row;
double rx,rz,rr,rw,rl;
double sx,sz,sr,sw,sl;
TFile* file;
TTree* tree;
void TPCResidualDisplay(){
	gStyle->SetOptStat(0);
	gTPCManager.InitializeHistograms();
	TString filename = Form("HSBeamthroughHistsFit.root",runnum);
//	dir+="Cor1st/";
	dir+="NoCor/";
//	dir+="NoCorL29/";
	file = new TFile(dir+filename);
	tree = (TTree*)file->Get("tree");
	tree->SetBranchAddress("layer",&layer);
	tree->SetBranchAddress("row",&row);
	tree->SetBranchAddress("residual_xCor",&rx);
	tree->SetBranchAddress("residual_zCor",&rz);
	tree->SetBranchAddress("residual_wCor",&rw);
	tree->SetBranchAddress("residual_l",&rl);
	tree->SetBranchAddress("residual_r",&rr);
	tree->SetBranchAddress("resolution_x",&sx);
	tree->SetBranchAddress("resolution_z",&sz);
	tree->SetBranchAddress("resolution_r",&sr);	
	tree->SetBranchAddress("resolution_w",&sw);	
	tree->SetBranchAddress("resolution_l",&sl);	
}
void DrawTPCResidualX(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResidualX");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			gTPCManager.SetPadContent(layer,row,-2);
			continue;
		}
//		if(abs(rx)>3) continue;
		gTPCManager.SetPadContent(layer,row,rx);
	}
	h->Draw("colz");
}
void DrawTPCResidualZ(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResidualZ");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			gTPCManager.SetPadContent(layer,row,-2);
			continue;
		}
//		if(abs(rx)>3) continue;
		gTPCManager.SetPadContent(layer,row,rz);
	}
	h->Draw("colz");
}
void DrawTPCResidualR(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResidualR");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			gTPCManager.SetPadContent(layer,row,-2);
			continue;
		}
		gTPCManager.SetPadContent(layer,row,rr);
	}
	h->Draw("colz");
}
void DrawTPCResidualW(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResidualW");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			gTPCManager.SetPadContent(layer,row,-2);
			continue;
		}
		gTPCManager.SetPadContent(layer,row,rw);
	}
	h->Draw("colz");
}
void DrawTPCResidualL(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResidualL");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			gTPCManager.SetPadContent(layer,row,-2);
			continue;
		}
		gTPCManager.SetPadContent(layer,row,rl);
	}
	h->Draw("colz");
}
void DrawTPCResolutionX(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResolutionX");
	h->Reset("");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			continue;
		}
		gTPCManager.SetPadContent(layer,row,sx);
	}
	h->Draw("colz");
}
void DrawTPCResolutionZ(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResolutionZ");
	h->Reset("");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			continue;
		}
		if(sz<0.3) continue;
		gTPCManager.SetPadContent(layer,row,sz);
	}
	h->Draw("colz");
}
void DrawTPCResolutionR(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResolutionR");
	h->Reset("");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			continue;
		}
		if(sr>3)continue;
		gTPCManager.SetPadContent(layer,row,sr);
	}
	h->Draw("colz");
}
void DrawTPCResolutionW(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResolutionW");
	h->Reset("");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			continue;
		}
		if(sr>3)continue;
		gTPCManager.SetPadContent(layer,row,sw);
	}
	h->Draw("colz");
}
void DrawTPCResolutionL(){
	auto h = gTPCManager.GetPadHistogram();
	int ent = tree->GetEntries();
	h->SetTitle("ResolutionL");
	h->Reset("");
	for(int i = 0; i < ent; ++i){
		tree->GetEntry(i);
		if(tpc::Dead(layer,row)){
			continue;
		}
		if(sr>3)continue;
		gTPCManager.SetPadContent(layer,row,sl);
	}
	h->Draw("colz");
}
