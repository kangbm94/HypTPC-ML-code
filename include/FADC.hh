#include "TPCManager.hh"
static int nhist=0;;
class FADC{
	protected:
		int layer=-1,row=-1;
		double ped=0,mean=0,sig=0;
		double rawped=0,rawmean=0,rawsig=0;
		TH1D* hist=nullptr;
		TH1D* rawhist=nullptr;
		TF1* func=nullptr;
		double t1 = 0,t2=170;
		int npeaks=0;
		vector<double>peaks;
	public:
		FADC(){}
		~FADC(){
//					delete hist;
//					delete rawhist;
		}
		FADC(vector<double>tb,vector<double>fadc,vector<double>rawfadc,int l,int r);
		int GetLayer(){return layer;}
		int GetRow(){return row;}
		TH1D* GetHistogram(){
			return hist;
		}
		TH1D* GetRawHistogram(){
			return rawhist;
		}
		int DoFit(double* par);
		double GetRMS(){return sig;}
};
class Baseline: public FADC{
	public:
		Baseline(vector<double>tb,vector<double>fadc);
};

FADC::FADC(vector<double>tb,vector<double>fadc,vector<double>rawfadc,int l,int r){
	
	int nb = tb.size();
	for(int i=0;i<nb;++i){
		if(i<20) ped+=fadc[i]/20;
		if(i<20) rawped+=rawfadc[i]/20;
		mean+=fadc[i]/nb;
		rawmean+=rawfadc[i]/nb;
		sig+=fadc[i]*fadc[i]/nb;
		rawsig+=rawfadc[i]*rawfadc[i]/nb;
	}
	sig=sqrt(sig-mean*mean);
	rawsig=sqrt(rawsig-rawmean*rawmean);
	mean-=ped;
	rawmean-=rawped;
	layer = l;row=r;
	hist=nullptr;
	hist = new TH1D(Form("hist%d",nhist),Form("Waveform(%d,%d)_%d",layer,row,nhist),nb,t1,t2);
	rawhist=nullptr;
	rawhist = new TH1D(Form("rawhist%d",nhist),Form("RawWaveform(%d,%d)_%d",layer,row,nhist),nb,t1,t2);
	for(int i=0;i<nb;++i){
		hist->SetBinContent(tb[i],fadc[i]-ped);
		rawhist->SetBinContent(tb[i],rawfadc[i]-rawped);
	}
	hist->GetYaxis()->SetRangeUser(-200,700);
	rawhist->GetYaxis()->SetRangeUser(-200,700);
	TSpectrum spec(20);
	double sigma=3,threshold=0.4;
	npeaks = spec.Search(hist,sigma,"goff",threshold); 
	//npeaks = spec.Search(hist,sigma,"",threshold); 
	double* x_peaks = spec.GetPositionX();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(x_peaks[i]);
	}
	stable_sort(peaks.begin(),peaks.end());
};
Baseline::Baseline(vector<double>tb,vector<double>fadc){
	
	int nb = tb.size();
	for(int i=0;i<20;++i){
		if(i<20) ped+=fadc[i]/20;
	}
	for(int i=0;i<nb;++i){
		mean+=fadc[i]/nb;
		sig+=fadc[i]*fadc[i]/nb;
	}
	sig=sqrt(sig-mean*mean);
	mean-=ped;
	hist=nullptr;
	hist = new TH1D(Form("hist%d",nhist),Form("Baseline_%d",nhist),nb,t1,t2);
	for(int i=0;i<nb;++i){
		hist->SetBinContent(tb[i],fadc[i]-ped);
	}
	hist->GetYaxis()->SetRangeUser(-200,700);
	TSpectrum spec(20);
	double sigma=3,threshold=0.4;
	npeaks = spec.Search(hist,sigma,"goff",threshold); 
	//npeaks = spec.Search(hist,sigma,"",threshold); 
	double* x_peaks = spec.GetPositionX();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(x_peaks[i]);
	}
	stable_sort(peaks.begin(),peaks.end());
};
int FADC::DoFit(double* par){
	double time[20]={0},de[20]={0};
	TString MultiGumbel;
	for(int i=0;i<npeaks;++i){
		time[i]=peaks[i];
		de[i]=hist->GetBinContent(hist->FindBin(peaks[i]));
		if(i!=0) MultiGumbel+="+";
		MultiGumbel += "["+to_string(i*3)+"]";
		MultiGumbel += "*exp(-(x-["+to_string(i*3+1)+"])"; 
		MultiGumbel += "/["+to_string(i*3+2)+"])";  
		MultiGumbel += "*exp(-exp(-(x-["+to_string(i*3+1)+"])/["+to_string(i*3+2)+"]))";
	}
	MultiGumbel+="+["+to_string(npeaks*3)+"]";
//	cout<<MultiGumbel<<endl;

	func = new TF1("func",MultiGumbel,0,170);
	for(int i=0;i<npeaks;++i){
		func->SetParameter(i*3,de[i]);
		func->SetParameter(i*3+1,time[i]);
		func->SetParLimits(i*3+1,time[i]-10,time[i]+10);
		func->SetParameter(i*3+2,3);
		func->SetParLimits(i*3+2,1,10);
	}
	auto* pm = new TPolyMarker(npeaks,time,de);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.3);
	hist->GetListOfFunctions()->Add(pm);
	func->SetParameter(npeaks*3,0);
	func->SetParLimits(npeaks*3,-0.5*sig,0.5*sig);
	func->SetRange(25,155);
//	cout<<Form("Fitting (%d,%d), npeaks = %d",layer,row,npeaks)<<endl;
	hist->Fit("func","QR");
	for(int i=0;i<10;++i){
		par[i]=0;
	}
	for(int i=0;i<npeaks;++i){
		par[i]=func->GetParameter(3*i);
	}
	return npeaks;
}


class FADCManager:public TPCManager{
	private:
		int evnum;
		vector<int>*layerTpc = new vector<int>;
		vector<int>*rowTpc = new vector<int>;
		vector<int>*fadclayerTpc = new vector<int>;
		vector<int>*fadcrowTpc = new vector<int>;
		vector<vector<double>>*tbTpc = new vector<vector<double>>;
		vector<vector<double>>*fadcTpc = new vector<vector<double>>;
		vector<vector<double>>*rawfadcTpc = new vector<vector<double>>;
		vector<vector<double>>*baselinetbTpc = new vector<vector<double>>;
		vector<vector<double>>*baselineTpc = new vector<vector<double>>;
		vector<FADC> FADCs;
		vector<Baseline> Base;
	public:
		FADCManager(){
		}
		virtual void LoadFile(TString FileName){ 
			DataFile = new TFile(FileName,"READ");
			LoadFADCChain("tpc");
		}
		virtual void SetEvent(int i);
		void LoadFADCChain(TString ChainName);
		bool CheckLayerRowHit(int layer, int row);
		void LoadWaveform();
		void Clear(){FADCs.clear();Base.clear();}
		int GetWaveNum(){
			return fadclayerTpc->size();
//			return FADCs.size();
		}
		int GetBaseNum(){
			cout<<"NBase = "<<Base.size()<<endl;
			return Base.size();
		}
		TH1D* GetWaveformHist(int layer=-1, int row=-1);
		TH1D* GetRawWaveformHist(int layer=-1, int row=-1);
		int DoFit(int layer, int row, double* par);
		TH1D* GetBaseline();
		TH1D* GetRawBaseline();
		bool Check();
};
void FADCManager::SetEvent(int i){
	layerTpc->clear();
	rowTpc->clear();
	fadclayerTpc->clear();
	fadcrowTpc->clear();
	baselineTpc->clear();
	baselinetbTpc->clear();
	DataChain->GetEntry(i);
}
void FADCManager::LoadFADCChain(TString ChainName){
	DataChain = (TChain*) DataFile->Get(ChainName);
	cout<<"entries: "<<DataChain->GetEntries()<<endl;
	DataChain->SetBranchAddress("evnum",&evnum);
	DataChain->SetBranchAddress("layerTpc",&layerTpc);
	DataChain->SetBranchAddress("rowTpc",&rowTpc);
	DataChain->SetBranchAddress("fadclayerTpc",&fadclayerTpc);
	DataChain->SetBranchAddress("fadcrowTpc",&fadcrowTpc);
	DataChain->SetBranchAddress("tbTpc",&tbTpc);
	DataChain->SetBranchAddress("fadcTpc",&fadcTpc);
	DataChain->SetBranchAddress("rawfadcTpc",&rawfadcTpc);
	DataChain->SetBranchAddress("baselinetbTpc",&baselinetbTpc);
	DataChain->SetBranchAddress("baselineTpc",&baselineTpc);
}
bool FADCManager::CheckLayerRowHit(int layer, int row){
	int nh = layerTpc->size();
	bool lay_flg=false,row_flg=false;
	for(int i=0;i<nh;++i){
//		lay_flg=false;row_flg=false;
		int hit_lay = layerTpc->at(i);
		int hit_row = rowTpc->at(i);
		if(layer==hit_lay) lay_flg=true;
		if(row==hit_row) row_flg=true;
		if(lay_flg and row_flg){cout<<Form("Hit (%d,%d)",layer,row)<<endl;break;}
	}
	return lay_flg and row_flg;
}
void FADCManager::LoadWaveform(){
	FADCs.clear();
	Base.clear();
	nhist=evnum;
	int nw = GetWaveNum();
	cout<<"Nwaves = "<<nw<<endl;
	int cnt=0;
	int bcnt=0;
	for(int iw=0;iw<nw;++iw){
		int w_layer = fadclayerTpc->at(iw);
		int w_row = fadcrowTpc->at(iw);
			FADCs.push_back(FADC(tbTpc->at(iw),fadcTpc->at(iw),rawfadcTpc->at(iw),w_layer,w_row));
			int l = FADCs[cnt].GetLayer();
			int r = FADCs[cnt].GetRow();
			cnt++;
		if(CheckLayerRowHit(w_layer,w_row) or 1){			cout<<Form("LR=(%d,%d)",l,r)<<endl;
		}
		/*
		else{
			Base.push_back(FADC(tbTpc->at(i),fadcTpc->at(i),rawfadcTpc->at(i),w_layer,w_row));
			cout<<Form("Baseline=(%d,%d)",w_layer,w_row)<<endl;
			bcnt++;
		}
		*/
	}
	Base.push_back(Baseline(baselinetbTpc->at(0),baselineTpc->at(0)));
}


TH1D* FADCManager::GetWaveformHist(int layer=-1, int row=-1){
	TH1D* h = nullptr;	
	int nw = FADCs.size();
	int nc = 0;
	for(int i=0;i<nw;++i){
		int l_=FADCs[i].GetLayer();
		int r_=FADCs[i].GetRow();
		if(l_==layer and r_==row){
			h=FADCs[i].GetHistogram();
		}
		else{
			cout<<Form("Hist(%d,%d):Wave(%d,%d) Does not exist!",l_,r_,layer,row)<<endl;
		}
	}
	return h;
}
TH1D* FADCManager::GetRawWaveformHist(int layer=-1, int row=-1){
	TH1D* h = nullptr;	
	int nw = FADCs.size();
	int nc = 0;
	for(int i=0;i<nw;++i){
		int l_=FADCs[i].GetLayer();
		int r_=FADCs[i].GetRow();
		if(l_==layer and r_==row){
			h=FADCs[i].GetRawHistogram();
		}
		else{
	//		cout<<Form("RawHist(%d,%d) Does not exist!",l_,r_)<<endl;
		}
	}
	return h;
}
TH1D* FADCManager::GetBaseline(){
	TH1D* h = nullptr;	
	int nw = GetBaseNum();
	if(nw) h=Base[0].GetHistogram();
	return h;
}
TH1D* FADCManager::GetRawBaseline(){
	TH1D* h = nullptr;	
	int nw = GetBaseNum();
	if(nw) h=Base[0].GetRawHistogram();
	return h;
}
int FADCManager::DoFit(int layer,int row,double* par){
	int nw = GetWaveNum();
	int nc = 0;
	int npeaks=0;
	for(int i=0;i<nw;++i){
		int l_=FADCs[i].GetLayer();
		int r_=FADCs[i].GetRow();
		if(l_==layer&&r_==row){
			npeaks= FADCs[i].DoFit(par);
			//			nc++;
		}
	}
	return npeaks;
}
bool FADCManager::Check(){
	int nb = Base.size();
	double par[10];
	for(int i=0;i<nb;++i){
		Base[i].DoFit(par);
		for(int j=0;j<10;j++){
			if(par[j]>100*exp(1)) return false;
		}
	}
	return true;
}
