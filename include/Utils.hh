#include "Math.hh"
#ifndef Utils_h
#define Utils_h
const double K18HS = 1318.9;
TString base_dir = "~/Desktop/JPARC/E42/MayRun/rootfiles/CH2/TPC/";
void Indicator(int i, int ent){
	if(i>ent){
		double tmp = ent;
		i=ent;
		tmp=i;
	}
	int tot =ent/100+1;
	bool stat= i%tot;
	if(!stat){
		cout<<setw(2)<<i/tot<<" \%"<<endl;
	}
	if(i==ent-1){
		cout<<"Loop Done"<<endl;
	}
}
void Indicator(int i, int ent, TString Loopname){
	int tot =ent/100+1;
	bool stat= i%tot;
	if(!stat){
		cout<<Loopname+": "<<setw(2)<<i/tot<<" \%"<<endl;
	}
	if(i==ent-1){
		cout<<Loopname<< ": 100%"<<endl;
	}
}
void PressAnyKey(){	
	gSystem->ProcessEvents();
	cout<<"Press Enter to continue"<<endl;
	cin.ignore();
}
void FillHist(TString title, TTree* tree,TString branch,TCut cut){
	tree->Draw(branch+">>"+title,cut);
}
void FillHistColz(TString title, TTree* tree,TString branch,TCut cut){
	tree->Draw(branch+">>"+title,cut,"colz");
}
bool ReadConfLine(ifstream& file, double* data){
	if(!file.is_open()){
		cout<<"file not open"<<endl;
		return false;
	}
	TString line;
	if(file.good()&&line.ReadLine(file)){
		line.ReplaceAll(",","");
		line.ReplaceAll("\"","");
		std::istringstream iss(line.Data());
		std::istream_iterator<std::string> begin(iss);
		std::istream_iterator<std::string> end;
		std::vector<TString> v(begin, end);
		int nl = v.size();
		for(int i=0;i<nl;++i){
			data[i]=v[i].Atoi();
		}
		return true;
	}
	else{
		return true;
	}
}


void ReadCSV(fstream& file, double* data){//
	string line,token;
	stringstream iss;
	int cnt=0;
	getline(file,line);
	iss<<line;
	while(getline(iss,token,',')){
		data[cnt]=stod(token);
		cnt++;
	}
}
void ReadTSV(fstream& file, double* data){//
	string line,token;
	stringstream iss;
	int cnt=0;

	getline(file,line);
	iss<<line;
	//	cout<<line<<endl;
	while(getline(iss,token,'\t')){
		data[cnt]=stod(token);
		cnt++;
	}
}

static int niter=0;
static int canv=0;
double GetPeakPosition(TH1* h){
	double peak=0;
	int nb = h->GetNbinsX();
	double bw = h->GetBinWidth(1);
	double x0 = h->GetBinCenter(1)-bw;
	double x1 = h->GetBinCenter(nb)+bw;
	if(h->GetEntries()>0){
		int mb = h->GetMaximum();
		peak = h->GetBinCenter(mb);
	}
	if(peak<x0||peak>x1){
		peak=0;
	}
	return peak;
}
int Bin(int nbin,double r1,double r2,double r){
	double dr = (r2-r1)/(nbin-1);
	r-=(r1-dr/2);
	int val =  (r+dr*1e6)/dr;
	return val-1e6;
}
double BinPos(int nbin, double r1,double r2, int bin){
	if(nbin==1){
		return (r1+r2)/2;
	}
	double dr = r2-r1;
	double bw = dr/(nbin-1);
	return r1+bin*bw;
}
#endif
