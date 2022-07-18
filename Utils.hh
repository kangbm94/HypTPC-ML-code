#include "Math.hh"
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
double SlewFunc(double x, double p1,double p2,double p3){
	return p1/sqrt(x-p2)-p3;
}
double SlewFunc(double* x, double* p){
	return SlewFunc(x[0],p[0],p[1],p[2]);
}
double SlewInv(double x, double p1,double p2,double p3){
	return p2-p1*p1/((x-p3)*(x-p3));
}
double SlewInv(double* x, double* p){
	return SlewInv(x[0],p[0],p[1],p[2]);
}
void PHCWritter(fstream& f, int* ids,double p1,double p2, double p3){
	//ids: CId, PlId, Seg, UorD, Type, nParam
	f<<ids[0]<<"\t"	;
	f<<ids[1]<<"\t"	;
	f<<ids[2]<<"\t"	;
	f<<ids[3]<<"\t"	;
	f<<ids[4]<<"\t"	;
	f<<ids[5]<<"\t" ;
	f<<p1<<"\t"	;
	f<<p2<<"\t";	
	f<<p3<<endl;
}
TString Filename(int run){
	if(run<10000){
		return Form("../k18-analyzer/rootfiles/Hodoscope0%d.root",run);
	}
	else{
		return Form("../k18-analyzer/rootfiles/Hodoscope%d.root",run);
	}
}

TString Detname[10]={"","bh1","bh2","bac","htof","bvh","tof","","wc","wcsum"};
int Segment[10]={0,0,0,0,0,0,0,0,0,0};
TString ToF(int a, int b){
	return Form("(%smt[%d][0]-%smt[%d][0])",Detname[a].Data(),Segment[a],Detname[b].Data(),Segment[b]);
}
TString Time(int det,int conf){
	if(det!=2){
		if(conf==0){
			return Form("(%sutime[%d][0])",Detname[det].Data(),Segment[det]-1);
		}
		else if (conf==1){
			return Form("(%sdtime[%d][0])",Detname[det].Data(),Segment[det]-1);
		}
		else{
			return Form("(%smt[%d][0])",Detname[det].Data(),Segment[det]-1);
		}
	}
	else{
		if(conf==0){
			return Form("(%sutime[%d][0]-0.5*(bh1utime+bh1dtime))",Detname[det].Data(),Segment[det]-1);
		}
		else if (conf==1){
			return Form("(%sdtime[%d][0]-0.5*(bh1utime+bh1dtime))",Detname[det].Data(),Segment[det]-1);
		}
		else{
			return Form("(%smt[%d][0])",Detname[det].Data(),Segment[det]-1);
		}
	
	}
}
TCut TimeCut(int det,double t1,double t2){
	TCut tcut1 = Form("%s>%f",Time(det,3).Data(),t1);
	TCut tcut2 = Form("%s<%f",Time(det,3).Data(),t2);
	if(det!=2){
		return tcut1&&tcut2;
	}
	else{
		TCut bh1nh= "bh1nhits==1";
		TCut bh2segc = Form("bh2hitpat[0]==%d",Segment[det]);
		return tcut1&&tcut2&&bh1nh&&bh2segc;
	}
}
TString UorD(int conf){
	if(conf==0){
		return "u";
	}
	else{
		return "d";
	}
}
TString DE(int a,int conf){//0,1
	TString val;
	if(conf==0){
		val=Form("(%sude[%d])",Detname[a].Data(),Segment[a]-1);	
	}
	else if(conf==1){
		val=Form("(%sdde[%d])",Detname[a].Data(),Segment[a]-1);	
	}
	else {
		cout<<"FUCK YOU!"<<endl;
		val = "FUCK";
	}
	return val;
}
TString BF2D(int det1, TString ht, int conf){
	TString val;
	val=Time(det1,conf)+":"+DE(det1,conf)+">>"+ht;
	return val;
}
TString AF2D(int det1,  TString ht,double p0,double p1,double p2,int conf){
	TString val;
	val=Time(det1,conf)+Form("-(%f/sqrt(%s-%f))+%f",p0,DE(det1,conf).Data(),p1,p2)+":"+DE(det1,conf)+">>"+ht;
	return val;
}
int RebinHist(TH1* hist){
	int peak=hist->GetBinContent(hist->GetMaximumBin());
	if(peak<51){
		int rebin_size=100/peak;
		if(rebin_size==7){
			return 6;
		}
		else if(rebin_size==11){
			return 10;
		}
		else
		{
			return rebin_size;
		}
	}
	else{
		return 1;
	}
}
