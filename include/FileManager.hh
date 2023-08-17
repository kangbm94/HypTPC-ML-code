#ifndef FileManager_h
#define FileManager_h
double XiMinusMass = 1.32171, XiMinusWidth = 1,XiStarMass = 1.530,XiStarWidth = 0.0099;

class FileManager{
	private:
	protected:
		TFile* DataFile;
		fstream ParameterFile;
		fstream OldParameterFile;
		TChain* DataChain;
		typedef vector<double> ParamArray;
		typedef map<TString,ParamArray> ParamMap;
		typedef ParamMap::const_iterator PIterator;	
		ParamMap Map;
		TFile* OutFile;	
		TTree* OutChain;
		double OutBufD[100];
		int OutBufI[100];
		vector<double>* OutBufV[20];
	public:
		FileManager(){}
		void WriteTag(TString title, TString comment);
		void LoadFile(TString FileName){
			if(DataFile) delete DataFile;
			DataFile = new TFile(FileName,"READ");}
		void LoadChain(TString ChainName){DataChain = (TChain*)DataFile->Get(ChainName);}
		int GetEntries(){
			return DataChain->GetEntries();
		}

		void LoadParameterFile(TString FileName){OldParameterFile.open(FileName,fstream::in);}
		TObject* GetHistogram(int HistNumber);
		TObject* DrawHistogram(TString Argument, TCut Cut, TString Options = "");
		void MakeParameterFile(TString FileName){ParameterFile.open(FileName,fstream::out);}
		void MakeParameterMap(TString FileName);
		void WriteParameter(vector<int> ID,vector<double> Param);
		void WriteComment(TString Comment){ParameterFile<<Comment<<endl;}
		TChain* MakePublicChain(){return DataChain;}
		TChain* GetPublicChain(){return DataChain;}
		void LoadParamMap(TString FileName);
		double GetParameter(TString key, int i=0);
	
		void CreateFile(TString filename){
			OutFile = new TFile(filename,"recreate");
			OutChain = new TTree("tree","tree");
		}
		void WriteFile(){
			OutFile ->Write();
		}
		void ClearBuffer(){
			for(int i=0;i<20;++i) if(OutBufV[i]) OutBufV[i]->clear();
		}
		void OutBranch(TString Title, int i, int conf ){
			if(conf==1)OutChain->Branch(Title,&OutBufD[i]);	
			else if(conf ==0) OutChain->Branch(Title,&OutBufI[i]);	
			else if(conf ==2) {
				OutBufV[i] = new vector<double>;
				OutChain->Branch(Title,&OutBufV[i]);	
			}
		}
		
		TTree* GetOutChain(){return OutChain;}

};

TObject* FileManager::GetHistogram(int HistNumber){
	TString HistName = Form("h%d",HistNumber);
	return (TObject*)DataFile->Get(HistName);
}

TObject* FileManager::DrawHistogram(TString Argument, TCut Cut, TString Options = ""){
//	TPad* DumPad;
//	DumPad->cd();
//	DataChain->Draw(Argument, Cut, Options+"goff");
	DataChain->Draw(Argument, Cut, Options);
	TH1* h = (TH1*)gPad->GetPrimitive("h");
	cout<<h->GetEntries()<<endl;
	return h;
}



void FileManager::WriteParameter(vector<int> ID={},vector<double> Param={}){
	for(int i=0;i<ID.size();++i){
		ParameterFile<<ID[i]<<"\t";
	}
	for(int i=0;i<Param.size()-1;++i){
		ParameterFile<<Param[i]<<"\t";
	}
	ParameterFile<<Param[Param.size()-1]<<endl;
}

void FileManager::WriteTag(TString title, TString comment){
	auto tag = new TNamed(title.Data(),comment.Data());
	tag->Write();
};
/*
void FileManager::LoadParamMap(TString FileName){
	ifstream ParamFile(FileName);
	if(!ParamFile.is_open()){
		cout<<"No Such File : "<<FileName<<endl;
	}
	
	TString Buf = "\n";
	TString Line;
	while(ParamFile.good() && Line.ReadLine(ParamFile)){
		Buf+=Line+"\n";
		if(Line.IsNull()||Line[0]=='#') continue;
		istringstream input_line(Line.Data());
		TString key;
		input_line>>key;
		cout<<"Key is : "<<key<<endl;
		ParamArray param_array;
		double param;
		while(input_line>>param){
			param_array.push_back(param);
		}
		Map[key] = param_array;		
	}

}
*/
#endif
