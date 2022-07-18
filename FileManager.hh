double XiMinusMass = 1.32171, XiMinusWidth = 1,XiStarMass = 1.530,XiStarWidth = 0.0099;

class FileManager{
	protected:
	protected:
		TFile* DataFile;
		fstream ParameterFile;
		fstream OldParameterFile;
		TChain* DataChain;
		typedef vector<double> ParamArray;
		typedef map<TString,ParamArray> ParamMap;
		typedef ParamMap::const_iterator PIterator;	
		ParamMap Map;

	public:
		FileManager(){}
		virtual void LoadFile(TString FileName){DataFile = new TFile(FileName,"READ");}
		virtual void LoadChain(TString ChainName){DataChain = (TChain*)DataFile->Get(ChainName);}
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



void FileManager::WriteParameter(vector<int> ID,vector<double> Param){
	for(int i=0;i<ID.size();++i){
		ParameterFile<<ID[i]<<"\t";
	}
	for(int i=0;i<Param.size()-1;++i){
		ParameterFile<<Param[i]<<"\t";
	}
	ParameterFile<<Param[Param.size()-1]<<endl;
}

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

