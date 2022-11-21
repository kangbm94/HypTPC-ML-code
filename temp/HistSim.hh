#include "../include/TPCPadHelper.hh"
class TestClass{
	private:
		TH2Poly* hist;
	public:
		TestClass(TH2Poly* in){
			hist = in;
		}
		void Draw(){
			hist->Draw();
		}
		void FillHist(int i){
			hist->SetBinContent(i,10);
		}
};
