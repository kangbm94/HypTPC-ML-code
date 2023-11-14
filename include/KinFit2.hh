#ifndef KinFit_h
#define KinFit_h
namespace{
}
class KinematicFitter{
	protected:
		// R -> P + Q;
		int step = 0;
		int best_step = 0;
		int nConst;
		int nMeas;
		int nUnkn;
		double Chi2_cut = 0.1;
		double Best_Chi2 = 1e9;
		double Best_MassDiff = 1e9;
		int MaxStep = 100;


		
		vector<double*> Pulls;
		vector<double> Lambdas;
		vector<double> Chi2s;
		vector<double> MassDiffs;
		vector<TMatrixD> Measurements;
		vector<TMatrixD> Unknowns;
		vector<TMatrixD> Variancies;//Actually, the inverse of the variance
		

		vector<TMatrixD> FMats;
		vector<TMatrixD> dFdMs;//\pdf{Constraint}{Measurement}
		vector<TMatrixD> dFdUs;//\pdf{Constraint }{Unknown }
		vector<TMatrixD> rMats;
		vector<TMatrixD> sMats;
	public:
		KinematicFitter(){};
		//Setters
		void SetVariance(double* var);
		void SetChi2DifCut(double cut){
			Chi2_cut = cut;
		}
		void SetMaximumStep(int max){
			MaxStep = max;
		}
		//Getters
		double GetLambda(int ent = -1){
			if(ent == -1) ent = step;
			return Lambdas.at(ent);
		}
		void Clear();
		virtual double DoKinematicFit(bool Do);
		int GetNStep(){
			return step;
		}
		int GetBestStep(){
			return best_step;
		}
		vector<double>GetStepChi2(){
			return Chi2s;
		}
		vector<double>GetStepMassDiff(){
			return MassDiffs;
		}
		int GetNDF(){
			return nConst - nUnkn;	
		}
	protected:
		//User Part: Set variables and constraints as you want.
		virtual void Initialize(){};
		virtual void SampleStepPoint(){};
		virtual void SetConstraints(){};
		virtual void Finalize(){};	
		TMatrixD TransposeMatrix(TMatrixD M);		

		//Core. Do not modify unless necessary
		void ProcessStep();
};
#endif
