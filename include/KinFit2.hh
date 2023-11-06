#ifndef KinFit_h
#define KinFit_h
namespace{
}
class KinematicFitter{
	private:
		// R -> P + Q;
		int step = 0;
		int best_step = 0;
		int nConst;
		int nMeas;
		int nUnkn;
		double Chi2_cut = 0.1;
		double Best_Chi2 = 1e9;
		int MaxStep = 100;
		bool MeasDir = false;


		TLorentzVector P;
		TVector3 Pres;
		TLorentzVector PCor;
		double mP;

		TLorentzVector Q;
		TVector3 Qres;
		TLorentzVector QCor;
		double mQ;
		
		TLorentzVector R;
		TVector3 Rres;
		TLorentzVector RCor;
		double mR;
		

		vector<double> Lambdas;
		vector<double> Chi2s;
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
		KinematicFitter(TLorentzVector P_,TLorentzVector Q_,TLorentzVector R_);
		//Setters
		void SetVariance(double* var);
		void SetChi2DifCut(double cut){
			Chi2_cut = cut;
		}
		void SetMaximumStep(int max){
			MaxStep = max;
		}
		void SetInvMass(double IM){
			mR = IM;
		}
		//Getters
		double GetLambda(int ent = -1){
			if(ent == -1) ent = step;
			return Lambdas.at(ent);
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {PCor,QCor,RCor};
			return ret;
		}
		void UseVertex(bool status,TVector3 Vert1,TVector3 Vert2);
		void Clear();
		double DoKinematicFit();
	private:
		//User Part: Set variables and constraints as you want.
		void Initialize();
		void SetConstraints();
		void Finalize();	
		TMatrixD TransposeMatrix(TMatrixD M);		

		//Core. Do not modify unless necessary
		void ProcessStep();
};
#endif
