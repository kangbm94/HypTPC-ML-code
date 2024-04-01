#ifndef XiFitter_h
#define XiFitter_h
#include "KinFit.hh"
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
class XiFitter: virtual public KinematicFitter{
	//R -> P + Q;
	protected:	
		double mP = 0.938272, mPi = 0.139570, mLd = 1.115683, mXi = 1.32171;
		TLorentzVector P;
		TVector3 Pres;
		TLorentzVector PCor;

		TLorentzVector Pi1;
		TVector3 Pi1res;
		TLorentzVector Pi1Cor;

		TLorentzVector Pi2;
		TVector3 Pi2res;
		TLorentzVector Pi2Cor;
		vector<double> MassDiffs;
		bool MeasDir = false;

		TLorentzVector XiCor;
	public:
		XiFitter(){}
		XiFitter(TLorentzVector P_,TLorentzVector Pi1,TLorentzVector Pi2);
		void SetInvMass(double IM){
			mR = IM;
		}
		vector<TLorentzVector> GetDaughterLV(){
			vector<TLorentzVector> ret = {PCor,Pi1Cor,Pi2Cor};
			return ret;
		}
		TLorentzVector GetLambda(){
			return PCor+Pi1Cor;
		}
		TLorentzVector GetXi(){
			return PCor+Pi1Cor+Pi2Cor;
		}
		void ToDecayPlane();
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif
