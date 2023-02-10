#include "TPCPadHelper.hh"
#ifndef BeamRemover_h
#define BeamRemover_h
class BeamRemover{
	private:
		vector<TVector3> rawhitarray;
		vector<vector<TVector3>> hitarray;
		vector<vector<int>> hough_flag;
		int MinBeamHit = 15;
		double MinZCut = -180;
		double MaxZCut = 150;
		std::vector<int> m_layer_info;
		std::vector<bool> m_Accidental;
		
		vector<double> cx;
		vector<double> cy;
		vector<double> z0;
		vector<double> r;
		vector<double> dz;
		
		vector<vector<double>> h_cx;
		vector<vector<double>> h_cy;
		vector<vector<double>> h_z0;
		vector<vector<double>> h_r;
		vector<vector<double>> h_dz;
		

		double y_min,y_max;
		double x_min,x_max;
		double Ywidth = 10; 
		double Ycut;
		double Xcut;
		bool enable;
	
		TH3D* Ci_hist;
		TH2D* histY;
		TH1D* Ydistrib;
		const Double_t ZTarget = -143.; // Target from center
  	const Int_t nBin_rdiff = 220;
 		const Double_t rdiff_min = -110.;
  	const Double_t rdiff_max = 110.;
  	const Int_t nBin_theta = 180;
  	const Double_t theta_min = (1.5)*cos(-1);
  	const Double_t theta_max = (2.5)*acos(-1);//Charge < -1 region.
  	const Int_t nBin_p = 100;
  	const Double_t p_min = 1600.;//MeV/c
  	const Double_t p_max = 2000.;//MeV/c
  	const int    thetaY_ndiv = 1080;
  	const double thetaY_min  =   60.;
  	const double thetaY_max  = 120.;
  	const int    r_ndiv =  2000;
  	const double r_min  = -5000.;
  	const double r_max  =  5000.;
  	std::vector<Double_t> hough_x;
  	std::vector<Double_t> hough_y;
  	std::vector<Double_t> hough_z;
  	std::vector<Double_t> hough_rd;
  	std::vector<Double_t> hough_theta;
  	std::vector<Double_t> hough_p;
		int MaxNBeam = 3;
		double MaxHoughWindow = 15;
		double MaxHoughWindowY = 5;
		const double Const = 0.299792458;
		double dMagneticField = 0.95;
		TH2D* ZXHist;
		TGraph* ArcGraph;
		TGraph* YThetaGraph;
		std::vector<TEllipse*> circles;
		TF1* Arc;
	public:
		BeamRemover(double xmin,double xmax,double ymin,double ymax);
		~BeamRemover();
		void LoadHit(TVector3 hit);
		TH2D* GetHoughHist(){
			return ZXHist;
		}
		vector<TEllipse*> GetCircle(int i);
		int SearchPeaks(TH1D* hist, std::vector<double> &peaks);
		void AssignHits(std::vector<double> peaks);
		void SetNPeaks(int np){
			hitarray.resize(np); } void AddHit(int i, TVector3 hit){
			hitarray[i].push_back(hit);
		}
		void DoCircleHough(int i);
		void DoYThetaHough(int i);
		double GetYwidth(){
			return Ywidth;
		}
		TH1D* GetYDistrib(){
			return Ydistrib;
		}
		int CompareHough(TVector3 pos, vector<double> hcx,vector<double>hcy,vector<double> hr); 
};
#endif
