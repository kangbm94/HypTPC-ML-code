#include "../include/Hough.hh"
#ifndef Hough_cc
#define Hough_cc
namespace{
	const auto qnan = TMath::QuietNaN();
	const double& valueHSCalc = 0.90607;
	const double& valueHSHall = 0.8738; 
	const double Const = 0.299792458;

	// HoughXZ binning
	//static const Int_t    Li_theta_ndiv = 180;
	static const Int_t    Li_theta_ndiv = 360;
	static const Double_t Li_theta_min  =   0;
	static const Double_t Li_theta_max  = 180;
	//static const Int_t    Li_rho_ndiv =  180;
	static const Int_t    Li_rho_ndiv =  720;
	static const Double_t Li_rho_min  = -720;
	static const Double_t Li_rho_max  =  720;

	// HoughCircle binning
	//  const Int_t    nBin_rdiff = 22;
	const Int_t    nBin_rdiff = 55;
	const Double_t rdiff_min = -110.;
	const Double_t rdiff_max = 110.;
//	const Int_t    nBin_theta = 1440;
	const Int_t    nBin_theta = 360;
	const Double_t theta_min = -1.*acos(-1);
	const Double_t theta_max = acos(-1);
	const Int_t    nBin_p = 900;
//	const Int_t    nBin_p = 300;
	const Double_t pmin = 50.;
	const Double_t pmax = 1550.; //MeV/c
}									 //const Double_t pmax = 2150.; //MeV/c
TH2D *histLinear = new TH2D("histLinear", "; #theta (deg.); #rho (mm)",
		Li_theta_ndiv, Li_theta_min, Li_theta_max,
		Li_rho_ndiv, Li_rho_min, Li_rho_max);
TH3D *histCircle = new TH3D("histCircle",";rd (mm); theta (rad); p(MeV/c)",
		nBin_rdiff, rdiff_min,  rdiff_max,
		nBin_theta, theta_min, theta_max,
		nBin_p, pmin, pmax);
TH3D *hc2 = new TH3D("hc2","f;rd (mm); theta (rad); p(MeV/c)",
		nBin_rdiff, rdiff_min,  rdiff_max,
		nBin_theta, theta_min, theta_max,
		nBin_p, pmin, pmax);
TH2D* hist_temp ;
TH2D* hist_temp2 ;
TH2D* hist_temp3 ;
//TH2D* hist_space = new TH2D("histSpace","histSpace",300,-150,400,300,-300,300);
TGraph* gr = new TGraph();
TGraph* gr2 = new TGraph();
// HoughY binning
//const Ynt_t    thetaY_ndiv = 360; //previous value
const Int_t    thetaY_ndiv = 1080;
const Double_t thetaY_min  =   0.;
const Double_t thetaY_max  = 180.;
const Int_t    r_ndiv =  5000;
const Double_t r_min  = -5000.;
const Double_t r_max  =  5000.;
TH2D *histY = new TH2D("histY",";theta (deg.);r (mm)", thetaY_ndiv, thetaY_min, thetaY_max, r_ndiv, r_min, r_max);

	Bool_t
HoughTransformLineXZ(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		Double_t *LinearPar, Int_t MinNumOfHits /*=8*/)
{


	Bool_t status = true;

	//Hough space
	//for TPC linear track
	//rho = z*cos(theta) + x*sin(theta)
	//x = p[0] + p[2]*z
	histLinear->Reset();

	for(Int_t i=0; i<gHitPos.size(); ++i){
		for(int ti=0; ti<thetaY_ndiv; ti++){
			Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
			Double_t rho = TMath::Cos(theta*TMath::DegToRad())*(gHitPos[i].Z() - ZTarget)
				+ TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].X();
			histLinear->Fill(theta, rho);
		} //ti
	} //i


	Int_t maxbin = histLinear->GetMaximumBin();
	Int_t mx, my, mz;
	histLinear->GetBinXYZ(maxbin, mx, my, mz);
	MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

	Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
	Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);

	LinearPar[0] = mr/TMath::Sin(mtheta); //m_x
	LinearPar[2]  = -TMath::Cos(mtheta)/TMath::Sin(mtheta); //m_u


	if(histLinear->GetMaximum() < 0.5*MinNumOfHits) status = false;
	return status;

}

	void
HoughTransformLineYZ(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		Double_t *LinearPar, Double_t MaxHoughWindowY)
{


	//Hough space
	//for TPC linear track
	//rho = z*cos(theta) + y*sin(theta)
	//y = p[1] + p[3]*z
	histLinear->Reset();
	for(Int_t i=0; i<gHitPos.size(); ++i){
		Double_t dist = TMath::Abs(LinearPar[2]*(gHitPos[i].Z() - ZTarget) - gHitPos[i].X() +
				LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
		if(dist < MaxHoughWindowY){
			for(int ti=0; ti<Li_theta_ndiv; ti++){
				Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
				Double_t rho = TMath::Cos(theta*TMath::DegToRad())*(gHitPos[i].Z() - ZTarget)
					+ TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].Y();
				histLinear->Fill(theta, rho);
			} //ti
		} //dist
	} //i


	Int_t maxbin = histLinear->GetMaximumBin();
	Int_t mx, my, mz;
	histLinear->GetBinXYZ(maxbin, mx, my, mz);
	MaxBin[0] = mx; MaxBin[1] = my;

	Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
	Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);

	LinearPar[1] = mr/TMath::Sin(mtheta); //m_y
	LinearPar[3] = -TMath::Cos(mtheta)/TMath::Sin(mtheta); //m_v


}
	void
HoughTransformLineYX(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		Double_t *LinearPar, Double_t MaxHoughWindowY)
{


	//Hough space
	//for TPC linear track
	//rho = x*cos(theta) + y*sin(theta)
	//y = tmp0 + tmp1*x
	histLinear->Reset();
	for(Int_t i=0; i<gHitPos.size(); ++i){
		Double_t dist = TMath::Abs(LinearPar[2]*(gHitPos[i].Z() - ZTarget) - gHitPos[i].X() +
				LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
		if(dist < MaxHoughWindowY){
			for(int ti=0; ti<Li_theta_ndiv; ti++){
				Double_t theta = histLinear->GetXaxis()->GetBinCenter(ti+1);
				Double_t rho = TMath::Cos(theta*TMath::DegToRad())*gHitPos[i].X()
					+ TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].Y();
				histLinear->Fill(theta, rho);
			} //ti
		} //dist
	} //i


	Int_t maxbin = histLinear->GetMaximumBin();
	Int_t mx, my, mz;
	histLinear->GetBinXYZ(maxbin, mx, my, mz);
	MaxBin[0] = mx; MaxBin[1] = my;

	Double_t mtheta = histLinear->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
	Double_t mr = histLinear->GetYaxis()->GetBinCenter(my);
	Double_t tmp0 = mr/TMath::Sin(mtheta);
	Double_t tmp1 = -TMath::Cos(mtheta)/TMath::Sin(mtheta);
	LinearPar[1] = tmp0 + LinearPar[0]*tmp1;
	LinearPar[3] = LinearPar[2]*tmp1;


}


	void
HoughTransformLineYTheta(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		Double_t *HelixPar, Double_t MaxHoughWindowY)
{

	delete gr2;
	gr2 = new TGraph();
	gr2->SetMarkerStyle(22);
	gr2->SetMarkerSize(2);

	//Hough space
	//mr = tmp_xval*cos(mtheta) + tmpz*sin(mtheta)
	//tmpz = p[2] + p[4]*tmp_xval
	histY -> Reset();
	std::vector<double>temp_tvec;
	int npts=0;
	bool crossing = false;
	bool negative = false;
	bool positive = false;
	for(Int_t i=0; i<gHitPos.size(); ++i){
		Double_t tmpx = -gHitPos[i].X();
		Double_t tmpy = gHitPos[i].Z() - ZTarget;
		Double_t tmpz = gHitPos[i].Y();
		Double_t r_cal = sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
		Double_t dist = fabs(r_cal - HelixPar[3]);
		if(dist < MaxHoughWindowY){
			Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
			temp_tvec.push_back(tmpt);
			npts++;
			if(tmpt > acos(-1) * 3/4) positive = true;
			if(tmpt < acos(-1) *-3/4) negative = true;
		}
	}
	if(positive and negative) crossing = true;
	std::sort(temp_tvec.begin(),temp_tvec.end());
	int nts = temp_tvec.size();
	double start = 0.;// 0 = -Pi, 1 = 2 Pi
	double end = 2./3.;
	if(crossing){
		start = 1./3.;
		end = 1;
	}
	for(Int_t i=0; i<gHitPos.size(); ++i){
		Double_t tmpx = -gHitPos[i].X();
		Double_t tmpy = gHitPos[i].Z() - ZTarget;
		Double_t tmpz = gHitPos[i].Y();
		Double_t r_cal = sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
		Double_t dist = fabs(r_cal - HelixPar[3]);
		if(dist < MaxHoughWindowY){
				Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
//				if(tmpt< 0 and crossing)tmpt += 2*acos(-1); 
			for(int ti=0; ti<thetaY_ndiv; ti++){
				Double_t theta = histY->GetXaxis()->GetBinCenter(ti+1);
				Double_t tmp_xval = HelixPar[3]*tmpt;
				histY->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval
						+sin(theta*acos(-1.)/180.)*tmpz);
			} //ti
			gr2->AddPoint(HelixPar[3]*tmpt,tmpz);
		} //dist
	}


	Int_t maxbin = histY->GetMaximumBin();
	Int_t mx, my, mz;
	histY->GetBinXYZ(maxbin, mx, my, mz);
	MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

	double mtheta = histY->GetXaxis()->GetBinCenter(mx)*acos(-1.)/180.;
	Double_t mr = histY->GetYaxis()->GetBinCenter(my);
	HelixPar[2] = mr/sin(mtheta);
	HelixPar[4] = -cos(mtheta)/sin(mtheta);

}
	Bool_t
HoughTransformCircleXZ(std::vector<TVector3> gHitPos,
		Int_t *MaxBin, Double_t *HelixPar,
		Int_t MinNumOfHits /*=8*/)
{

	// Equation
	// (x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
	// p = r * Const * dMagneticField;
	Double_t dMagneticField = HS_field_0*(valueHSHall/valueHSCalc);
	Bool_t status = true;

	delete gr;
	gr = new TGraph();
	gr->SetMarkerStyle(22);
	gr->SetMarkerSize(2);
	//for TPC circle track
	//Hough-transform
	histCircle -> Reset();
	delete histCircle;
	hc2 -> Reset();
	delete hc2;
	histCircle = new TH3D("histCircle",";rd (mm); theta (rad); p(MeV/c)",
		nBin_rdiff, rdiff_min,  rdiff_max,
		nBin_theta, theta_min, theta_max,
		nBin_p, pmin, pmax);
	hc2 = new TH3D("hc2","f;rd (mm); theta (rad); p(MeV/c)",
		nBin_rdiff, rdiff_min,  rdiff_max,
		nBin_theta, theta_min, theta_max,
		nBin_p, pmin, pmax);
	for(Int_t i=0; i<gHitPos.size(); ++i){
		Double_t x = -gHitPos[i].X();
		Double_t y = gHitPos[i].Z() - ZTarget;
		gr->AddPoint(y,x);
		for(Int_t ird=0; ird<nBin_rdiff; ++ird){
			Double_t rd = histCircle->GetXaxis()->GetBinCenter(ird+1);
			for(Int_t ip=0; ip<nBin_p; ++ip){
				Double_t p = histCircle->GetZaxis()->GetBinCenter(ip+1);
				Double_t r = p/(Const*dMagneticField);

				//a*sin(theta) + b*cos(theta) + c = 0
				Double_t a = 2.*(r+rd)*y;
				Double_t b = 2.*(r+rd)*x;
				Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

				Double_t r0 = sqrt(a*a + b*b);
				if(fabs(-1.*c/r0)>1.) continue;
				Double_t theta1_alpha = asin(-1.*c/r0);

				Double_t theta2_alpha;
				if(theta1_alpha>0.) theta2_alpha = acos(-1.) - theta1_alpha;
				else theta2_alpha = -1.*acos(-1.) - theta1_alpha;

				Double_t theta_alpha = TMath::ATan2(b, a);

				Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
				Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
				Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

				Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
				Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
				Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

				Double_t theta1 = TMath::ATan2(ycenter1, xcenter1);
				Double_t theta2 = TMath::ATan2(ycenter2, xcenter2);
				if(TMath::IsNaN(theta1)){
					std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
						<<"rd="<<rd<<", r"<<r<<std::endl;
				}
				if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
					std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
					std::cout<<"x:"<<x<<", y:"<<y
						<<", theta1:"<<theta1<<", theta2:"<<theta2
						<<", theta_alpha:"<<theta_alpha<<std::endl;
				}
				histCircle->Fill(rd, theta1, p);
				histCircle->Fill(rd, theta2, p);
				hc2->Fill(rd,theta1,p);
				hc2->Fill(rd,theta2,p);
			} // p bins
		} // rdiff bins
	} //i


	Int_t maxbin = histCircle->GetMaximumBin();
	Int_t mx, my, mz;
	histCircle->GetBinXYZ(maxbin, mx, my, mz);
	MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

	Double_t hough_rd = histCircle->GetXaxis()->GetBinCenter(mx);
	Double_t hough_theta = histCircle->GetYaxis()->GetBinCenter(my);
	Double_t hough_p = histCircle->GetZaxis()->GetBinCenter(mz);
	HelixPar[3] = hough_p/(Const*dMagneticField); //helix r
	HelixPar[0] = (HelixPar[3] + hough_rd)*cos(hough_theta); //helix cx
	HelixPar[1] = (HelixPar[3] + hough_rd)*sin(hough_theta); //helix cy
	
//	histCircle->GetXaxis()->SetRange(mx,mx);
	hc2->GetZaxis()->SetRange(800,800);
	hist_temp = (TH2D*)hc2->Project3D("xy");
	hc2->GetZaxis()->SetRange(0,nBin_p);
	hc2->GetYaxis()->SetRange(my,my);
	hist_temp2 = (TH2D*)hc2->Project3D("xz");
	hc2->GetYaxis()->SetRange(0,nBin_theta);
	histCircle->GetZaxis()->SetRange(450,450);
	hist_temp3 = (TH2D*)histCircle->Project3D("xy");

	cout<<Form("Theta, p = (%f, %f)",hough_theta, hough_p)<<endl;

	if(histCircle->GetMaximum() < 0.5*MinNumOfHits) status = false;
	return status;
}



#endif
