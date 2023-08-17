#ifndef Hough_hh
#define Hough_hh
Bool_t HoughTransformLineXZ(std::vector<TVector3> gHitPos,
			    Int_t *MaxBin, Double_t *LinearPar,
			    Int_t MinNumOfHits=8);
//HS-ON
Bool_t HoughTransformCircleXZ(std::vector<TVector3> gHitPos, Int_t *MaxBin, Double_t *HelixPar,
			      Int_t MinNumOfHits=8);


void HoughTransformLineYZ(std::vector<TVector3> gHitPos,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);
void HoughTransformLineYX(std::vector<TVector3> gHitPos,
			  Int_t *MaxBin, Double_t *LinearPar,
			  Double_t MaxHoughWindowY);
//HS-ON
void HoughTransformLineYTheta(std::vector<TVector3> gHitPos,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Double_t MaxHoughWindowY);

#endif
