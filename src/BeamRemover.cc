#include "../include/BeamRemover.hh"
BeamRemover::BeamRemover(double xmin,double xmax,double ymin ,double ymax ){
	Ci_hist = new TH3D("hist_circle","rd(mm);theta(rad);p(MeV/c)",
			nBin_rdiff,rdiff_min,rdiff_max,
			nBin_theta,theta_min,theta_max,
			nBin_p,p_min,p_max);
	histY = new TH2D("histY","theta(deg),r(mm)",
			thetaY_ndiv,thetaY_min,thetaY_max,
			r_ndiv,r_min,r_max);
	Ydistrib = new TH1D("Ydistribution","Ydistribution",140,-350,350);
	x_min=xmin,x_max=xmax;
	y_min=ymin,y_max=ymax;
	ZXHist = new TH2D("CircleHist","CircleHist",100,-250,250,100,-250,250);
	Arc = new TF1("Arc","-TMath::Sqrt([2]*[2]-(x-[1])*(x-[1]))+[0]",-150,400);
}
void BeamRemover::DoCircleHough(int i){
	int nh = hitarray[i].size();
	bool prev_add = true;
	MaxNBeam = 3;
	hough_rd.resize(MaxNBeam);
	hough_theta.resize(MaxNBeam);
	hough_p.resize(MaxNBeam);
	hough_flag[i].resize(nh,0);
	//	TCanvas* c1 = new TCanvas("chough","chough",800,800);
	for(int ib = 0 ;ib<MaxNBeam;++ib){
		Ci_hist->Reset();
		for(int ih = 0; ih<nh;++ih){
			if(hough_flag[i].at(ih)>0) continue;
			auto pos = hitarray[i].at(ih);
			ZXHist->Fill(pos.Z(),pos.x());
			for(int ird = 0;ird<nBin_rdiff;++ird){
				double rd = Ci_hist->GetXaxis()->GetBinCenter(ird+1);
				for(int ip=0;ip<nBin_p;++ip){
					double x = -pos.x();
					double y = pos.z()-ZTarget;
					double p = Ci_hist->GetZaxis()->GetBinCenter(ip+1);
					double r = p/(Const*dMagneticField);
					double a = 2.*(r+rd)*y;	
					double b = 2.*(r+rd)*x;	
					double c = -1.*(rd*rd+2.*r*rd+x*x+y*y);
					double r0 = sqrt(a*a+b*b);
					if(fabs(-1*c/r0)>1){
						continue;
					}
					double theta1_alpha = asin(-1.*c/r0);
					double theta2_alpha;
					if(theta1_alpha>0.)
						theta2_alpha = acos(-1.) - theta1_alpha;
					else
						theta2_alpha = -1.*acos(-1.) - theta1_alpha;
					double theta_alpha = atan2(b, a);
					double xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
					double ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
					double r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

					double xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
					double ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
					double r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));
					double theta1 = atan2(ycenter1, xcenter1)+2*acos(-1);
					double theta2 = atan2(ycenter2, xcenter2)+2*acos(-1);
					Ci_hist->Fill(rd, theta1, p);
					Ci_hist->Fill(rd, theta2, p);
				}// ip
			}//ird
		}//ih
		int maxbin = Ci_hist->GetMaximumBin();
		int mx,my,mz;
		Ci_hist->GetBinXYZ(maxbin, mx, my, mz);
		/*
			 for(int i=0;i<hough_x.size();++i){
			 int bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i])+fabs(mz-hough_z[i]);
			 if(bindiff<=2) hough_flag = false;
			 }
			 */
		hough_x.push_back(mx);
		hough_y.push_back(my);
		hough_z.push_back(mz);
		//    if(!hough_flag) continue;
		hough_rd[ib] = Ci_hist->GetXaxis()->GetBinCenter(mx);
		hough_theta[ib] = Ci_hist->GetYaxis()->GetBinCenter(my);
		hough_p[ib] = Ci_hist->GetZaxis()->GetBinCenter(mz);
		double hough_r = hough_p[ib]/(Const*dMagneticField);
		double hough_cx = (hough_r + hough_rd[ib])*cos(hough_theta[ib]);
		double hough_cy = (hough_r + hough_rd[ib])*sin(hough_theta[ib]);
		int hough_count = 0;
		double minz = 1115.683;
		double maxz = -1115.683;
		vector<double> zdist;
		for( int ih = 0; ih < nh; ++ih){	
			//			if(hough_flag[i].at(ih)>0) continue;
			auto pos = hitarray[i].at(ih);
			double x = -pos.X();
			double y = pos.Z() - ZTarget;
			double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
			double dist = abs(r_cal-hough_r);
			if(dist < MaxHoughWindow){
				if(pos.Z()<minz) minz = pos.Z();
				if(pos.Z()>maxz) maxz = pos.Z();
				hough_flag[i].at(ih) += pow(2,ib);
				hough_count++;
				zdist.push_back(pos.Z());
			}
		}
		sort(zdist.begin(),zdist.end());
		int ngap = 0;
		for(int iz = 0; iz<hough_count-1;++iz){
			double dist = abs(zdist.at(iz)-zdist.at(iz+1));
			if(dist>50){
				cout<<"Zdist : "<<dist<<endl;
	//			ngap++;
				if(dist>125) hough_count = 0;
				//		break;
			}
		}
		if(ngap>3) hough_count = 0;

		if(hough_count > MinBeamHit and minz < MinZCut and maxz > MaxZCut){
			cout<<"Minz = "<<minz<<endl;
			h_cx[i].push_back(hough_cx);
			h_cy[i].push_back(hough_cy);
			h_r[i].push_back(hough_r);
			for( int ih = 0; ih < nh; ++ih){
				auto pos = hitarray[i].at(ih);
				bool beamflag = hough_flag[i].at(ih)%int(pow(2,ib+1))/int(pow(2,ib));
				if(beamflag){
					for(int ihist = 0;ihist<100*(ib+1);ihist++)ZXHist->Fill(pos.Z(),pos.X());
				}
			}
		}
		//		histY->Reset();
	}//ib
	ZXHist->Reset("ice");
	for( int ih = 0; ih < nh; ++ih){
		int beamid = -1;
		auto pos = hitarray[i].at(ih);
		ZXHist->Fill(pos.Z(),pos.X());
		if(hough_flag[i].at(ih)>0){
			beamid = CompareHough(pos,h_cx[i],h_cy[i],h_r[i]);
		}
		if(beamid > -1){
			hough_flag[i].at(ih) = pow(2,beamid); 
			for(int ib=0;ib<100*(beamid+1);++ib)ZXHist->Fill(pos.Z(),pos.X());
		}
	}
	int nbeam = h_cx[i].size();
	for(int ib = 0;ib<nbeam;++ib){
		double hcx = h_cx[i].at(ib);	
		double hcy = h_cy[i].at(ib);	
		double hr = h_r[i].at(ib);	
		ArcGraph = new TGraph(); 	
		for(int ih = 0; ih<nh; ++ih){
			if(hough_flag[i].at(ih) == pow(2,ib)){
				auto pos = hitarray[i].at(ih);
//				ArcGraph->AddPoint(pos.Z()-ZTarget,-pos.X());
				ArcGraph->SetPoint(ArcGraph->GetN(),pos.Z()-ZTarget,-pos.X());
			}
		}
		cout<<Form("Params : (%f,%f,%f)",hcx,hcy,hr)<<endl;
		Arc->SetParameter(0,hcx);
		Arc->SetParameter(1,hcy);
		Arc->SetParameter(2,hr);
		Arc->SetParLimits(2,0.5*hr,1.5*hr);
		ArcGraph->Fit("Arc","QR0");
		h_cx[i].at(ib)=Arc->GetParameter(0);
		h_cy[i].at(ib)=Arc->GetParameter(1);
		h_r[i].at(ib)=Arc->GetParameter(2);
		delete ArcGraph;
	}
}

void
BeamRemover::DoYThetaHough(int i ){
	int nb = h_cx[i].size();
	int nh = hitarray[i].size();
	for(int ib = 0; ib<nb;++ib){
		double hcx = h_cx[i].at(ib);
		double hcy = h_cy[i].at(ib);
		double hr = h_r[i].at(ib);
		histY->Reset();
		for(int ih = 0;ih<nh;++ih){
			if(hough_flag[i].at(ih) != pow(2,ib)) continue;
			auto pos = hitarray[i].at(ih);
			double x = -pos.X();
			double y = pos.Z()-ZTarget;
	  		for(int ti=0; ti<thetaY_ndiv; ti++){
	    		double theta = thetaY_min+ti*(thetaY_max-thetaY_min)/thetaY_ndiv;
					double tmpx = -pos.X();
					double tmpy = pos.Z()-ZTarget;
					double tmpz = pos.Y();
	    		double tmp_t = atan2(tmpy - hcy,tmpx-hcx);
	    		double tmp_xval = hr * tmp_t;
	    	
					histY->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval+sin(theta*acos(-1.)/180.)*tmpz);
				}
		}//ih
    int maxbinY = histY->GetMaximumBin();
    int mxY,myY,mzY;
    histY->GetBinXYZ(maxbinY, mxY, myY, mzY);
    double mtheta = histY->GetXaxis()->GetBinCenter(mxY)*acos(-1.)/180.;
    double mr = histY->GetYaxis()->GetBinCenter(myY);
    cout<<"rad = "<<mr<<", theta "<<mtheta<<endl;
		double p0 = mr/sin(mtheta);
    double p1 = -cos(mtheta)/sin(mtheta);
		h_z0[i].push_back(p0);
		h_dz[i].push_back(p1);
		for(int ih = 0;ih<nh;++ih){
			if(hough_flag[i].at(ih) != pow(2,ib)) continue;
			auto pos = hitarray[i].at(ih);
	  	double tmpx = -pos.x();
	  	double tmpy = pos.z() - ZTarget;
	  	double tmpz = pos.y();
	  	double tmp_t = atan2(tmpy - hcy,tmpx - hcx);
	  	double tmp_xval = hr * tmp_t;
	  	double distY = fabs(p1*tmp_xval - tmpz + p0)/sqrt(pow(p1, 2) + 1);
			if(distY > MaxHoughWindowY) hough_flag[i].at(i) = 0;
		}
	}
}



vector<TEllipse*>
BeamRemover::GetCircle(int i){
	circles.clear();
	int nb = h_cx[i].size();
	for(int ib=0;ib<nb;++ib){
		double x = -h_cx[i].at(ib);
		double y = h_cy[i].at(ib)+ZTarget;
		double r = h_r[i].at(ib);
		auto cir = new TEllipse(y,x,r);
		double mom = r* Const*dMagneticField;
		cout<<"Mom = " << mom<<endl ;
		cir->SetFillStyle(0);
		circles.push_back(cir);
	}
	return circles;
}

int
BeamRemover::SearchPeaks(TH1D* hist,std::vector<double> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.2;
	int npeaks = spec.Search(hist,sig,"",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}
void 
BeamRemover::AssignHits(std::vector<double>peaks){
	int np = peaks.size();
	hitarray.resize(np);
	h_cx.resize(np);
	h_cy.resize(np);
	h_z0.resize(np);
	h_r.resize(np);
	h_dz.resize(np);
	hough_flag.resize(np);
	for(auto hitp:rawhitarray){
		double y = hitp.Y();
		for(int ipeak = 0; ipeak<np;ipeak++){
			double peak = peaks.at(ipeak);
			if(abs(y-peak)<Ywidth){
				hitarray[ipeak].push_back(hitp);
				break;
			}
		}
	}
}
int
BeamRemover::CompareHough(TVector3 pos, vector<double> hcx,vector<double>hcy,vector<double> hr){ 
			
	int beamid = -1;
	double x = -pos.X();
	double y = pos.Z() - ZTarget;
	int nbeam = hcx.size();
	double min_dist = MaxHoughWindow;
	for(int ib=0;ib<nbeam;++ib){
		double hough_cx = hcx.at(ib),hough_cy = hcy.at(ib),hough_r = hr.at(ib);
		double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
		double dist = abs(r_cal-hough_r);
		if(dist<min_dist){
			min_dist = dist;
			beamid = ib;
		}
	}
	return beamid;
}








void 
BeamRemover::LoadHit(TVector3 hit){
	double x = hit.X(),y=hit.Y();
	if(x_min<x and x<x_max and y_min<y and y<y_max){}
	else{
		Ydistrib->Fill(y);
		rawhitarray.push_back(hit);
	} 
}
BeamRemover::~BeamRemover(){};

















