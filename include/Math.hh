#include "Constants.hh"
#ifndef Math_h
#define Math_h
double nand = std::numeric_limits<double>::quiet_NaN();
double square(double a){ 
	return a*a;
}
double Norm(double x1,double x2=0,double x3=0,double x4=0){
	return sqrt(square(x1)+square(x2)+square(x3)+square(x4));
}
int PM(int a){
	if(a%2==0){
		return 1;
	}
	else{
		return -1;
	}
}
double Min(double a, double b){
	if(a>b){
		return b;
	}
	else{
		return a;
	}
}
double Max(double a, double b){
	if(a>b){
		return a;
	}
	else{
		return b;
	}
}
double Power(double x,int n){
	double value=1;
	for(int i=0;i<n;i++){
		value=value*x;
	}
	return value;
}
double Polynomial(double x, double* p,int n){
	double value=0;
	for(int i=0;i<n+1;i++){
		value+=p[i]*Power(x,i);
	}
	return value;
}
double SquareSum(double a, double b){
	return sqrt(a*a+b*b);
}
double NormGaussian(double x, double mean, double sigma, double amplitude){
	double par=(x-mean)/sigma;
	double val = amplitude*exp(-par*par/2)/sigma/sqrt(2*Pi());
	return val;
}
double Gaussian(double x, double mean, double sigma, double peak){
	double par=(x-mean)/sigma;
	double val = peak*exp(-par*par/2);
	return val;
}
double fGaussian(double* x,double* p){
	return Gaussian(x[0],p[0],p[1],p[2]);
}
double Step(double a){
	if(a>0){
		return 1;
	}
	else{
		return 0;
	}
}
//
double QuadRoot(double a,double b,double c, int conf){
	if(conf<0){
		return  (-b -sqrt(b*b-4*a*c))/(2*a);
	}
	else{
		return (-b +sqrt(b*b-4*a*c))/(2*a);
	}
}
double T1(double a,double b, double c, double d){
	double p0=1/(3*a);
	double p1=2*b*b*b-9*a*b*c+27*a*a*d;
	double p2 = b*b-3*a*c;
	double val= -b*p0-p0*(
			pow((p1+sqrt(p1*p1-4*p2*p2*p2))/2,1./3)
			+pow((p1-sqrt(p1*p1-4*p2*p2*p2))/2,1./3)
			);
	return val;
}



//Randoms
double Rndm(double r1=0,double r2=1){
	double tmp = r2;
	if(r1>r2){
		r2=r1;
		r1=tmp;
	}
	double range = r2-r1;
	return r1+range*gRandom->Rndm();
}
double GenUniformRandom(double range = 1.){
	double rnd=gRandom->Rndm();
	return range*rnd;
}
void GenSphericalRandom(double &theta, double &phi){
	phi=GenUniformRandom(2*Pi());
	double cth= GenUniformRandom(2)-1;
	theta=acos(cth);
}
TVector2 GenCircleRandom(double r, double phi1,double phi2){
	double pi = Pi();
	double phi=pi*Rndm(phi1/pi, phi2/pi);
	double x = r*cos(phi),y=r*sin(phi);
	return TVector2(x,y);
}
TVector2 Generate2DGaus(double z0,double x0,double sig){
	//double r = abs(gRandom->Gaus(0,sig));
	double dz = gRandom->Gaus(0,sig);
	double dx = gRandom->Gaus(0,sig);
	TVector2 V(z0,x0);
	TVector2 dV(dz,dx);
	return V+dV;
	//	return V+GenCircleRandom(r,0,2*acos(-1));
}



double chebyshev(double x, int n){
	double val;
	switch(n){
		case 0:
			val=1;
			break;
		case 1: 
			val=x;
			break;
		case 2: 
			val= 2*x*x-1;
			break;
		case 3:
			val= 4*x*x*x-3*x;
			break;
		case 4: 
			val= 8*x*x*x*x-8*x*x+1;
			break;
		case 5: 
			val= 16*x*x*x*x*x-20*x*x*x+5*x;
			break;
		case 6: 
			val= 32*x*x*x*x*x*x-48*x*x*x*x+18*x*x-1;
			break;
	}
	return val;
}
void CircleFit(std::vector<TVector3> posarr,double* param){
	int n =posarr.size();
	double Sumx=0,Sumy=0;
	double Sumx2=0,Sumy2=0,Sumxy=0;
	double Sumx3=0,Sumy3=0,Sumx2y=0,Sumy2x=0;
	for (Int_t i=0;i<n;i++) {
		double x = posarr.at(i).X();
		double y = posarr.at(i).Y();
		Sumx+=x;Sumy+=y;
		Sumx2+=x*x;Sumy2+=y*y;Sumxy+=x*y;
		Sumx3+=x*x*x;Sumy3+=y*y*y;Sumx2y+=x*x*y;Sumy2x+=y*y*x;
	}
	double a_1 = Sumx3 + Sumy2x, a_2 = Sumx2, a_3 = Sumx;
	double b_1 = Sumy3 + Sumx2y, b_2 = Sumy2, b_3 = Sumy;
	double c_1 = a_2+b_2,c_2 = Sumx,c_3 = Sumy;
	double A = (a_1*(b_2*n-b_3*c_3)+a_3*(b_1*c_3-b_2*c_1)+Sumxy*(b_3*c_1-b_1*n))
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double B = (b_1*(a_2*n-a_3*c_2)+b_3*(a_1*c_2-a_2*c_1)+Sumxy*(a_3*c_1-a_1*n))      
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double C = (c_1*(a_2*b_2-Sumxy*Sumxy)+c_2*(Sumxy*b_1-a_1*b_2)+c_3*(Sumxy*a_1-b_1*a_2))
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double cx = -0.5 * A;
	double cy = -0.5 * B;
	double rad = sqrt(cx*cx+cy*cy-C);	
	param[0] = cx;
	param[1] = cy;
	param[2] = rad;
	//	std::cout<<Form("(cx,cy,rad) = (%f,%f,%f)",cx,cy,rad)<<std::endl;
}
void LinearFit(std::vector<TVector3> posarr,double* param){
	// y = ax + b;
	int n =posarr.size();
	double Sumx=0,Sumy=0;
	double Sumx2=0,Sumxy=0;
	for (Int_t i=0;i<n;i++) {
		double x = posarr.at(i).X();
		double y = posarr.at(i).Y();
		Sumx+=x;Sumy+=y;
		Sumx2+=x*x;Sumxy+=x*y;
	}
	double b = (Sumx*Sumxy-Sumx2*Sumy)/( Sumx*Sumx-n*Sumx2);	
	double a = (Sumx*Sumy-n*Sumxy)/( Sumx*Sumx-n*Sumx2);	
	param[0]=b;//b = z0;
	param[1]=a;//a = r*dz;
}






bool SortY(TVector3 a, TVector3 b){
	return a.Y()<b.Y();
}

static bool CircleFitWithRadius(std::vector<TVector3> posarr,double r,double* param, double charge = 0){
	/*
		 Minimize the variance, Sum[((x-a)^2+(y-b)^2-r^2)^2].
		 For numetical stabability, move the frame to CoM.
		 Then, Sum[x] = Sum[y] = 0, some terms cancel out.
		 Sum[
		 a^4  + 2a^2b^2 +b^4

		 -4a^3x -> 0  
		 -4b^3y -> 0

		 -4b^2ax -> 0
		 + 8abxy
		 -4a^2by -> 0

		 +a^2(-2r^2+6x^2+2*y^2) +b^2(-2r^2+2x^2+6y^2)

		 +4a(r^2x-x^3-xy^2)
		 +4b(r^2y-x^2y-y^3)
	//
	constant to a and b. ->neglected
	+r^4
	-2r^2x^2
	-2r^2y^2
	+x^4+2x^2y^2+y^4
	//
	]
	*/
	double gx,gy=0;
	double n = posarr.size();
	//cout<<"CircleFit::n = "<<n<<endl;
	if(n<3) return false;
	for(auto pos:posarr){
		gx+= pos.X()/n;
		gy+= pos.Y()/n;
	}
	double x,y,xx,xy,yy,xxx,xxy,xyy,yyy,xxxx,xxyy,yyyy = 0;
	for(auto pos:posarr){
		double X = pos.X() - gx;
		double Y = pos.Y() - gy;
		xx+=X*X/n;xy+=X*Y/n;yy+=Y*Y/n;
		xxx+=X*X*X/n;xxy+=X*X*Y/n;xyy+=X*Y*Y/n;yyy+=Y*Y*Y/n;
		xxxx+=X*X*X*X/n;xxyy+=X*X*Y*Y/n;yyyy+=Y*Y*Y*Y/n;
	}
//	TString form = Form("pow(x,4)+2*pow(x*y,2)+pow(y,4) +8*x*y*%f+pow(x,2)*%f+pow(y,2)*%f+4*x*%f+4*y*%f" ,
	TString form = Form("x*x*x*x+2*x*x*y*y+y*y*y*y +8*x*y*%f+x*x*%f+y*y*%f+4*x*%f+4*y*%f" ,
			xy,
			-2*r*r+6*xx+2*yy,
			-2*r*r+6*yy+2*xx,
			-xxx-xyy,
			-yyy-xxy
			);
	TF2 fcirc("fcirc",form,-r-200,-r-200,r+200,r+200);
	if(charge == -1){
		fcirc.SetRange(-200,-r-200,r+200,r+200);
	}
	else if (charge == 1){
		fcirc.SetRange(-r-200,-r-200,200,r+200);
	}
//	cout<<Form("Charge = %g",charge)<<endl;
	double cx,cy;
	fcirc.GetMinimumXY(cx,cy);
	cx+=gx;
	cy+=gy;
	param[0]=cx;
	param[1]=cy;
	param[2]=r;
	return true;
}
static bool LinearFitWithSlope(std::vector<TVector3> posarr,double dz,double* param){
		bool crossing=0,positive=0,negative=0;
		double cx = param[0];
		double cy = param[1];
		double rad = param[3];
		for(auto pos : posarr){
			double dx = pos.X()-cx;
			double dy = pos.Y()-cy;
			double t = atan2(dy,dx);
			if( t > acos(-1)*3/4) positive = true;
			if( t < -acos(-1)*3/4) negative = true;
		}
		if(positive and negative) crossing = true;
		double t_sum,z_sum;
		for(auto pos : posarr){
			double dx = pos.X()-cx;
			double dy = pos.Y()-cy;
			double t = atan2(dy,dx);
			if(crossing) t+=2*acos(-1);
			
//			cout<<"( "<<pos.Y()-143<<" ,  "<<-pos.X()<<") angle : "<<t<<endl;
			t_sum+= t;
			z_sum+=pos.Z();	
		}
		double n = posarr.size();
		double z0= - (rad *dz* t_sum - z_sum)/n;
		param[2] = z0;
		param[4] = dz;
		return true;
}


#endif
