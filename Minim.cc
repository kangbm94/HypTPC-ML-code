#include "./src/Minimizer.cc"

class LC{
	private:
		double a;
	public:
		LC(double a_){
			a=a_;
		}
		string oneStr = "hello";
		function<double(double*)> func= [this](double* b)
  {
			return a+b[1]*b[0];
  };
    
};
void LambdaExp(){
	auto cl = LC(2);
	double b[2]={1,7};
//	auto rt = cl.filePath;
}
void Minim(){
	auto cl = LC(2);
	double b[2]={5,7};
	double g = cl.func(b);
	cout<<g<<endl;
//	auto rt = cl.filePath;
	Minimizer min(2,cl.func,b);
	cout<<min.dF(0,0.1)<<endl;
}
