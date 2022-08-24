#include "Math.hh"
#ifndef Minimizer_h
#define Minimizer_h
class Minimizer{
	private:
		static const int max_param=10;
		double m_value=0;
		double m_pars[max_param]={0};
		double m_parMin[max_param]={0};
		double m_parMax[max_param]={0};
		int m_npar=0;
		double cut = 0.1;
		int max_loop = 1000;
		function<double(double*)>m_f;
	public:
		Minimizer(int npar,function<double(double*)>&f,double* pars){
			m_npar=npar;
			m_f=f;
			for(int i=0;i<m_npar;++i){
				m_pars[i]=pars[i];
			}
		}
		bool CheckParamLimits(int par);
		double dF(int par,double step);
		void Minimize(double* pars,double* steps);
		
};




#endif
