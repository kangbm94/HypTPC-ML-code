#include "../include/Minimizer.hh"
double Minimizer::dF(int par, double step){
	double par_asc[max_param]={0};
	double par_dsc[max_param]={0};
	for(int i=0;i<m_npar;++i){
		par_asc[i]=m_pars[i];
		par_dsc[i]=m_pars[i];
	}
	par_asc[par]+=step/2;
	par_dsc[par]-=step/2;
	return (m_f(par_asc)-m_f(par_dsc))/step;
}
bool Minimizer::CheckParamLimits(int par){
	if (m_pars[par]<m_parMax[par]&&m_parMin[par]<m_pars[par]) return true;
	else return false;
}
void Minimizer::Minimize(double* pars,double* steps){
	double f0 = m_f(m_pars);	
	double df[max_param]={0};
	bool ReCalc=true;
	int cnt=0;
	while(ReCalc&&cnt<max_loop){
		double f_before = m_f(m_pars);
		for(int i=0;i<m_npar;++i){
			df[i]=dF(i,steps[i]);
		}
		for(int i=0;i<m_npar;++i){
			m_pars[i]+=df[i]*steps[i];
		}
		double f_after = m_f(m_pars);
		if(f_before-f_after<cut){
			ReCalc=false;
		}
		cnt++;
	}
	for(int i=0;i<m_npar;++i){
		pars[i]=m_pars[i];
	}
}
