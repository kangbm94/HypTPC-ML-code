#include <gmpxx.h>

void test(){
	mpf_set_default_prec(256);
	double a = 3e20;
	mpf_class ma(a);
	double b = 3e-20;
	mpf_class mb(b);
//	ma+mb;
}
