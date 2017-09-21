#include "WeibullR.h"	

	using namespace Rcpp;
	
SEXP CallgetCCC2(SEXP arg1, SEXP arg2)  {

int F=Rcpp::as<int>(arg1);
int model=Rcpp::as<int>(arg2);

return wrap(getCCC2(F, model));
}
