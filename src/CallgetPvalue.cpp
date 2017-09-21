#include "WeibullR.h"	

	using namespace Rcpp;
	
SEXP CallgetPvalue(SEXP arg1, SEXP arg2, SEXP arg3)  {

int F=Rcpp::as<int>(arg1);
double R2=Rcpp::as<double>(arg2);
int model=Rcpp::as<int>(arg3);


struct AbPval Zzstruct;
Zzstruct=getPvalue(F,R2,model);

NumericVector Outvals(2);

Outvals[0]=Zzstruct.Pval;
Outvals[1]=Zzstruct.CCC2;

return Outvals;
}
