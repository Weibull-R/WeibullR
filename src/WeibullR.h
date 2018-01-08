#ifndef _WeibullR_H
#define _WeibullR_H

struct AbPval{
	double Pval,CCC2;
};


#ifdef __cplusplus

#include <RcppArmadillo.h>

// abremDebias code
RcppExport SEXP MLEloglike(SEXP arg1, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6);
RcppExport SEXP MLEsimplex(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5);
// used in secant method for determination of 3rd parameter optimization
RcppExport SEXP MLEdMaxLLdx(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4);
//******************//
// end abremDebias


// abremPivotals code
RcppExport SEXP LSLR(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4);
RcppExport SEXP pivotalMC(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6,SEXP arg7,SEXP arg8);
RcppExport SEXP CallgetCCC2(SEXP arg1, SEXP model);
RcppExport SEXP CallgetPvalue(SEXP arg1, SEXP arg2, SEXP arg3);

extern "C" double getCCC2( int F, int model);
extern "C" struct AbPval getPvalue(int F, double R2, int model);
// end abremPivotals

// contour code
RcppExport SEXP getContour(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7);

extern "C" void R_init_WeibullR(DllInfo* info);

#endif
#endif
