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
// exported for testing purposes only. Never called from R.
RcppExport SEXP MLEdMaxLLdx(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4);
//******************//
// end abremDebias


// abremPivotals code
//RcppExport SEXP LSLRw2pXonY (SEXP arg1, SEXP arg2);
//RcppExport SEXP LSLRw2pYonX (SEXP arg1, SEXP arg2);
//RcppExport SEXP LSLRln2pXonY (SEXP arg1, SEXP arg2);
//RcppExport SEXP LSLRln2pYonX (SEXP arg1, SEXP arg2);
//RcppExport SEXP LSLRg2pXonY (SEXP arg1, SEXP arg2);
//RcppExport SEXP LSLRg2pYonX (SEXP arg1, SEXP arg2);

RcppExport SEXP LSLR(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4);
RcppExport SEXP pivotalMC(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6,SEXP arg7,SEXP arg8);	
RcppExport SEXP CallgetCCC2(SEXP arg1, SEXP model);
RcppExport SEXP CallgetPvalue(SEXP arg1, SEXP arg2, SEXP arg3);

extern "C" double getCCC2( int F, int model);
extern "C" struct AbPval getPvalue(int F, double R2, int model);
// end abremPivotals

#endif
#endif
