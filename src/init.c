/*
library(tools)
package_native_routine_registration_skeleton("C:/Users/Dad_laptop/Documents/Github/WeibullR")
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP adjustedRank(SEXP);
//extern SEXP CallgetCCC2(SEXP, SEXP);
//extern SEXP CallgetPvalue(SEXP, SEXP, SEXP);
extern SEXP getContour(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LSLR(SEXP, SEXP, SEXP, SEXP);
extern SEXP MLEdMaxLLdx(SEXP, SEXP, SEXP, SEXP);
extern SEXP MLEloglike(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MLEsimplex(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pivotalMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP plotData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"adjustedRank",  (DL_FUNC) &adjustedRank,  1},
//    {"CallgetCCC2",   (DL_FUNC) &CallgetCCC2,   2},
//    {"CallgetPvalue", (DL_FUNC) &CallgetPvalue, 3},
    {"getContour",    (DL_FUNC) &getContour,    7},
    {"LSLR",          (DL_FUNC) &LSLR,          4},
    {"MLEdMaxLLdx",   (DL_FUNC) &MLEdMaxLLdx,   4},
    {"MLEloglike",    (DL_FUNC) &MLEloglike,    5},
    {"MLEsimplex",    (DL_FUNC) &MLEsimplex,    5},
    {"pivotalMC",     (DL_FUNC) &pivotalMC,     8},
    {"plotData",      (DL_FUNC) &plotData,      7},
    {NULL, NULL, 0}
};

void R_init_WeibullR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

