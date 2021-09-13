/*
library(tools)
 package_native_routine_registration_skeleton("C:/Rpack/WeibullR")
*/


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP adjustedRank(SEXP);
extern SEXP callMLE3p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getContour(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LSLR(SEXP);
extern SEXP MLEloglike(SEXP, SEXP, SEXP, SEXP);
extern SEXP MLEsimplex(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pivotalMC(SEXP);
extern SEXP plotData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"adjustedRank", (DL_FUNC) &adjustedRank, 1},
    {"callMLE3p",    (DL_FUNC) &callMLE3p,    6},
    {"getContour",   (DL_FUNC) &getContour,   6},
    {"LSLR",         (DL_FUNC) &LSLR,         1},
    {"MLEloglike",   (DL_FUNC) &MLEloglike,   4},
    {"MLEsimplex",   (DL_FUNC) &MLEsimplex,   5},
    {"pivotalMC",    (DL_FUNC) &pivotalMC,    1},
    {"plotData",     (DL_FUNC) &plotData,     8},
    {NULL, NULL, 0}
};

void R_init_WeibullR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
