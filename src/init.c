// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
//File 'WeibullR/libs/x64/WeibullR.dll':
//  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
void R_init_WeibullR(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}