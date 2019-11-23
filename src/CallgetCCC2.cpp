#include "WeibullR.h"
#include <math.h>

	using namespace Rcpp;

SEXP CallgetCCC2(SEXP arg1, SEXP arg2)  {

int F=Rcpp::as<int>(arg1);
int model=Rcpp::as<int>(arg2);

// The following code was originally preprepared as a separate ANSI C compliant function, getCCC2(int,int).
// The thinking was to facilitate wider use of the function in other software.
// It has been modified here as C++ due to apparant bugs in the MinGW compiler when adding such C code
// errors like "cc1plus.exe: out of memory allocating 139178 bytes" would occur randomly
// and only on the R-Forge build system for 64-bit Windows.

double T1w[23]={0.792235,0.7990604,0.8076126,0.8204102,0.832331,0.8425375,0.8514909,0.8593213,0.8662665,
 0.8724075,0.8779149,0.8828711,0.887337,0.8914107,0.8951199,0.8985272,0.9016913,0.9045873,
 0.9073225,0.9098415,0.9121919,0.9144172,0.9164957};
double T2w[26]={2.482857,2.593721,2.689116,2.773029,2.847969,2.915728,2.977454,3.034384,3.087205,3.136568,
 3.182761,3.226223,3.267092,3.306225,3.343156,3.378561,3.412309,3.444714,3.475747,3.505642,
 3.534433,3.562236,3.588777,3.614618,3.639701,3.663851};
double T3w[28]={3.663851,3.872307,4.037513,4.17472,4.292113,4.394901,4.486445,4.568626,4.643578,4.712588,
 4.776053,4.83515,4.890495,4.942337,4.991377,5.037647,5.081707,5.12344,5.163087,5.2011,
 5.237467,5.272338,5.305946,5.338348,5.369513,5.399487,5.428722,5.456727};
double T1l[23]={0.7938923,0.7992166,0.8143357,0.8286594,0.8416131,0.8531055,0.863076,0.8717764,0.8794219,
    0.8862083,0.8921895,0.8975986,0.9024265,0.9068011,0.9107908,0.9144347,0.9177708,0.9208458,
    0.9236726,0.9262948,0.9287454,0.931017,0.9331573};
double T2l[26]={2.705413,2.847212,2.969813,3.077389,3.173427,3.260117,3.339296,3.412094,3.479407,3.542081,
    3.600777,3.655789,3.707801,3.756996,3.803559,3.847988,3.890183,3.93063,3.969467,4.006691,
    4.042462,4.076862,4.109992,4.142034,4.173005,4.202877};
double T3l[28]={4.202877,4.458735,4.659025,4.823861,4.963904,5.085855,5.193705,5.290567,5.378335,5.45858,
    5.532547,5.601127,5.6651,5.725016,5.781398,5.834598,5.884945,5.932743,5.978272,6.0218,
    6.06339,6.103111,6.141497,6.178083,6.213509,6.2477,6.280481,6.312331};

double T1[23]={0.0},T2[26]={0.0},T3[28]={0};
double CCC2=0.0;
int i=0;
int Fbl=0;
int Fbu=0;
double qwl=0.0;
double qwu=0.0;
double qwccc2=0.0;

if(model==1)  {
    for( i=0; i<23; T1[i]=T1w[i], i++);
    for( i=0; i<26; T2[i]=T2w[i], i++);
    for( i=0; i<28; T3[i]=T3w[i], i++);
}
else   {
    for( i=0; i<23; T1[i]=T1l[i], i++);
    for( i=0; i<26; T2[i]=T2l[i], i++);
    for( i=0; i<28; T3[i]=T3l[i], i++);
}
 if(F<26)  {
  CCC2=T1[F-3];
  }else{
   if(F<151)  {
    i=5;
    if(F%i==0)  {
/* The qweibull of CCC2 can be taken directly from T2
 * Offset value is 25/i-1 (for R) Will be 25/i for C++   */
     CCC2=1-1/exp(T2[F/i-5]);
    }else{
/* The qweibull of CCC2 will have to be interpolated from T2
 * establish  F and qweibull bounds   */
     Fbl=i*(F/i);
     Fbu=Fbl+i;
     qwl=T2[Fbl/i-5];
     qwu=T2[Fbu/i-5];
/* Then interpolate using log(F) and log(Fbounds) */
/* casting integers to double for arguments to the log function */
/* for complieance using platform: i386-pc-solaris2.10 (32-bit) */
     qwccc2=qwl+((log((double) F)-log((double) Fbl))/(log((double) Fbu)-log((double) Fbl))*(qwu-qwl));
     CCC2=1-1/exp(qwccc2);
    }
   }else{
    if(F<1501)   {
     i=50;
     if(F%i==0)  {
/* The qweibull of CCC2 can be taken directly from T3
 * Note there is a difference in the F/i offset for element selection!!!
 * In this case the offset = 150/i-1 (for R),will be 150/i for C++  */
      CCC2=1-1/exp(T3[F/i-3]);
     }else{
/* The qweibull of CCC2 will have to be interpolated from T3 */

/* establish  F and qweibull bounds  */
      Fbl=i*(F/i);
      Fbu=Fbl+i;
      qwl=T3[Fbl/i-3];
      qwu=T3[Fbu/i-3];
/* Then interpolate using log(F) and log(Fbounds) */
      qwccc2=qwl+((log((double) F)-log((double) Fbl))/(log((double) Fbu)-log((double) Fbl))*(qwu-qwl));
      CCC2=1-1/exp(qwccc2);
     }
    }
   }
  }



//return wrap(getCCC2(F, model));
return wrap(CCC2);
}
