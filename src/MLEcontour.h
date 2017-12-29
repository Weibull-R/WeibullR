#ifndef _MLEcontour_H
#define _MLEcontour_H

class MLEcontour {

MLEmodel model;
arma::colvec par_hat;
int dist_num;
double MLLx;
double RatioLL;
double RadLimit;




public:
MLEcontour(SEXP, arma::colvec, int, double, double, double);
arma::rowvec getContourPt( double);

};
// end of class declaration


#endif