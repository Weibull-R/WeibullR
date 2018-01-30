#ifndef _MLEmodel_H
#define _MLEmodel_H


    using namespace Rcpp ;

class MLEmodel {

arma::colvec time;
arma::colvec qty;
Rcpp::NumericVector N;

//double failcomp;
//double suscomp;
//double discomp;
//double intcomp;
int endf;
int ends;
int endd;
int endil;
int endir;
arma::colvec fail;
arma::colvec nf;
arma::colvec susp;
arma::colvec ns;
arma::colvec disc;
arma::colvec nd;
arma::colvec left;
arma::colvec right;
arma::colvec ni;

public:
MLEmodel();
MLEmodel(SEXP);
double LogLike(arma::colvec, int, int, double);
double tryLL(arma::colvec, int);
SEXP MLE_Simplex(SEXP, arma::colvec, double, int);
SEXP dMaxLLdx( SEXP, arma::colvec, double);
};
// end of class declaration


#endif