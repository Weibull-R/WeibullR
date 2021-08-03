#ifndef _MLEmodel_H
#define _MLEmodel_H


    using namespace Rcpp ;

class MLEmodel {

//arma::colvec time;
//arma::colvec qty;
// the N vector is used in Loglike and tryLL
Rcpp::NumericVector N;


//int endf;
//int ends;
//int endd;
//int endil;
//int endir;

// nf, ns, dn and ni are the qty vectors
arma::colvec fail;
arma::colvec nf;
arma::colvec susp;
arma::colvec ns;
arma::colvec disc;
arma::colvec nd;
arma::colvec left;
arma::colvec right;
arma::colvec ni;

// dist_num is required for any MLEmodel
int dist_num;

// members holding simplex control
double limit;
int maxit;

public:
MLEmodel();
MLEmodel(SEXP);
double LogLike(arma::colvec, int, double);
//double tryLL(arma::colvec, int);
double tryLL(arma::colvec);
double tryLL2(arma::colvec, double);
//SEXP MLE_Simplex(SEXP, arma::colvec, double, int);
SEXP MLE_Simplex(arma::colvec, double, int);
SEXP MLE3p(SEXP arg3, SEXP arg4, SEXP arg5); 
//new support for 3p by seek
std::vector<double> genTzvec(double start, double end, double maxtz, int npts) ;
//old support for 3p by secant
//SEXP dMaxLLdx( SEXP, arma::colvec, double);

// these permit more descriptive code outside of primary constructor
void setSimplexLimit(double);
void setSimplexMaxit(int);

// important access
int getDistNum();

};
// end of class declaration


#endif
