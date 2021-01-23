#ifndef _LSLRmodel_H
#define _LSLRmodel_H

class LSLRmodel {
private:
Rcpp::NumericVector data;
Rcpp::NumericVector positions;
int regression_order;
int dist_num;
int npar;
double limit;

public:
LSLRmodel(SEXP);
LSLRmodel( 	
	Rcpp::NumericVector _data,
	Rcpp::NumericVector _positions,
	int _regression_order,
	int _dist_num,
	int _npar,
	double _limit
);

std::vector<double> LSLRfit();
double dR2dx(double X);

};

#endif 