#ifndef _Pivotal_H
#define _Pivotal_H

class Pivotal {
private:
Rcpp::NumericVector positions;
Rcpp::NumericVector event;
double R2test;
double CItest;
double P1;
double P2;
unsigned int S;
int seed;
Rcpp::NumericVector dp;
int regression_order;
int dist_num;
int npar;
double limit;
bool ProgRpt;

public:
Pivotal(SEXP);
SEXP Execute();

};

#endif
