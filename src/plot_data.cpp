#include "WeibullR.h"

SEXP plotData (SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7)
{
    using namespace Rcpp ;
	NumericVector left  = arg1;
	NumericVector right  = arg2;
	IntegerVector qty = arg3;
	NumericVector tmean  = arg4;
	arma::colvec ptime  = Rcpp::as<arma::colvec >(arg5);
	arma::colvec ppp  = Rcpp::as<arma::colvec >(arg6);
	arma::colvec  adj_rank  = Rcpp::as<arma::colvec >(arg7);

	int n = sum(qty);
	arma::colvec time(n);
	arma::colvec ptpos(n);
	arma::colvec ptrank(n);
	arma::uvec ptwt(n);
	arma::colvec t1(n);
	arma::colvec t2(n);
	arma::colvec lnpos(n);
	arma::colvec lnrank(n);
	arma::uvec lnwt(n);

	arma::uvec p_index_av;
	int p_index=0;
	int pt_row=0;
	int ln_row=0;

	for(int x_index=0; x_index<left.size(); x_index++)  {
		if(left[x_index] == right[x_index])  {
			for(int q=0; q<qty[x_index]; q++)  {
				p_index_av=arma::find(ptime == tmean[x_index]);
				p_index = arma::as_scalar(p_index_av(0));
				time(pt_row) = left[x_index];
				ptpos(pt_row) = arma::as_scalar(ppp(p_index));
				ptrank(pt_row)= arma::as_scalar(adj_rank(p_index));
				ptwt(pt_row) = 1;
				ptime.shed_row(p_index);
				ppp.shed_row(p_index);
				adj_rank.shed_row(p_index);
				pt_row++;
			}
		}else{
			for(int q=0; q<qty[x_index]; q++)  {
				p_index_av=arma::find(ptime == tmean[x_index]);
				p_index = arma::as_scalar(p_index_av(0));
				t1(ln_row) = left[x_index];
				t2(ln_row) = right[x_index];
				lnpos(ln_row) = arma::as_scalar(ppp(p_index));
				lnrank(ln_row)= arma::as_scalar(adj_rank(p_index));
				lnwt(ln_row) = 1;
				ptime.shed_row(p_index);
				ppp.shed_row(p_index);
				adj_rank.shed_row(p_index);
				ln_row++;
			}
		}
	}
	if(pt_row<n) {
		time.shed_rows(pt_row, n-1);
		ptpos.shed_rows(pt_row, n-1);
		ptrank.shed_rows(pt_row, n-1);
		ptwt.shed_rows(pt_row, n-1);
	}
	if(ln_row<n) {
		t1.shed_rows(ln_row, n-1);
		t2.shed_rows(ln_row, n-1);
		lnpos.shed_rows(ln_row, n-1);
		lnrank.shed_rows(ln_row, n-1);
		lnwt.shed_rows(ln_row, n-1);
	}




	return List::create(Rcpp::Named("dpoints")=
		DataFrame::create(
			Rcpp::Named("time") = wrap( time),
			Rcpp::Named("ppp")=wrap(ptpos),
			Rcpp::Named("adj_rank")=wrap(ptrank),
			Rcpp::Named("weight")=wrap(ptwt)
		),
		Rcpp::Named("dlines") =
		DataFrame::create(
			Rcpp::Named("t1")=wrap(t1),
			Rcpp::Named("t2")=wrap(t2),
			Rcpp::Named("ppp")=wrap(lnpos),
			Rcpp::Named("adj_rank")=wrap(lnrank),
			Rcpp::Named("weight")=wrap(lnwt))
	);
}
