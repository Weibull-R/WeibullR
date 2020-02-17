#include "WeibullR.h"	


SEXP adjustedRank (SEXP arg1)			
{			
    using namespace Rcpp ;			
			
	arma::colvec event = Rcpp::as<arma::colvec>(arg1);		
	int N = event.n_rows;		
// abundance of caution regarding mixed-type math (int and double)			
	double Ndbl = (double) N;		
	int F = arma::as_scalar(sum(event));		
	arma::colvec adj_rank(N+1);		
	adj_rank.fill(0.0);		
// median_rank is only used to accumulate the elements			
// of the final return vector, so NumericVector is best choice			
	Rcpp::NumericVector order_num(F);		
	for(int i=1,j=0; i<N+1; i++)		
	{		
		double rr=(double)(N-i)+1.0;	
		if(event(i-1)>0)	
		{	
			adj_rank(i)= (rr*adj_rank(i-1)+Ndbl+1.0)/(rr+1.0);
			if(j<F) {
			order_num[j]=adj_rank(i);
			j++; }
		}	
		else	
		{	
			adj_rank(i)=adj_rank(i-1);
		}	
	}		
			
	return order_num;		
}
