/* pivotalMC.cpp					
 *					
 * This program is free software; you can redistribute it and/or modify it					
 * under the terms of the GNU General Public License as published by the					
 * Free Software Foundation; either version 3, or (at your option) any					
 * later version.					
 *					
 * These functions are distributed in the hope that they will be useful,					
 * but WITHOUT ANY WARRANTY; without even the implied warranty of					
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the					
 * GNU General Public License for more details.					
 *					
 *  You should have received a copy of the GNU General Public License					
 *  along with this program; if not, a copy is available at					
 *  http://www.r-project.org/Licenses/					
 *					
 * This collection of functions implementing a fast re-sampling engine for least squares linear regression
 * models. The covered distributions are Weibull, log-normal, and Gumbel (extreme value type 1) due to
 * wide use in reliability analysis. Probability plotting positions are a required argument vector 
 * therefore permitting	any methodology for point estimation to be employed in pre-processing.
 * Alternate minimization functions for X or Y axis variances are provided.
 *					
 * These functions are consistent with "The New Weibull Handbook, Fifth Edition" and have been checked
 * against SuperSMITH and Mintab commercial software packages.					
 *									
 * Author: Jacob Ormerod					
 *     Copyright (C) 2013-2014 OpenReliability.org					
 */	
#include "WeibullR.h"	
			
    using namespace Rcpp ;
	
SEXP pivotalMCw2pXonY(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){				
				
 //	src1<-'			
// ppp Must be determined in calling code, choices include several point estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();
	
// need the event vector for masking			
// length of event vector will identify the total number of occurrences.			
	Rcpp::NumericVector event(arg2);		
	int N=event.size();		
			
// establish output to be provided and prepare case records			
	Rcpp::NumericVector SimControl(arg3);		
	double R2test= SimControl[0];		
	double  CItest=SimControl[1];		
	int prrout=0;		
	int pivout=0;		
			
	unsigned int S = as<unsigned int>(arg4);		
	unsigned int Spct = S/100;		
	int seed = as<int>(arg5);		
// get the descriptive quantiles for pivotals			
	Rcpp::NumericVector dq(arg6);		
	int ndq = dq.size();		
			
// variables to control a progress output			
	bool ProgRpt = as<bool>(arg7);		
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Beta or Eta				
// but prr output would vary according to plotting position differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Beta and Eta				
// but should ultimately be able to be transformed so that the median will conform to Eta=Beta=1				
	double Eta = SimControl[2];			
	double Beta = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];
	if(seed>0) {
	SetSeed(seed);
	}
				
// Fill a matrix appropriate for application of linear fit XonY				
 // explicit loop to fill arma object from Rcpp object				
	arma::mat X(F,2);			
	X.fill(1.0);			
	for (int i=0; i<F; i++) {			
		X(i,1)=log(log(1/(1-ppp[i])));		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=log(log(1/(1-dq[i])));		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);				
				
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rweibull(N, Beta, Eta));		
		y = arma::sort(y);
		for(int j=0,k=0; j<N; j++)  {		
			if(event[j]==0)  {y.shed_row(j-k);}	
			k++;	
		}
		y=log(y);				
				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "X over Y" regression of the Weibull				
		coef = arma::solve(X, y);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = y - X*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(y-mean(y))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/b_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
				
		qpiv.row(i)=arma::trans((CBq-coef(0))/coef(1));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
		            return wrap(0.0);		
				
	}			
 //	 '			
}

	SEXP pivotalMCln2pXonY(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){				
			
 //	src2<-'			
// ppp Must be determined in calling code, choices include exact and Benard's estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();	
	
// need the event vector for masking			
// length of event vector will identify the total number of occurrences.			
	Rcpp::NumericVector event(arg2);		
	int N=event.size();		
			
// establish output to be provided and prepare case records			
	Rcpp::NumericVector SimControl(arg3);		
	double R2test= SimControl[0];		
	double  CItest=SimControl[1];		
	int prrout=0;		
	int pivout=0;		
			
	unsigned int S = as<unsigned int>(arg4);		
	unsigned int Spct = S/100;		
	int seed = as<int>(arg5);		
// get the descriptive quantiles for pivotals			
	Rcpp::NumericVector dq(arg6);		
	int ndq = dq.size();		
			
// variables to control a progress output			
	bool ProgRpt = as<bool>(arg7);		
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Mulog orSigmalog				
// but prr output would vary according to rank position differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Mulog and Sigmalog				
// but should ultimately be able to be transformed so that the median will conform to Mulog=0;Sigmalog=1				
	double Mulog = SimControl[2];			
	double Sigmalog = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];			
	if(seed>0) {
	SetSeed(seed);
	}			
				
// Fill a matrix appropriate for application of linear fit XonY				
				
	arma::mat X(F,2);			
	X.fill(1.0);			
	for (int i=0; i<F; i++) {			
		X(i,1)=Rf_qnorm5(ppp[i],0.0,1.0,1,0);		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=Rf_qnorm5(dq[i],0.0,1.0,1,0);		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);			
				
				
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rnorm(N, Mulog, Sigmalog));		
		y = arma::sort(y);
		for(int j=0,k=0; j<N; j++)  {	
			if(event[j]==0)  {y.shed_row(j-k);}
			k++;
		}	
		
//  for lognormal the data at y is already log(y)					
//y=log(y);				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "X over Y" regression of the Weibull				
		coef = arma::solve(X, y);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = y - X*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(y-mean(y))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/s_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
// need to confirm coeficient use here (x still equals (y-b)/m, if y=mx+b)				
		qpiv.row(i)=arma::trans((CBq-coef(0))/coef(1));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
	            return wrap(0.0);			
				
	}			
 //	 '			
}				


SEXP pivotalMCw2pYonX(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){			
		
 //	src3<-'			
// ppp Must be determined in calling code, choices include exact and Benard's estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();
	
// need the event vector for masking			
// length of event vector will identify the total number of occurrences.			
	Rcpp::NumericVector event(arg2);		
	int N=event.size();		
			
// establish output to be provided and prepare case records			
	Rcpp::NumericVector SimControl(arg3);		
	double R2test= SimControl[0];		
	double  CItest=SimControl[1];		
	int prrout=0;		
	int pivout=0;		
			
	unsigned int S = as<unsigned int>(arg4);		
	unsigned int Spct = S/100;		
	int seed = as<int>(arg5);		
// get the descriptive quantiles for pivotals			
	Rcpp::NumericVector dq(arg6);		
	int ndq = dq.size();		
			
// variables to control a progress output			
	bool ProgRpt = as<bool>(arg7);		
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Beta or Eta				
// but prr output would vary according to plotting position differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Beta and Eta				
// but should ultimately be able to be transformed so that the median will conform to Eta=Beta=1				
	double Eta = SimControl[2];			
	double Beta = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];			
	if(seed>0) {
	SetSeed(seed);
	}			
				
// Fill a matrix appropriate for application of linear fit XonY				
				
	arma::colvec x(F);			
				
	for (int i=0; i<F; i++) {			
		x(i)=log(log(1/(1-ppp[i])));		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=log(log(1/(1-dq[i])));		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);			
	arma::mat Y(F,2);			
	Y.fill(1.0);			
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rweibull(N, Beta, Eta));		
		y = arma::sort(y);
		for(int j=0,k=0; j<N; j++)  {	
			if(event[j]==0)  {y.shed_row(j-k);}
			k++;
		}		
		y=log(y);		
		for(int j=0; j<F; j++) {		
			Y(j,1)=y(j);	
		}		
				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "YonX" regression of the Weibull				
		coef = arma::solve(Y, x);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = x - Y*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(x-mean(x))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/b_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
// need to confirm coeficient use here with x and y reversed x=my+b)				
		qpiv.row(i)=arma::trans(coef(1)*CBq+coef(0));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
		            return wrap(0.0);		
				
	}			
 //	 '			
}				


SEXP pivotalMCln2pYonX(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){				
				
 //	src4<-'			
// ppp Must be determined in calling code, choices include exact and Benard's estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();			
				
// need the event vector for masking		
// length of event vector will identify the total number of occurrences.		
	Rcpp::NumericVector event(arg2);	
	int N=event.size();	
		
// establish output to be provided and prepare case records		
	Rcpp::NumericVector SimControl(arg3);	
	double R2test= SimControl[0];	
	double  CItest=SimControl[1];	
	int prrout=0;	
	int pivout=0;	
		
	unsigned int S = as<unsigned int>(arg4);	
	unsigned int Spct = S/100;	
	int seed = as<int>(arg5);	
// get the descriptive quantiles for pivotals		
	Rcpp::NumericVector dq(arg6);	
	int ndq = dq.size();	
		
// variables to control a progress output		
	bool ProgRpt = as<bool>(arg7);	
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Mulog orSigmalog				
// but prr output would vary according to mrank differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Mulog and Sigmalog				
// but should ultimately be able to be transformed so that the median will conform to Mulog=0;Sigmalog=1				
	double Mulog = SimControl[2];			
	double Sigmalog = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];			
	if(seed>0) {
	SetSeed(seed);
	}			
				
// Fill a matrix appropriate for application of linear fit XonY				
				
	arma::colvec x(F);			
				
	for (int i=0; i<F; i++) {			
		x(i)=Rf_qnorm5(ppp[i],0.0,1.0,1,0);		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=Rf_qnorm5(dq[i],0.0,1.0,1,0);		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);			
	arma::mat Y(F,2);			
	Y.fill(1.0);			
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rnorm(N, Mulog, Sigmalog));		
		y = arma::sort(y);	
		for(int j=0,k=0; j<N; j++)  {	
			if(event[j]==0)  {y.shed_row(j-k);}
			k++;
		}		
//  for lognormal the data at y is already log(y)				
		for(int j=0; j<F; j++) {		
			Y(j,1)=y(j);	
		}		
				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "YonX" regression of the Weibull				
		coef = arma::solve(Y, x);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = x - Y*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(x-mean(x))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/s_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
// need to confirm coefficient use here with x and y reversed x=my+b)				
		qpiv.row(i)=arma::trans(coef(1)*CBq+coef(0));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
	            return wrap(0.0);			
				
	}			
 //	 '			
}


SEXP pivotalMCg2pXonY(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){				
			
 //	src5<-'			
// ppp Must be determined in calling code, choices include exact and Benard's estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();			
				
// need the event vector for masking	
// length of event vector will identify the total number of occurrences.	
	Rcpp::NumericVector event(arg2);
	int N=event.size();
	
// establish output to be provided and prepare case records	
	Rcpp::NumericVector SimControl(arg3);
	double R2test= SimControl[0];
	double  CItest=SimControl[1];
	int prrout=0;
	int pivout=0;
	
	unsigned int S = as<unsigned int>(arg4);
	unsigned int Spct = S/100;
	int seed = as<int>(arg5);
// get the descriptive quantiles for pivotals	
	Rcpp::NumericVector dq(arg6);
	int ndq = dq.size();
	
// variables to control a progress output	
	bool ProgRpt = as<bool>(arg7);
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Beta or Eta				
// but prr output would vary according to plotting position differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Beta and Eta				
// but should ultimately be able to be transformed so that the median will conform to Eta=Beta=1				
	double Eta = SimControl[2];			
	double Beta = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];			
	if(seed>0) {
	SetSeed(seed);
	}			
				
// Fill a matrix appropriate for application of linear fit XonY				
 // explicit loop to fill arma object from Rcpp object				
	arma::mat X(F,2);			
	X.fill(1.0);			
	for (int i=0; i<F; i++) {			
		X(i,1)=log(log(1/(1-ppp[i])));		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=log(log(1/(1-dq[i])));		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);			
				
				
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rweibull(N, Beta, Eta));		
		y = arma::sort(y);
		for(int j=0,k=0; j<N; j++)  {	
			if(event[j]==0)  {y.shed_row(j-k);}
			k++;
		}			
//  for gumbel the data at y is already log(y)				
 // y=log(y);				
				
				
				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "X over Y" regression of the Weibull				
		coef = arma::solve(X, y);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = y - X*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(y-mean(y))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/b_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
				
		qpiv.row(i)=arma::trans((CBq-coef(0))/coef(1));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
		            return wrap(0.0);		
				
	}			
 //	 '			
}				


SEXP pivotalMCg2pYonX(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7){				
			
 //	src6<-'			
// ppp Must be determined in calling code, choices include exact and Benard's estimation methods				
// quantity of ppp will identify the number of complete failures.				
	Rcpp::NumericVector ppp(arg1);			
	int F=ppp.size();			
				
// need the event vector for masking		
// length of event vector will identify the total number of occurrences.		
	Rcpp::NumericVector event(arg2);	
	int N=event.size();	
		
// establish output to be provided and prepare case records		
	Rcpp::NumericVector SimControl(arg3);	
	double R2test= SimControl[0];	
	double  CItest=SimControl[1];	
	int prrout=0;	
	int pivout=0;	
		
	unsigned int S = as<unsigned int>(arg4);	
	unsigned int Spct = S/100;	
	int seed = as<int>(arg5);	
// get the descriptive quantiles for pivotals		
	Rcpp::NumericVector dq(arg6);	
	int ndq = dq.size();	
		
// variables to control a progress output		
	bool ProgRpt = as<bool>(arg7);	
			
int ProgPct=0;				
int LastPct=0;				
				
// testing has suggested that prr output does not depend on value of Beta or Eta				
// but prr output would vary according to plotting position differences (as with treatment of censoring)				
// The pivotal quantities for confidence bounds will be effected by the sampled Beta and Eta				
// but should ultimately be able to be transformed so that the median will conform to Eta=Beta=1				
	double Eta = SimControl[2];			
	double Beta = SimControl[3];			
				
	RNGScope scope;			
	Environment base("package:base");			
	Function SetSeed = base["set.seed"];			
	if(seed>0) {
	SetSeed(seed);
	}			
				
// Fill a matrix appropriate for application of linear fit XonY				
				
	arma::colvec x(F);			
				
	for (int i=0; i<F; i++) {			
		x(i)=log(log(1/(1-ppp[i])));		
	}			
				
// establish the quantiles for evaluation of confidence bounds				
	arma::colvec CBq(ndq);			
	for(int i = 0; i<ndq; i++)  {			
		CBq(i)=log(log(1/(1-dq[i])));		
	}			
				
// The main iteration loop to generate the population of R-square values				
// for random samples of the 2-parameter weibull distribution				
	arma::colvec y, coef, res;			
	double Residual, TVar, pvalue, CCC2;			
	arma::colvec R2(S);			
	arma::mat qpiv(S,ndq);			
	arma::mat Y(F,2);			
	Y.fill(1.0);			
	for(unsigned int i=0; i<S; i++)  {			
		y = Rcpp::as<arma::colvec>(rweibull(N, Beta, Eta));		
		y = arma::sort(y);
		for(int j=0,k=0; j<N; j++)  {	
			if(event[j]==0)  {y.shed_row(j-k);}
			k++;
		}		
//  for gumbel the data at y is already log(y)				
		for(int j=0; j<F; j++) {		
			Y(j,1)=y(j);	
		}		
				
// solve the linear equation and extract the R-square value using Armadillo's solve function				
// this method applies the "YonX" regression of the Weibull				
		coef = arma::solve(Y, x);		
				
// the prr vector is built here only if called for by output control				
	if(R2test>0.0) {			
		res  = x - Y*coef;		
		Residual = arma::as_scalar(sum(square(res)));		
		TVar = arma::as_scalar(sum(square(x-mean(x))));		
		R2(i) = (TVar-Residual)/TVar;		
	}			
				
// Here the pivotals for confidence bounds are built				
// note that this pivotal is composed of (yp-u_hat)/b_hat, which is negative of Lawless' pivotal				
// the pivotals matrix is only built if called for by output control				
	if(CItest>0.0)  {			
// need to confirm coefficient use here with x and y reversed x=my+b)				
		qpiv.row(i)=arma::trans(coef(1)*CBq+coef(0));		
	}			
				
				
	if(ProgRpt) {			
// Progress report				
ProgPct = (i+1)/Spct;				
if(ProgPct > LastPct)  {				
Rprintf("%3d%% completion",(i+1)/Spct);				
Rprintf("\r");				
R_FlushConsole();				
R_ProcessEvents();				
}				
LastPct = ProgPct;				
		}		
//close main loop				
	}			
				
	int LCB=0;			
	arma::rowvec LBpiv(ndq);			
	arma::rowvec HBpiv(ndq);			
	arma::rowvec median(ndq);			
				
				
				
				
// process the prr vector according to output control				
	if(R2test>0.0) {			
	prrout = 1;			
	R2=arma::sort(R2);			
	if(R2test< 1.0) {			
		prrout=2;		
		arma::colvec Absolute(1);		
		Absolute(0)=1.0;		
		R2=join_cols(R2,Absolute);		
// pve_u is an integer representation  of the percentile of R-square				
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");		
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation				
		pvalue= (double) (pvalue_u(0)) /S*100;		
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10				
		CCC2= (double) R2(S/10-1);		
				
				
	}}			
				
	if(CItest>0.0) {			
		pivout=1;		
		for(int i=0; i<ndq; i++)  {		
			qpiv.col(i)=arma::sort(qpiv.col(i));	
		}		
	if(CItest< 1.0) {			
		pivout=2;		
		LCB=(int) S*(1-CItest)/2;		
		LBpiv=qpiv.row(LCB-1);		
		HBpiv=qpiv.row(S-LCB-1);		
		median=qpiv.row(S/2-1);		
				
	}}			
				
	int outputcase=prrout+4*pivout;			
	switch(outputcase)			
	{			
		case 1:		
			return wrap(R2);	
			break;	
		case 2:		
		return DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 );	
			break;	
				
	case 4:			
		return wrap(qpiv);		
		break;		
// this is the unlikely case that both extended output objects are called for				
	case 5:			
		return List::create(		
		Rcpp::Named("prr")=wrap(R2),		
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 6:			
		return List::create(		
		Rcpp::Named("prrCCC2")=		
		DataFrame::create(		
			Rcpp::Named("Pvalue")=pvalue,	
			Rcpp::Named("CCC2")=CCC2 ),	
		Rcpp::Named("pivotals")=wrap(qpiv) );		
		break;		
				
	case 8:			
		return DataFrame::create(		
			Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
			Rcpp::Named("Median")=wrap(arma::trans(median)),	
			Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );	
		break;		
				
	case 9:			
		return List::create(		
			Rcpp::Named("prr")=wrap(R2),	
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
	case 10:			
		return List::create(		
			Rcpp::Named("prrCCC2")=	
			 DataFrame::create(	
				Rcpp::Named("Pvalue")=pvalue,
				Rcpp::Named("CCC2")=CCC2 ),
			Rcpp::Named("pivotals")=	
			DataFrame::create(	
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),
				Rcpp::Named("Median")=wrap(arma::trans(median)),
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )
			);	
		break;		
				
		default:		
		            return wrap(0.0);		
				
	}			
 //	 '			
}

//result<-.Call("pivotalMC", x$ppp, event, c(R2,CI,P1,P2), S, seed, unrel, ProgRpt, casenum , package="WeibullR")				
			
SEXP pivotalMC(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6,SEXP arg7,SEXP arg8)		
{				

int casenum =as<int>(arg8);		
		
switch(casenum) {				
case 0 :				
return pivotalMCw2pXonY(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
break;				
				
case 1 :				
return pivotalMCw2pYonX(arg1, arg2, arg3, arg4, arg5, arg6, arg7);		
break;				
				
case 2 :								
return pivotalMCln2pXonY(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
break;				
				
case 3 :								
return pivotalMCln2pYonX(arg1, arg2, arg3, arg4, arg5, arg6, arg7);		
break;				
				
case 4 :				
return pivotalMCg2pXonY(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
break;				
				
case 5 :				
return pivotalMCg2pYonX(arg1, arg2, arg3, arg4, arg5, arg6, arg7);		
break;					
				
default:				
return pivotalMCw2pXonY(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
break;				
}				
				
}	
