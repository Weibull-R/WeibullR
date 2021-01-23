 /* LSLRmodel.cpp
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
 *
 * Author: David J. Silkworth
 *     Copyright (C) 2013-2021 OpenReliability.org
 */
#include "WeibullR.h"
#include "LSLRmodel.h"
#include <math.h>

    using namespace Rcpp ;

LSLRmodel::LSLRmodel( SEXP arg1) {
// unpack the arg list and set class members
//parlist<-list(fail=x$time, ppp=x$ppp, limit=limit, reg_order=reg_order, npar=npar, dist_num=dist_num )
	Rcpp::List L(arg1);
	data = L["fail"];
	positions = L["ppp"];
	regression_order = Rcpp::as<int>(L["reg_order"]);
	dist_num = Rcpp::as<int>(L["dist_num"]);
	npar = Rcpp::as<int>(L["npar"]);
	limit = Rcpp::as<double>(L["limit"]);
}

LSLRmodel::LSLRmodel(
// This constructor is used for 3-parameter code and the dR2dx method
	Rcpp::NumericVector _data,
	Rcpp::NumericVector _positions,
	int _regression_order,
	int _dist_num,
	int _npar,
	double _limit
) {	
	
	data = _data;
	positions = _positions;
	regression_order = _regression_order;
	dist_num = _dist_num;
	npar = _npar;
	limit = _limit;
}

std::vector<double> LSLRmodel::LSLRfit() {
	std::vector<double> outvec;

if(npar == 2) {
if(dist_num == 0) {
	
	int F=data.size();		
	arma::mat X(F,2);
	arma::colvec y(F);
		
// fill the arma objects		
    for(int i=0; i<F; i++)  {		
		X(i,0)=1.0;	
		X(i,1)=log(log(1/(1-positions[i])));	
		y(i)=log(data[i]);	
    }

// solve the linear equation and extract the R-square value using Armadillo's solve function	

	arma::colvec coef = arma::solve(X, y);
	arma::colvec res  = y - X*coef;	
	double Residual = arma::as_scalar(sum(square(res)));	
	double TVar = arma::as_scalar(sum(square(y-mean(y))));	
	double R2 = (TVar-Residual)/TVar;	


// Finally prepare a single vector with each coefficient and the variance (R2)
	if(regression_order == 0) {
	// this method applies the "X over Y" regression of the Weibull	
        outvec.push_back(exp(coef(0)));	
        outvec.push_back(1/coef(1));	
        outvec.push_back(R2);
	}else{
	// this method applies the "Y over X" regression of the Weibull		
        outvec.push_back(exp(-coef(0)/coef(1)));
        outvec.push_back(coef(1));
        outvec.push_back(R2);
	}
	
} // close dist_num == 0

if(dist_num == 1) {
	int F=data.size();		
	arma::mat X(F,2);
	arma::colvec y(F);
		
// fill the arma objects		
    for(int i=0; i<F; i++)  {		
		X(i,0)=1.0;	
		X(i,1)=Rf_qnorm5(positions[i],0.0,1.0,1,0);
		y(i)=log(data[i]);	
    }

// solve the linear equation and extract the R-square value using Armadillo's solve function	

	arma::colvec coef = arma::solve(X, y);
	arma::colvec res  = y - X*coef;	
	double Residual = arma::as_scalar(sum(square(res)));	
	double TVar = arma::as_scalar(sum(square(y-mean(y))));	
	double R2 = (TVar-Residual)/TVar;	


// Finally prepare a single vector with each coefficient and the variance (R2)
	if(regression_order == 0) {
	// this method applies the "X over Y" regression of the Lognormal
        outvec.push_back(coef(0));
        outvec.push_back(coef(1));
        outvec.push_back(R2);
	}else{
	// this method applies the "Y over X" regression of the Lognormal	
        outvec.push_back(-coef(0)/coef(1));
        outvec.push_back(1/coef(1));
        outvec.push_back(R2);
	}	
} // close dist_num == 1	
} // close npar == 2

if(npar == 3) {
// data, positions, regression_order, dist_num and limit are known by  the present class object
	int F=data.size();	
	int warning=0;	
	int maxit=100;	
	double DL = limit;	
		
 // just get the first failure, they were required to have been sorted		
	double C1=data[0];	
		
 // Tao Pang's original variable labels from FORTRAN are used where possible		
 // initial step is based on limit*10,000		
 // aruguments to pow must be double, not int for sparc build		
	double DX=DL*pow((double) 10.0,(double) 4.0);	
	double X0=0.0;	
	int istep=0;	
	double X1=X0+DX;	
	if(X1>C1) {X1=X0+0.9*(C1-X0);}	
		
//	double FX0=dR2dx(X0,fail,ppp,DL, LSLR2p);	
//	double FX1=dR2dx(X1,fail,ppp,DL, LSLR2p);	
	double FX0=dR2dx(X0);	
	double FX1=dR2dx(X1);	
 // FX1 will contain slope sign information to be used only one time to find X2		
	double D=fabs(FX1-FX0);	
	double X2=X1+fabs(X1-X0)*FX1/D;	
	if(X2>C1) {X2=X1+0.9*(C1-X1);}	
	X0=X1;	
	X1=X2;	
		
 // This is the start of the main loop		
	while(fabs(DX)>DL && istep<maxit)  {	
	FX0=FX1;	
//	 FX1=dR2dx(X1,fail,ppp,DL, LSLR2p);	
	FX1=dR2dx(X1);	
		
	if(FX1!=FX1)  {	
		FX1=FX0;
		warning=1;
		break;
	}	
		
 // FX1 will contain slope sign information to be used only one time to find X2		
	D=fabs(FX1-FX0);	
	X2=X1+fabs(X1-X0)*FX1/D;	
	if(X2>C1) {X2=X1+0.9*(C1-X1);}	
	X0=X1;	
	X1=X2;	
	DX=X1-X0;	
	istep=istep+1;	
	}	
		
	Rcpp::NumericVector mdata(F);	
	for(int i=0; i<F; i++) {mdata[i]=data[i]-X0;}	
	std::unique_ptr<LSLRmodel> final(new LSLRmodel(	
		mdata, positions, regression_order,
		dist_num, npar-1, limit)
		);
	std::vector<double> finalfit = final->LSLRfit();	
		
	outvec.push_back(finalfit[0]);	
	outvec.push_back(finalfit[1]);	
	outvec.push_back(X0);	
	outvec.push_back(finalfit[2]);	
	outvec.push_back(warning);	
	
}  // close npar == 3
	
//	return Rcpp::wrap(X);
	return outvec;
} 

double LSLRmodel::dR2dx(double X) {	
// This method is a workhorse in 3p optimization it creates new, temporary, LSLRmodel objects 
// using modified data for 2p fitting, establishing the slope of the secant
	int F=data.size();	
	Rcpp::NumericVector mdata(F);	
		
	for(int i=0; i<F; i++) {mdata[i]=data[i]-X;}	
	std::unique_ptr<LSLRmodel> LM1(new LSLRmodel(	
		mdata, positions, regression_order,
		dist_num, npar-1, limit)
		);
	std::vector<double> fit1 = LM1->LSLRfit();	
		
	for(int i=0; i<F; i++) {mdata[i]=data[i]-(X-0.1*limit);}	
	std::unique_ptr<LSLRmodel> LM2(new LSLRmodel(	
		mdata, positions, regression_order,
		dist_num, npar-1, limit)
		);
	std::vector<double> fit2 = LM2->LSLRfit();	
	double slope=(fit1[2]-fit2[2])/(0.1*limit);	
		
	return slope;	
}		


