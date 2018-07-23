 /* LSLR.cpp
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
 * This collectioin of functions implement least squares linear regression (LSLR).
 * The covered distributions are Weibull, lognormal, and Gumbel (extreme value type 1) due to wide use
 * in reliability analysis. Probability plotting positions are a required argument vector
 * therefore permitting	any methodology for point estimation to be employed in pre-processing.
 * Alternate minimization functions for X or Y axis variances are provided.
 *
 * Two arguments are required: a vector of data values (often recorded as time), and an equally
 * corresponding vector of plotting positions.  These are accepted without checking of any kind.

 * These functions are consistent with "The New Weibull Handbook, Fifth Edition" and have been checked
 * against SuperSMITH and Mintab commercial software packages.
 *
 * Author: Jacob Ormerod
 *     Copyright (C) 2013-2014 OpenReliability.org
 */

#include "WeibullR.h"
#include <math.h>

    using namespace Rcpp ;

SEXP LSLRw2pXonY (SEXP arg1, SEXP arg2)
{
 //   using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
	X(i,1)=log(log(1/(1-position[i])));
	y(i)=log(fail[i]);
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
        outvec[0]=exp(coef(0));
        outvec[1]=1/coef(1);
        outvec[2]=R2;

        return outvec;
        }


SEXP LSLRw2pYonX (SEXP arg1, SEXP arg2)
{
 //   using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
	X(i,1)=log(fail[i]);
	y(i)=log(log(1/(1-position[i])));
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
        outvec[0]=exp(-coef(0)/coef(1));
        outvec[1]=coef(1);
        outvec[2]=R2;

        return outvec;
        }



SEXP LSLRg2pXonY (SEXP arg1, SEXP arg2)
{
 //   using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
	X(i,1)=log(log(1/(1-position[i])));
	y(i)=log(fail[i]);
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
         outvec[0]=exp(coef(0));
        outvec[1]=1/coef(1);
        outvec[2]=R2;

        return outvec;
        }



SEXP LSLRg2pYonX (SEXP arg1, SEXP arg2)
{
//    using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
    X(i,1)=log(fail[i]);
    y(i)=log(log(1/(1-position[i])));
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
        outvec[0]=exp(-coef(0)/coef(1));
        outvec[1]=coef(1);
        outvec[2]=R2;

        return outvec;
        }



        SEXP LSLRln2pXonY (SEXP arg1, SEXP arg2)
{
 //   using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
	X(i,1)=Rf_qnorm5(position[i],0.0,1.0,1,0);
	y(i)=log(fail[i]);
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
        outvec[0]=coef(0);
        outvec[1]=coef(1);
        outvec[2]=R2;

        return outvec;
        }



    SEXP LSLRln2pYonX (SEXP arg1, SEXP arg2)
    {
 //       using namespace Rcpp ;

        Rcpp::NumericVector fail(arg1);
        Rcpp::NumericVector position(arg2);


        int F=fail.size();
// declare the arma objects with n_rows = F for bounds checking
        arma::mat X(F,2);
        arma::colvec y(F);
// fill the arma objects
        for(int i=0; i<F; i++)  {
	X(i,0)=1.0;
	X(i,1)=log(fail[i]);
	y(i)=Rf_qnorm5(position[i],0.0,1.0,1,0);
        }
        arma::colvec coef, res;
        double Residual, TVar, R2;

// solve the linear equation and extract the R-square value using Armadillo's solve function
// this method applies the "X over Y" regression of the Weibull
        coef = arma::solve(X, y);
        res  = y - X*coef;
        Residual = arma::as_scalar(sum(square(res)));
        TVar = arma::as_scalar(sum(square(y-mean(y))));
        R2 = (TVar-Residual)/TVar;
// Finally prepare a single vector with each coefficient and the variance (R2)
        Rcpp::NumericVector outvec(3);
         outvec[0]= -coef(0)/coef(1);
        outvec[1]=1/coef(1);
        outvec[2]=R2;

        return outvec;
        }

// numeric 2-point numerical derivative for secant in discrete Newtonian root finding method
double dR2dx(double X, NumericVector &data, NumericVector &ppp,double limit, SEXP (*LSLR2p)(SEXP, SEXP) )  {

	int F=data.size();
	Rcpp::NumericVector mdata(F);
	for(int i=0; i<F; i++) {mdata[i]=data[i]-X;}
	Rcpp::NumericVector fit1=LSLR2p(mdata,ppp);
	for(int i=0; i<F; i++) {mdata[i]=data[i]-(X-0.1*limit);}
	Rcpp::NumericVector fit2=LSLR2p(mdata,ppp);
	double slope=(fit1[2]-fit2[2])/(0.1*limit);

	return slope;
}

SEXP LSLR3p(SEXP arg1, SEXP arg2, SEXP arg3, SEXP (*LSLR2p)(SEXP, SEXP) )
{
	        Rcpp::NumericVector fail(arg1);
	        Rcpp::NumericVector ppp(arg2);
	double DL=as<double>(arg3);
	        int F=fail.size();
			int warning=0;



	int maxit=100;
 // just get the first failure, they were required to have been sorted
	double C1=fail[0];








 // Tao Pang's original variable labels from FORTRAN are used where possible
 // initial step is based on limit*10,000
 // aruguments to pow must be double, not int for sparc build
	double DX=DL*pow((double) 10.0,(double) 4.0);
	double X0=0.0;
	int istep=0;
	double X1=X0+DX;
	if(X1>C1) {X1=X0+0.9*(C1-X0);}

	double FX0=dR2dx(X0,fail,ppp,DL, LSLR2p);
	double FX1=dR2dx(X1,fail,ppp,DL, LSLR2p);
 // FX1 will contain slope sign information to be used only one time to find X2
	double D=fabs(FX1-FX0);
	double X2=X1+fabs(X1-X0)*FX1/D;
	if(X2>C1) {X2=X1+0.9*(C1-X1);}
	X0=X1;
	X1=X2;
 // development diagnostic code
 //	arma::rowvec DFrow(4);
 //	DFrow(0)=istep;
 //	DFrow(1)=X0;
 //	DFrow(2)=DX;
 //	DFrow(3)=FX1;
 //	arma::mat DF;
 //	DF=DFrow;


 // This is the start of the main loop
//	while(fabs(DX)>DL&& istep<maxit&&warning==0)  {
	while(fabs(DX)>DL && istep<maxit)  {
	FX0=FX1;
	 FX1=dR2dx(X1,fail,ppp,DL, LSLR2p);

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
 // development diagnostic code
 //	DFrow(0)=istep;
 //	DFrow(1)=X0;
 //	DFrow(2)=DX;
 //	DFrow(3)=FX1;
 //	DF=join_cols(DF,DFrow);
	}

	Rcpp::NumericVector mdata(F);
	for(int i=0; i<F; i++) {mdata[i]=fail[i]-X0;}
	Rcpp::NumericVector finalfit=LSLR2p(mdata,ppp);
	Rcpp::NumericVector outvec(5);
	outvec[0]=finalfit[0];
	outvec[1]=finalfit[1];
	outvec[2]=X0;
	outvec[3]=finalfit[2];
	outvec[4]=warning;

	return outvec;


	}

SEXP LSLR(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4)
{

int casenum =as<int>(arg4);

Rcpp::NumericVector outvec2p(3);
Rcpp::NumericVector outvec3p(4);



switch(casenum) {
case 0 :
outvec2p=LSLRw2pXonY(arg1, arg2);
return outvec2p;
break;

case 1 :
outvec2p=LSLRw2pYonX(arg1, arg2);
return outvec2p;
break;

case 2 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRw2pXonY);
return outvec3p;
break;

case 3 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRw2pYonX);
return outvec3p;
break;

case 4 :
outvec2p=LSLRln2pXonY(arg1, arg2);
return outvec2p;
break;

case 5 :
outvec2p=LSLRln2pYonX(arg1, arg2);
return outvec2p;
break;

case 6 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRln2pXonY);
return outvec3p;
break;

case 7 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRln2pYonX);
return outvec3p;
break;

case 8 :
outvec2p=LSLRg2pXonY(arg1, arg2);
return outvec2p;
break;

case 9 :
outvec2p=LSLRg2pYonX(arg1, arg2);
return outvec2p;
break;

case 10 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRg2pXonY);
return outvec3p;
break;

case 11 :
outvec3p=LSLR3p(arg1, arg2, arg3, &LSLRg2pYonX);
return outvec3p;
break;

default:
outvec2p=LSLRw2pXonY(arg1, arg2);
return outvec2p;
break;
}

}







