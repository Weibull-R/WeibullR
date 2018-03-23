// traditional header file content - class declaration
#include "WeibullR.h"
#include "MLEmodel.h"
#include <math.h>

    using namespace Rcpp ;

// class implementation
MLEmodel::MLEmodel(){}

MLEmodel::MLEmodel( SEXP arg1) {
	Rcpp::List L(arg1);
	time=Rcpp::as<arma::colvec>(L["fsdi"]);
	qty=Rcpp::as<arma::colvec>(L["q"]);
	N=L["N"];
// Provide a first non-sense element to front of time vector
// so that position math works when Nf is zero.
	time.insert_rows(0,1);
	qty.insert_rows(0,1);

	endf=N[0];
	ends=endf+N[1];
	endd=ends+N[2];
	endil=endd+N[3];
	endir=endil+N[3];

	if(N[0]>0)  {
		fail=time.rows(1,endf);
		nf=qty.rows(1,endf);
	}

	if(N[1]>0)  {
		susp=time.rows(endf+1,ends);
		ns=qty.rows(endf+1,ends);
	}

	if(N[2]>0)  {
		disc=time.rows(ends+1,endd);
		nd=qty.rows(ends+1,endd);
	}
	if(N[3]>0)  {
		left=time.rows(endd+1,endil);
		right=time.rows(endil+1,endir);
		ni=qty.rows(endd+1,endil);
	}

}


////************* Method LogLike ***************/
double MLEmodel::LogLike(arma::colvec par, int sign, int dist_num, double tz)  {

		double failcomp =0.0;
		double suscomp =0.0;
		double discomp =0.0;
		double intcomp =0.0;

		if(dist_num==1) {
			if(N[0]>0)  {
				for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dweibull(fail(i)-tz,par(0),par(1),1);
				}
			}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
//  Need to exempt any zero or negative values
					if(susp(i)-tz>0)  {
						suscomp=suscomp+ns(i)*R::pweibull(susp(i)-tz,par(0),par(1),0,1);
					}
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::pweibull(disc(i)-tz,par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
// if left(i) should become zero or less, then this data needs to be handled as if left censored on right(i)
					if(left(i)-tz>0)  {
						intcomp=intcomp+ni(i)*log(
						R::pweibull(left(i)-tz,par(0),par(1),0,0) -
						R::pweibull(right(i)-tz,par(0),par(1),0,0)
						);
					}else{
						intcomp=intcomp+ni(i)*log(1-R::pweibull(right(i)-tz,par(0),par(1),0,0));
					}
				}
			}

		}
		else if(dist_num==2)  {
				if(N[0]>0)  {
					for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dlnorm(fail(i)-tz,par(0),par(1),1);
					}
				}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
//  Need to exempt any zero or negative values
					if(susp(i)-tz>0)  {
						suscomp=suscomp+ns(i)*R::plnorm(susp(i)-tz,par(0),par(1),0,1);
					}
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::plnorm(disc(i)-tz,par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
// if left(i) should become zero or less, then this data needs to be handled as if left censored on right(i)
					if(left(i)-tz>0)  {
						intcomp=intcomp+ni(i)*log(
						R::plnorm(left(i)-tz,par(0),par(1),0,0) -
						R::plnorm(right(i)-tz,par(0),par(1),0,0)
						);
					}else{
						intcomp=intcomp+ni(i)*log(1-R::plnorm(right(i)-tz,par(0),par(1),0,0));
					}
				}
			}
		}
	return sign*(failcomp+suscomp+discomp+intcomp);
}


////************* Method tryLL ***************/
//
//  This method is specifically intended to be called by the getContourPt method of class contour
//  It will return zero given any negative parameter argument, or other condition producing a non-finite result.
//
//double MLEmodel::tryLL(arma::colvec par, int sign, int dist_num, double tz)  {
double MLEmodel::tryLL(arma::colvec par, int dist_num)  {
		double failcomp =0.0;
		double suscomp =0.0;
		double discomp =0.0;
		double intcomp =0.0;
		double value =0.0;
// elimintate needless processing with negative parameter
	if(par(0)>0 && par(1)>0) {
		if(dist_num==1) {
			if(N[0]>0)  {
				for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dweibull(fail(i),par(0),par(1),1);
				}
			}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
					suscomp=suscomp+ns(i)*R::pweibull(susp(i),par(0),par(1),0,1);
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::pweibull(disc(i),par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
					intcomp=intcomp+ni(i)*log(
					R::pweibull(left(i),par(0),par(1),0,0) -
					R::pweibull(right(i),par(0),par(1),0,0)
					);
				}
			}

		}
		else if(dist_num==2)  {
				if(N[0]>0)  {
					for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dlnorm(fail(i),par(0),par(1),1);
					}
				}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
					suscomp=suscomp+ns(i)*R::plnorm(susp(i),par(0),par(1),0,1);
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::plnorm(disc(i),par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
					intcomp=intcomp+ni(i)*log(
					R::plnorm(left(i),par(0),par(1),0,0) -
					R::plnorm(right(i),par(0),par(1),0,0)
					);
				}
			}
		}
		
		value = failcomp+suscomp+discomp+intcomp;
// removing this important stability feature, which is the only C++11 dependency in package
// due to a bug in mingw64 that prevented compiling of CallgetCCC2.cpp or CallgetPvalue.cppdue to 
// cc1plus.exe: out of memory allocating 65536 bytes
//https://sourceforge.net/p/mingw-w64/mailman/message/33182613/
//See if http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56746 applies to your
//case. If the answer is affirmative, passing -ftrack-macro-expansion=0
//will reduce memory usage.
//		if(!std::isfinite(value)) {
//			value=0.0;
//		}
	}
	return value;
}


////************* Method MLE_Simplex ***************/
//
//	An implementation of the Nelder-Meade simplex algorithm 
//  for optimizing the likelihood, specific to two parameters

SEXP MLEmodel::MLE_Simplex(SEXP arg1, arma::colvec vstart, double tz, int listout)  {
		Rcpp::List L(arg1);
		int dist_num=Rcpp::as<int>(L["dist_num"]);
// vstart is now brought in as an argument
		double limit=Rcpp::as<double>(L["limit"]);
		int maxit=Rcpp::as<int>(L["maxit"]);

// this algorithm is optimized specificity for the two parameter case
// variables to hold number of  parameters and number of vertices are vestigial
		int npar=2, k= npar+1;
 // coefficients for reflection, expansion, and contraction
		double ALPHA=1.0, BETA=0.5, GAMMA=2.0;
 // set the sign for minimization of negative LogLikelihood
		int sign= -1;

 // construct the initial simplex
		arma::mat v(2,3);
		for(int i=0;i<k;i++)  {
		v.col(i)=vstart;
		}
		v(0,1)=v(0,1)*1.2;
		v(1,2)=v(1,2)*1.2;

		arma::colvec funval(3);
		for(int i=0;i<k;i++)  {
		funval(i)=LogLike(v.col(i), sign, dist_num, tz );
		}



 // assign vertex order variables
		arma::uvec ndx=sort_index(funval);
		int vs=(int) ndx(0);
		int vh=(int) ndx(1);
		int vg=(int) ndx(2);
 // generate the initial error measure
		double P2avg=arma::as_scalar(mean(v.row(1)));
		double error=arma::as_scalar(sum(pow((v.row(1)-P2avg),2.0)/npar));
 // initialize the output dataframe for construction as a matrix
		arma::rowvec DFrow(4);
		DFrow(0)=v(0,vs);
		DFrow(1)=v(1,vs);
		DFrow(2)=funval(vs);
		DFrow(3)=error;
		arma::mat DF=DFrow;

 // initialization of variables used in the loop
		arma::colvec vm, vr, ve, vc;
		double fr, fe, fc;
		int loopcount=1;
		int warn=0;

 // This is the main loop for the minimization
	while(error>limit)  {
 // calculate the centroid
		vm=(v.col(vs)+v.col(vh))/2.0;
 // reflect vg to new vertex vr
		vr=(vm+ALPHA*(vm-v.col(vg)));
		fr=LogLike(vr, sign, dist_num, tz );
 // depending on success, save reflected vertex in place of vg
		if (fr < funval(vh) && fr >= funval(vs)) {
		v.col(vg)=vr;
		funval(vg)=fr;
		}else{
			if(funval(vs)<fr)  {

 // test for an outside contraction
			if (fr < funval(vg) && fr >= funval(vh))   {
				vc=vm+BETA*(vr-vm);
				}else{
 // this is an inside contraction
					vc=vm-BETA*(vr-vm);
				}
				fc=LogLike(vc, sign, dist_num, tz );
 // upon acceptance replace vg with contracton
				if (fc < funval(vg)) {
					v.col(vg) = vc;
					funval(vg)=fc;
				}else{
 // upon rare-to-never case shrink the simplex
					v.col(vh)=v.col(vs)+(v.col(vh)-v.col(vs))/2.0;
					v.col(vg)=v.col(vs)+(v.col(vg)-v.col(vs))/2.0;
 // This case results in two function calculations and simply replacing vh and vg points
					funval(vh)=LogLike(v.col(vh), sign, dist_num, tz );
					funval(vg)=LogLike(v.col(vg), sign, dist_num, tz );
				}
			}else{
// now we make an expansion
				ve=vm+GAMMA*(vr-vm);
				fe=LogLike(ve, sign, dist_num, tz );
 // store the better of reflection or expansion in place of vg
				if (fe < fr) {
					v.col(vg) =ve;
					funval(vg)=fe;
				}else{
					v.col(vg) = vr;
					funval(vg)=fr;
				}
			}
		}
 //re-assign vertex order variables
		ndx=sort_index(funval);
		vs=(int) ndx(0);
		vh=(int) ndx(1);
		vg=(int) ndx(2);

		P2avg=arma::as_scalar(mean(v.row(1)));
		error=arma::as_scalar(sum(pow((v.row(1)-P2avg),2.0)/npar));
		DFrow(0)=v(0,vs);
		DFrow(1)=v(1,vs);
		DFrow(2)=funval(vs);
		DFrow(3)=error;
		DF=join_cols(DF,DFrow);

		loopcount=loopcount+1;
		if(loopcount>maxit)  {
			warn=1;
			break;
		}


// then close main iteration loop
	}


		Rcpp::NumericVector outvec(4);
		if(dist_num==1) {
		outvec[0]=v(1,vs);
		outvec[1]=v(0,vs);
		}else{
		outvec[0]=v(0,vs);
		outvec[1]=v(1,vs);
		}
// note: multiplication by sign here assures positive log-likelihood is delivered
		outvec[2]=sign*funval(vs);
		outvec[3]=warn;

		if(listout==0)  {
			return outvec;
		}else {
			return List::create(outvec,wrap(DF));

		}

}


////************* Method dMaxLLdx ***************/
//
// A method for getting the secant slope used during 3rd parameter optimization.

	SEXP MLEmodel::dMaxLLdx( SEXP arg1, arma::colvec vstart, double tz)  {
		Rcpp::List L(arg1);
		int dist_num=Rcpp::as<int>(L["dist_num"]);
// vstart is now brought in as an argument
		double limit=Rcpp::as<double>(L["limit"]);
// int maxit=Rcpp::as<int>(L["maxit"]);

		Rcpp::NumericVector fit1(MLE_Simplex(arg1, vstart, tz, 0));
		Rcpp::NumericVector fit2(MLE_Simplex(arg1, vstart, tz-0.1*limit, 0));
		Rcpp::NumericVector dLLdx(3);
		 dLLdx[0]=(fit1[2]-fit2[2])/(0.1*limit);
// need to return current fit parameters for next start, switching position of parameters for weibull
		if(dist_num==1)  {
			dLLdx[1]=fit1[1];
			dLLdx[2]=fit1[0];
		}else{
			dLLdx[1]=fit1[0];
			dLLdx[2]=fit1[1];
		}
		return dLLdx;
	}


	// Exported Functions

	SEXP MLEloglike(SEXP arg1, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6)  {
		MLEmodel mymodel(arg1);
		arma::colvec par=Rcpp::as<arma::colvec>(arg3);
		int dist_num=Rcpp::as<int>(arg4);
		int sign=Rcpp::as<int>(arg5);
		double tz=Rcpp::as<double>(arg6);
		return wrap(mymodel.LogLike(par, sign, dist_num, tz));
	}

	SEXP MLEsimplex(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5)  {
		MLEmodel mymodel(arg1);
		arma::colvec vstart=Rcpp::as<arma::colvec>(arg3);
		double tz=Rcpp::as<double>(arg4);
		int listout=Rcpp::as<int>(arg5);
		return mymodel.MLE_Simplex(arg2, vstart, tz, listout);
		}

// used in secant method for determination of 3rd parameter optimization
	SEXP MLEdMaxLLdx(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4)  {
		MLEmodel mymodel(arg1);
		arma::colvec vstart=Rcpp::as<arma::colvec>(arg3);
		double tz=Rcpp::as<double>(arg4);
		return mymodel.dMaxLLdx(arg2, vstart, tz);
	}



