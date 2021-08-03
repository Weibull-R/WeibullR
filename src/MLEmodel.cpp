// traditional header file content - class declaration
#include "WeibullR.h"
#include "MLEmodel.h"
#include <math.h>
//#include <cmath.h> // for abs (cmath.h may need linker option -lstdc++ for GCC)

    using namespace Rcpp ;

// class implementation
MLEmodel::MLEmodel(){}

MLEmodel::MLEmodel( SEXP arg1) {
	Rcpp::List L(arg1);
	arma::colvec time=Rcpp::as<arma::colvec>(L["fsdi"]);
	arma::colvec qty=Rcpp::as<arma::colvec>(L["q"]);
	N=L["N"];
	dist_num=L["dist_num"];
// Provide a first non-sense element to front of time vector
// so that position math works when Nf is zero.
	time.insert_rows(0,1);
	qty.insert_rows(0,1);

	int endf=N[0];
	int ends=endf+N[1];
	int endd=ends+N[2];
	int endil=endd+N[3];
	int endir=endil+N[3];

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

void MLEmodel::setSimplexLimit(double _limit) {
	limit = _limit;
}
void MLEmodel::setSimplexMaxit(int _maxit) {
	maxit = _maxit;
}

int MLEmodel::getDistNum() {
	return dist_num;
}


////************* Method LogLike ***************/
double MLEmodel::LogLike(arma::colvec par, int sign, double tz)  {

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
//double MLEmodel::tryLL(arma::colvec par, int dist_num)  {
double MLEmodel::tryLL(arma::colvec par)  {
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
// The following code restores stability should NaN be encountered
// NaN has the property of reporting false for any equality comparison
		if(value!=value) {
			value=0.0;
		}
	}
	return value;
}

////************* Method tryLL2 ***************/
//
//  This method is specifically intended to be called by the MLE_Simplex method of class MLEmodel
//  It will return zero given any negative parameter argument, or other condition producing a non-finite result.
//
//  tryLL2 takes a tz argument as double. I suppose this could simply be a polymorphism of tryLL
//  
//  as with tryLL there is no sign argument. The sign value is 'baked-in' for a minimization with the MLE_Simplex calls
//double MLEmodel::tryLL(arma::colvec par, int sign, int dist_num, double tz)  {
//double MLEmodel::tryLL(arma::colvec par, int dist_num)  {
double MLEmodel::tryLL2(arma::colvec par, double tz)  {
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
// The following code restores stability should NaN be encountered
// NaN has the property of reporting false for any equality comparison
		if(value!=value) {
			value=0.0;
		}
	}
	return value;
}

////************* Method MLE_Simplex ***************/
//
//	An implementation of the Nelder-Meade simplex algorithm 
//  for optimizing the likelihood, specific to two parameters

SEXP MLEmodel::MLE_Simplex(arma::colvec vstart, double tz, int listout)  {
// dist_num, limit, and maxit have now been set as members of the MLEmodel class

// this algorithm is optimized specificity for the two parameter case
// variables to hold number of  parameters and number of vertices are vestigial
		int npar=2, k= npar+1;
 // coefficients for reflection, expansion, and contraction
		double ALPHA=1.0, BETA=0.5, GAMMA=2.0;
 // set the sign for minimization of negative LogLikelihood
		int sign= -1;
// declare objects that must exist beyond dist_num if blocks		
		Rcpp::NumericVector outvec(4);	
		arma::colvec funval(3);	
		int vs = 0;	
		arma::mat DF;
		int warn=0;	
			
// separate weibull from lognormal treatment here						
// starting with weibull						
if(dist_num == 1)  {						
						
 // construct the initial simplex						
		arma::colvec transvstart(2);				
		arma::colvec expv1(2);				
		transvstart(0) = log(vstart(0));				
		transvstart(1) = log(vstart(1));				
		arma::mat v(2,3);				
		for(int i=0;i<k;i++)  {				
		v.col(i)=transvstart;				
		}				
		v(0,1)=v(0,1)*1.2;				
		v(1,2)=v(1,2)*1.2;				
						
//		arma::colvec funval(3);				
		for(int i=0;i<k;i++)  {				
//		funval(i)=LogLike(v.col(i), sign, dist_num, tz );				
		expv1(0) = exp(v(0,i));				
		expv1(1) = exp(v(1,i));				
//		funval(i)=LogLike(expv1, sign, dist_num, tz );	
		funval(i)=LogLike(expv1, sign, tz );
		}				
						
						
						
 // assign vertex order variables						
		arma::uvec ndx=sort_index(funval);				
		vs=(int) ndx(0);				
		int vh=(int) ndx(1);				
		int vg=(int) ndx(2);				
 // generate the initial error measure						
		double P2avg=arma::as_scalar(mean(v.row(1)));				
		double error=arma::as_scalar(sum(pow((v.row(1)-P2avg), 2.0)/npar));				
 // initialize the output dataframe for construction as a matrix						
	arma::rowvec DFrow(4);				
		DFrow(0)=v(0,vs);				
		DFrow(1)=v(1,vs);				
		DFrow(2)=funval(vs);				
		DFrow(3)=error;				
		DF=DFrow;				
						
 // initialization of variables used in the loop						
		arma::colvec vm, vr, ve, vc;				
		double fr, fe, fc;				
		int loopcount=1;				
//		int warn=0;				
						
 // This is the main loop for the minimization						
	while(error>limit)  {					
 // calculate the centroid						
		vm=(v.col(vs)+v.col(vh))/2.0;				
 // reflect vg to new vertex vr						
		vr=(vm+ALPHA*(vm-v.col(vg)));				
//		fr=LogLike(vr, sign, dist_num, tz );				
//		fr=sign*tryLL(vr, dist_num);				
		expv1(0) = exp(vr(0));				
		expv1(1) = exp(vr(1));				
//		fr=sign*tryLL(expv1, dist_num);	
		fr=sign*tryLL2(expv1, tz);	
		if(fr==0) {				
			warn=2;			
			break;			
		}				
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
//				fc=LogLike(vc, sign, dist_num, tz );		
				expv1(0) = exp(vc(0));		
				expv1(1) = exp(vc(1));		
//				fc=LogLike(expv1, sign, dist_num, tz );	
				fc=LogLike(expv1, sign, tz );	
 // upon acceptance replace vg with contracton						
				if (fc < funval(vg)) {		
					v.col(vg) = vc;	
					funval(vg)=fc;	
				}else{		
 // upon rare-to-never case shrink the simplex						
					v.col(vh)=v.col(vs)+(v.col(vh)-v.col(vs))/2.0;	
					v.col(vg)=v.col(vs)+(v.col(vg)-v.col(vs))/2.0;	
 // This case results in two function calculations and simply replacing vh and vg points						
//					funval(vh)=LogLike(v.col(vh), sign, dist_num, tz );	
					expv1(0) = exp(v(0,vh));	
					expv1(1) = exp(v(1,vh));	
//					funval(vh)=LogLike(expv1, sign, dist_num, tz );	
					funval(vh)=LogLike(expv1, sign, tz );	
//					funval(vg)=LogLike(v.col(vg), sign, dist_num, tz );	
					expv1(0) = exp(v(0,vg));	
					expv1(1) = exp(v(1,vg));	
//					funval(vg)=LogLike(expv1, sign, dist_num, tz );
					funval(vg)=LogLike(expv1, sign, tz );	
				}		
			}else{			
// now we make an expansion						
				ve=vm+GAMMA*(vr-vm);		
//				fe=LogLike(ve, sign, dist_num, tz );		
				expv1(0) = exp(ve(0));		
				expv1(1) = exp(ve(1));		
//				fe=sign*tryLL(expv1, dist_num);	
				fe=sign*tryLL2(expv1, tz);
				if(fe==0) {		
					warn=3;	
					break;	
				}		
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
		error=arma::as_scalar(sum(pow((v.row(1)-P2avg), 2.0)/npar));				
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
						
						
//		Rcpp::NumericVector outvec(4);				
//		if(dist_num==1) {				
		outvec[0]=exp(v(1,vs));				
		outvec[1]=exp(v(0,vs));				
//		}else{				
//		outvec[0]=v(0,vs);				
//		outvec[1]=v(1,vs);				
//		}				
// end of weibull						
						
}else{						
// start the lognormal treatment						
						
 // construct the initial simplex						
//		arma::colvec transvstart;				
//		arma::colvec expv1;				
//		transvstart(0) = vstart(0);				
//		transvstart(1) = log(vstart(1));				
		arma::mat v(2,3);				
		for(int i=0;i<k;i++)  {				
		v.col(i)=vstart;				
		}				
		v(0,1)=v(0,1)*1.2;				
		v(1,2)=v(1,2)*1.2;				
						
//		arma::colvec funval(3);				
		for(int i=0;i<k;i++)  {				
//		funval(i)=LogLike(v.col(i), sign, dist_num, tz );
		funval(i)=LogLike(v.col(i), sign, tz );
//		expv1(0) = v(0,i);				
//		expv1(1) = exp(v(1,i));				
//		funval(i)=LogLike(expv1, sign, dist_num, tz );				
		}				
						
						
						
 // assign vertex order variables						
		arma::uvec ndx=sort_index(funval);				
		vs=(int) ndx(0);				
		int vh=(int) ndx(1);				
		int vg=(int) ndx(2);				
 // generate the initial error measure						
		double P2avg=arma::as_scalar(mean(v.row(1)));				
		double error=arma::as_scalar(sum(pow((v.row(1)-P2avg), 2.0)/npar));				
 // initialize the output dataframe for construction as a matrix						
		arma::rowvec DFrow(4);				
		DFrow(0)=v(0,vs);				
		DFrow(1)=v(1,vs);				
		DFrow(2)=funval(vs);				
		DFrow(3)=error;				
		DF=DFrow;				
						
 // initialization of variables used in the loop						
		arma::colvec vm, vr, ve, vc;				
		double fr, fe, fc;				
		int loopcount=1;				
//		int warn=0;				
						
 // This is the main loop for the minimization						
	while(error>limit)  {					
 // calculate the centroid						
		vm=(v.col(vs)+v.col(vh))/2.0;				
 // reflect vg to new vertex vr						
		vr=(vm+ALPHA*(vm-v.col(vg)));				
//		fr=LogLike(vr, sign, dist_num, tz );				
//		fr=sign*tryLL(vr, dist_num);
		fr=sign*tryLL2(vr, tz);
//		expv1(0) = vr(0);				
//		expv1(1) = exp(vr(1));				
//		fr=sign*tryLL(expv1, dist_num);				
		if(fr==0) {				
			warn=2;			
			break;			
		}				
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
//				fc=LogLike(vc, sign, dist_num, tz );
				fc=LogLike(vc, sign, tz );
//				expv1(0) = vc(0);		
//				expv1(1) = exp(vc(1));		
//				fc=LogLike(expv1, sign, dist_num, tz );		
 // upon acceptance replace vg with contracton						
				if (fc < funval(vg)) {		
					v.col(vg) = vc;	
					funval(vg)=fc;	
				}else{		
 // upon rare-to-never case shrink the simplex						
					v.col(vh)=v.col(vs)+(v.col(vh)-v.col(vs))/2.0;	
					v.col(vg)=v.col(vs)+(v.col(vg)-v.col(vs))/2.0;	
 // This case results in two function calculations and simply replacing vh and vg points						
//					funval(vh)=LogLike(v.col(vh), sign, dist_num, tz );	
					funval(vh)=LogLike(v.col(vh), sign, tz );	
//					expv1(0) = v(0,vh);	
//					expv1(1) = exp(v(1,vh));	
//					funval(vh)=LogLike(expv1, sign, dist_num, tz );	
//					funval(vg)=LogLike(v.col(vg), sign, dist_num, tz );	
					funval(vg)=LogLike(v.col(vg), sign, tz );
//					expv1(0) = v(0,vg);	
//					expv1(1) = exp(v(1,vg));	
//					funval(vg)=LogLike(expv1, sign, dist_num, tz );	
				}		
			}else{			
// now we make an expansion						
				ve=vm+GAMMA*(vr-vm);		
//				fe=LogLike(ve, sign, dist_num, tz );
				fe=LogLike(ve, sign, tz );
//				expv1(0) = ve(0);		
//				expv1(1) = exp(ve(1));		
//				fe=sign*tryLL(expv1, dist_num);		
				if(fe==0) {		
					warn=3;	
					break;	
				}		
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
		error=arma::as_scalar(sum(pow((v.row(1)-P2avg), 2.0)/npar));				
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
						
						
//		Rcpp::NumericVector outvec(4);				
//		if(dist_num==1) {				
//		outvec[0]=exv(p(1,vs));				
//		outvec[1]=v(0,vs);				
//		}else{				
		outvec[0]=v(0,vs);				
		outvec[1]=v(1,vs);				
//		}				
// end of lognormal treatment						
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

std::vector<double> MLEmodel::genTzvec(double start, double end, double maxtz, int npts)  {		
double spacing;		
std::vector<double> v(npts);		
	if(end == maxtz)  {	
		spacing = (end-start)/npts;  // protecting maxtz
	}else{	
	// or we get a sequence including the end point	
	// spacing will always take the sign of the end point	
		spacing =(end-start)/(npts-1);
	}	
	int n=0;	
	std::generate(v.begin(), v.end(), [ &n, &spacing, &start]() mutable { return start+(n++ * spacing);});	
return v;		
}		


SEXP MLEmodel::MLE3p(SEXP arg3, SEXP arg4, SEXP arg5) {									
	arma::colvec vstart=Rcpp::as<arma::colvec>(arg3);								
	double maxtz=Rcpp::as<double>(arg4);								
	Rcpp::List L(arg5);								
	int num_pts=Rcpp::as<int>(L["num_points"]);								
	double err_t0_limit=Rcpp::as<double>(L["err_t0_limit"]);								
	double err_gof_limit=Rcpp::as<double>(L["err_gof_limit"]);								
									
// Declare objects to retain scope within the loop								
	Rcpp::NumericVector out_vec(4);								
	arma::colvec try_list;								
	bool positive_runout = false;								
	bool negative_runout = false;								
	bool rebound = false;								
	double rebound_value= 0.0;								
									

	
//	std::vector<Rcpp::NumericVector>DF(num_pts);								
									
//set up for first trial with items that will be accessed/updated inside the loop									
	double start=0.0;								
	double end = maxtz;								
	bool t0_found = false;								
	int trial = 1;								
	int max_ind;								
	double err_t0=1.0;								
	double err_gof=1.0;


// Declare the workhorse objects outside the while loop, they are redefined on each pass									
	Rcpp::NumericVector Simplex_out(4);								

	arma::colvec tzvec(num_pts);	
	arma::mat DFmat(num_pts, 4);
	arma::colvec P1(num_pts);
	arma::colvec P2(num_pts);
	arma::colvec gof_vec(num_pts);	
									
while(!t0_found)  {
//for( int i=0; i<2; i++) {	
									
	tzvec = genTzvec(start, end, maxtz, num_pts);


	for( int index=0; index<num_pts; index++) {								
		Simplex_out =  Rcpp::as<Rcpp::NumericVector>(MLE_Simplex(vstart, tzvec[index], 0));							
		gof_vec(index) = Simplex_out[2];
		P1(index) = Simplex_out[0];
		P2(index) = Simplex_out[1];
//Rcout<<"trial "<<trial<<"  index "<<index<<"  vstart[0]  "<<vstart[0]<<"  vstart[1]  "<<vstart[1]<<"   gof_vec(index) "<<gof_vec(index)<<"\n"; 	
//Rcout<<"P1 "<<P1(index)<<"  P2 "<<P2(index)<<"  tz  "<<tzvec[index]<<"  gof  "<<gof_vec(index)<<"\n"; 	
		DFmat = arma::join_rows(P1, P2, tzvec, gof_vec);
		// update vstart here					
		if(dist_num == 1) {				
			vstart[0] = P2(index);			
			vstart[1] = P1(index);			
		}else{				
			vstart[0] = P1(index);			
			vstart[1] = P2(index);	
		}
	}
		
// this code appears to be the problem, can't seem to concatinate the tzvec onto a try_list

		if(trial == 1) {
			try_list = tzvec;
		}else{
			try_list = arma::join_cols(try_list, tzvec);
		}
/*
// diagnostic for try_list		
		for(unsigned int i=0; i<try_list.size(); i++) {
Rcout<<"try "<<i<<" tz "<<try_list(i)<<"\n";
		}
*/

	max_ind = (int) arma::index_max(gof_vec);
//	Rcout<<"next trial "<<trial+1<<"     max_ind  "<<max_ind<<"\n"; 
	
								
// update vstart here					
	if(dist_num == 1) {				
		vstart[0] = P2(max_ind);			
		vstart[1] = P1(max_ind);			
	}else{				
		vstart[0] = P1(max_ind);			
		vstart[1] = P2(max_ind);			
	}				
					
//	Rcout<<"next trial "<<trial+1<<"     vstart[0]  "<<vstart[0]<<"     vstart[1]  "<<vstart[1]<<"\n"; 
				
	
								
									
	if(max_ind != 0)  {								
// get error measures for tz and gof by comparison with max_ind-1, for use in various locations as well as output									
		err_t0 = abs((tzvec[max_ind]-tzvec[max_ind-1])/(tzvec[max_ind]));							
		err_gof = abs((gof_vec[max_ind] - gof_vec[max_ind-1])/(gof_vec[max_ind]));							
	}								
									
	if(max_ind != 0 && max_ind != (num_pts-1)) {								
		if( err_t0 > err_t0_limit) {							
// establish simplex optcontrol for next trial if necessary									
			if(err_gof/num_pts < 1e-5) setSimplexLimit(err_gof/num_pts);						
			start = tzvec[max_ind-1];						
			end = tzvec[max_ind+1];	
//Rcout<<"trial "<<trial<<"     start  "<<start<<"     end "<<end<<"\n";
		
		}else{							
			t0_found = true;						
		}							
	}else{								
		if(trial==1 && max_ind==0)  {							
		// reverse the seek for negative t0							
			start = tzvec[1];		// make sure to cover all early positve values				
			end = (-1) * maxtz;						
		}else{							
			if(tzvec[num_pts-1] > 0) {						
				if(max_ind == num_pts-1) {					
					// get err_t0, then check if less than default limit				
					if(err_t0 < err_t0_limit) {				
						t0_found = true;			
						positive_runout = true;			
					}else{				
						// establish simplex optcontrol for next trial if necessary			
						if(err_gof/num_pts < 1e-5) setSimplexLimit(err_gof/num_pts);			
					// check for rebound case
/*					
						arma::vec shift_gof = join_cols((arma::colvec) gof_vec[0], gof_vec.head(num_pts-1));			
						arma::vec progression = gof_vec-shift_gof;			
						int rebound_ind = (int) arma::index_min(progression); 
*/
std::vector<double> progression;	
progression.push_back(0.0);	
for(int item=1; item<num_pts; item++)  {	
	progression.push_back(gof_vec[item] - gof_vec[item-1]);
}	
int rebound_ind = std::min_element(progression.begin(), progression.end())-progression.begin();	

						if(rebound_ind!=0)  {			
						// if so, next trial is a rework of this trial with end set to rebound point			
							// start is unchanged		
							end = tzvec[rebound_ind];		
							// decrement the trial so it will replace last and flag this unusual event (for further study?).		
							trial = trial-1;		
							rebound = true;		
							rebound_value = tzvec[rebound_ind];		
						}else{			
							start = tzvec[num_pts-1];		
							end = maxtz;		
						}			
					}				
				}else{					
				// max_ind must be 0, the seek is still positive t0					
					start = try_list[(int) try_list.n_rows - (num_pts+2)];
					// go back one increment from tzvec[0]
					//start = tzvec[0]-(tzvec[1]-tzvec[0]);
					end = tzvec[1];				
				}					
			}else{						
			// get err_gof, then check if less than default limit						
				if(err_gof < err_gof_limit) {					
					t0_found = true;				
					negative_runout = true;				
				}else{					
					// establish next mlefit limit if necessary				
					if(err_gof/num_pts < 1e-5) setSimplexLimit(err_gof/num_pts);				
					start = end;				
					end = end*10;				
				}					
			}						
		}							
	}

//	if(trial>9) {t0_found = TRUE;}
//Rcout<<"trial "<<trial<<"     start  "<<start<<"     end "<<end<<"\n";	
	if(!t0_found)  {								
		trial++;							
	}
//out_vec = {DF[max_ind][0], DF[max_ind][1], tzvec[max_ind], gof_vec[max_ind]};	
//out_vec = DFmat.row(max_ind);
//close the main loop finding t0									
}

// this arma::mat has been returning as a SEXP without the Rcpp::wrap			
out_vec = Rcpp::wrap(DFmat.row(max_ind));
									
Rcpp::List Lout = Rcpp::List::create(	
	Rcpp::Named("outvec")=out_vec,
	Rcpp::Named("positive_runout")=Rcpp::wrap(positive_runout),
	Rcpp::Named("negative_runout")=Rcpp::wrap(negative_runout),
	Rcpp::Named("rebound")=Rcpp::wrap(rebound),
	Rcpp::Named("rebound_value")=Rcpp::wrap(rebound_value),
	Rcpp::Named("try_list")=Rcpp::wrap(try_list)
);	

return Lout;	
								
}									


	// Exported Functions

	SEXP MLEloglike(SEXP arg1, SEXP arg3, SEXP arg5, SEXP arg6)  {
		std::unique_ptr<MLEmodel> mymodel(new MLEmodel(arg1));
		arma::colvec par=Rcpp::as<arma::colvec>(arg3);
		//int _dist_num=Rcpp::as<int>(arg4);
		int sign=Rcpp::as<int>(arg5);
		double tz=Rcpp::as<double>(arg6);
		return wrap(mymodel->LogLike(par, sign, tz));
		//return wrap(mymodel.getDistNum());
	}

	SEXP MLEsimplex(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5)  {		
		std::unique_ptr<MLEmodel> mymodel(new MLEmodel(arg1));	
		Rcpp::List L(arg2);	
		mymodel->setSimplexLimit(Rcpp::as<double>(L["limit"]));	
		mymodel->setSimplexMaxit(Rcpp::as<double>(L["maxit"]));	
		arma::colvec vstart=Rcpp::as<arma::colvec>(arg3);	
		double tz=Rcpp::as<double>(arg4);	
		int listout=Rcpp::as<int>(arg5);	
		return mymodel->MLE_Simplex(vstart, tz, listout);	
	}
	
	SEXP callMLE3p(SEXP arg1,  SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5) {		
		std::unique_ptr<MLEmodel> model(new MLEmodel(arg1)); 	
		Rcpp::List L(arg2);	
		model->setSimplexLimit(Rcpp::as<double>(L["limit"]));	
		model->setSimplexMaxit(Rcpp::as<double>(L["maxit"]));			
		return model->MLE3p(arg3, arg4, arg5);
		//return wrap(model->getDistNum());
// Note: listout is handled in R code upon return		
}		
