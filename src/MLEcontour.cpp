

#include "WeibullR.h"
#include "MLEmodel.h"
#include "MLEcontour.h"
#include <math.h>

    using namespace Rcpp ;

// class implementation

MLEcontour::MLEcontour( SEXP arg1, arma::colvec arg2, int arg3, double arg4, double arg5, double arg6) {
	model = MLEmodel(arg1);
	par_hat = arg2;
	dist_num = arg3;
	MLLx =  arg4;
	RatioLL = arg5;
	RadLimit = arg6;
}


arma::rowvec MLEcontour::getContourPt( double theta)  {

	double Rincr = 2.5;
//	int CorrLevel = -1;
//	int IncrLevel = 0;
	double r_test = 0;
	arma::colvec innerPt = par_hat;
	double innerR = 0;
	double innerLL = MLLx;
	bool outer_found = false;
	double LL_test = 0.0;
	arma::colvec outerPt(2);
	arma::colvec testPt(2);
	double outerR = 0.0;
	double outerLL = 0.0;

//  main loop seeking outer point
	while(!outer_found)  {
		Rincr = Rincr/5;
		//CorrLevel++;
		r_test = r_test + Rincr;
		if(dist_num==1) {
			testPt(0) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(0));
			testPt(1) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(1));
		}else{
			testPt(0) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(0));
			testPt(1) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(1));
		}

// This is the initial attempt to identify the outerPt, it is primarily successful
// So it is handled first before any incrementing loop
		LL_test = model.tryLL(testPt, dist_num);
		if(LL_test < 0) {
			if(LL_test < RatioLL) {
				outerPt = testPt;
				outerR = r_test;
				outerLL = LL_test;
				outer_found = true;
//			}else{
// innerPt should remain unchanged here
// outerPt has been found - r_test will be set from innerR based on gradient_fraction in follow-on code
//				r_test = innerR;
			}
		}else{
// in the rare case that the initial Rincr resulted in excess parameter
// innerPt is already set at par_hat
// but r_test needs to reset to zero (innerR)
			r_test = innerR;
// Ready to return to main while(!outer_found) loop
// at that point Rincr will be reduced and CorrLevel will be incremented
		}


		while(LL_test > RatioLL && Rincr>RadLimit)  {
			//IncrLevel++;
			r_test = r_test + Rincr;
			if(dist_num==1) {
				testPt(0) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(0));
				testPt(1) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(1));
			}else{
				testPt(0) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(0));
				testPt(1) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(1));
			}


			LL_test = model.tryLL(testPt, dist_num);

			if(LL_test < 0) {
				if(LL_test < RatioLL) {
					outerPt = testPt;
					outerR = r_test;
					outerLL = LL_test;
					outer_found = true;
				}else{
					innerPt = testPt;
					innerR = r_test;
					innerLL = LL_test;
				}
			}else{
	// LL failed to calculate, most likely due to negative parameter
	// possibly due to excess parameter
	// break out of this Rincr level and move to the main loop with last found inner_r
				r_test = innerR;
				break;
			}

	// Here is the opportunity to test for Rincr<RadLimit to establish a close approximation of an outer point find
			if(Rincr<RadLimit)  {
	// the test  point is just as good as an identified contourPt
	// outerPt will be same a innerPt, later to be determined to be contourPt
				outerPt = innerPt;
				outerR = innerR;
				outerLL = innerLL;
				outer_found = true;
			}

	// close the IncrLevel while loop
		}

// close this CorrLevel while loop (it could be the final with Rincr < RadLimit)
	}


// completion of getOuter routine
// start of final seek of contourPt

	double gradient_fraction = 0.0;
	double delta_r = outerR-innerR;
	double delta_LL = innerLL - outerLL;
	//int test_i = 0;
// the unstablePt is a flag, but must be double to return in an arma::colvec
	double unstablePt = 0.0;
	double contourR = 0.0;
	arma::rowvec contourPt(3);
	double Beta = 0.0;
	double Eta = 0.0;
	double Mulog = 0.0;
	double Sdlog = 0.0;

	while(delta_r > RadLimit)  {
		//test_i++;
		gradient_fraction = (innerLL - RatioLL)/ delta_LL;
		if(gradient_fraction < .25) {
			r_test = innerR + .2*delta_r;
		}else{
			if(gradient_fraction < .5) {
				r_test = innerR + .4*delta_r;
			}else{
				if(gradient_fraction > .75) {
					r_test = innerR + .8*delta_r;
				}else{
					r_test = innerR + .6*delta_r;
				}
			}
		}

		if(dist_num==1) {
			testPt(0) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(0));
			testPt(1) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(1));
		}else{
			testPt(0) = (1+ r_test*sin(theta))*arma::as_scalar(par_hat(0));
			testPt(1) = (1+ r_test*cos(theta))*arma::as_scalar(par_hat(1));
		}


		LL_test = model.tryLL(testPt, dist_num);

		if(LL_test < 0) {
			if(LL_test < RatioLL) {
				outerR = r_test;
				outerLL = LL_test;
			}else{
				innerR = r_test;
				innerLL = LL_test;
			}
		}else{
// some unexpected instability has occurred set outer pt as inner pt to end function
// this should never occur, but needed to assure stable run
			unstablePt = 1.0;
			outerR = innerR;
			outerLL = innerLL;
		}

		delta_r = outerR-innerR;
		delta_LL = innerLL - outerLL;
// closure of final contourPt search loop
	}

	if(LL_test > RatioLL) {
		contourR = innerR;
	}else{
		contourR = outerR;
	}

	if(dist_num == 1)  {
		Beta = (1+ contourR*cos(theta))*arma::as_scalar(par_hat(0));
		Eta = (1+ contourR * sin(theta))*arma::as_scalar(par_hat(1));

		contourPt(0) = Eta;
		contourPt(1) = Beta;
		contourPt(2) = unstablePt;
	}else{
		Mulog =  (1+ contourR * sin(theta))*arma::as_scalar(par_hat(0));
		Sdlog = (1+ contourR*cos(theta))*arma::as_scalar(par_hat(1));

		contourPt(0) = Mulog;
		contourPt(1) = Sdlog;
		contourPt(2) = unstablePt;
	}



	return contourPt;

// Close getContourPt method
}



	// Exported Functions

SEXP getContour(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6, SEXP arg7) {
	arma::colvec par=Rcpp::as<arma::colvec>(arg2);
	int dist_num=Rcpp::as<int>(arg3);
	double MLLx = Rcpp::as<double>(arg4);
	double RatioLL = Rcpp::as<double>(arg5);
	double RadLimit = Rcpp::as<double>(arg6);
	int ptDensity=Rcpp::as<int>(arg7);
	MLEcontour mycontour(arg1, par, dist_num, MLLx, RatioLL, RadLimit);

	arma::mat contourpts(ptDensity+1,3);
	const double pi = 3.14159265358979323846;

	double theta = pi;
	for(int row_num=0; row_num<ptDensity+1; row_num++) {
		contourpts.row(row_num) = mycontour.getContourPt(theta);
		theta = theta + 2*pi/( (double) ptDensity);
	}

	return wrap(contourpts);
}

