# weibayes.r
#
# WeiBayes is a term for a fitting method applicable to data sets with suspensions but very few failures - even one, that
# certainly would defy linear regression methods and present a challenge even for mle fitting.
# Although not a true Bayesian solution the name is given rise due to the input of a "known" (or "prior")
# value for one of the Weibull parameters. Most often, the Beta.
#
# Reliasoft takes issue with the WeiBayes terminology, preferring to reference it as a 1-parameter Weibull method.
# Since one parameter is given as a "known" the resulting fit is a determinant solution. Whereas a true Bayesian solution
# would represent the "prior" as a probability distribution, such that resulting fit must be presented as an uncertainty.
#
# In reality when there is scant failure information uncertainty abounds. Practicioners have confidence using the
# 1-parameter methods, when the character of the item (system, component, or material) in question can be related to previous
# experience. In particular past experience is expectated to relate to the beta parameter.
#
# The first method is very simplistic mathematically, but uncanny in its ability to identify an optimal likelihood.
# This method function is given the simple name weibayes. This method requires a beta value argument.

# The second method uses a true optimization of  the likelihood based on input of either beta or eta.
# This method function is termed weibayes.mle requiring either a beta value or an eta value as input.
#
# Finally a true Bayesian-Weibull will have to express the output as an uncertainty, expected to be similar to confidence bounding.
# A little more research will be made before this implementation is made.

weibayes<-function(x, s=NULL, beta) {
## this is to be an internal function for now
	expand_qty<-function(x,q) {
		outvec<-NULL
		for(line in 1:length(x)) {
			outvec<-c(outvec, rep(x[line], q[line]))
		}
		outvec
	}

## could x be a dataframe with time and event columns??
	suspensions<-NULL
	if(is.vector(x))  {
		if(anyNA(x))  {
		stop("NA in failure data")
		}
		if(any(x<=0))  {
		stop("non-positive values in failure/occurrence data")
		}


## I'm not convinced this needs to be sorted here, but it doesn't hurt
		fail_vec<-sort(x)



		if(length(s)>0)  {
		if(anyNA(s))  {
		stop("NA  in suspension data")
		}
		if(any(s<=0))  {
		stop("non-positive values in suspension data")
		}
		susp_vec<-sort(s)

		}
## end pure vector argument processing
	times<-c(fail_vec, susp_vec)
	nfail<-length(fail_vec)

	}else{
## here a time-event dataframe can be evaluated, if provided as x
## This is the support for a time-event dataframe
		if (class(x) == "data.frame") {

		## this test is drawn from Abrem.R
			if(is.null(x$time) || is.null(x$event)){
				stop(': Argument \"x\" is missing $time and/or ","$event columns...')
			}

## verify positive time values
			if (anyNA(x$time)) {
				stop("NA in failure or suspension data")
			}
			if (any(x$time<= 0)) {
				stop("non-positive values in failure or suspension data")
			}
## verify 1's and 0's only in event
## using Jurgen's validation code
			ev_info <- levels(factor(x$event))
			if(identical(ev_info,c("0","1")) || identical(ev_info,"1")){
			# okay x$event is indeed holding event indicators
			}else{
			stop("event column not '1' or '0' ")
			}

			if(length(s)>0)  {
			warning("argument 's' ignored when time-event dataframe provided")
			}

			if(!is.null(x$qty)) {
## now that we know there is a qty field we know that the fail and event vectors need to be expanded
## But let's be sure the qty field is all integer, else future havoc could ensue
				if(any(!is.integer(x$qty))) x$qty<-ceiling(x$qty)
				times<-expand_qty(x$time, x$qty)
				events<-expand_qty(x$event, x$qty)

			}else{
				times<-x$time
				events<-x$event
			}
			nfail<-sum(events)

		}  ## close dataframe processing
	}  ## close input argument processing

## after all that data input manipulation
## here is the simple 1-parameter, weibayes method
	 t_eta<-(times^beta)/nfail
	out_val<-sum(t_eta)^(1/beta)

	out_val
}

