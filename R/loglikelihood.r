 # loglikelihood.r file
 #
 #  Functions returning the log likelihood for weibull and log-normal fitted regressions with special attention to removal
 #  of suspension data that may have become negative by the explicit subtraction of the t0, threashold, value
 #  from original data.  The purpose of the log likelihood in a median rank regression environment is to specifically
 #  use the likelihood ratio test to validate the improvement of fit provided by third parameter optimization.
 #  A valid improvement of fit is judged to require an LRT-P value greater than 0.5.
 #
 #  Functions LLw and LLln only handle datasets with failure and suspension data.
 #  For data including intervals wblrLikelihoodmust be used.
 #
 # Author: Jacob T. Ormerod
 #   (c)2014-2021 OpenReliability.org


LLw<-function(x,s=NULL,Eta,Beta)  {
	suscomp<-0
	failcomp<-sum(dweibull(x,Beta,Eta,log=TRUE))
	if(length(s)>0)  {
		if(any(s<=0))  {
			s2<-NULL
			for(i in 1:length(s))  {
				if(s[i]>0) {s2<-c(s2,s[i])}
			}
			s<-s2
			if(length(s)>0)  {
				suscomp<-sum(pweibull(s,Beta,Eta,lower.tail=FALSE,log.p=TRUE))
			}
		}else{
			suscomp<-sum(pweibull(s,Beta,Eta,lower.tail=FALSE,log.p=TRUE))
		}
	}
	value<-failcomp+suscomp
	value
}

LLln<-function(x,s=NULL,Mulog,Sigmalog)  {
	suscomp<-0
	failcomp<-sum(dlnorm(x,Mulog,Sigmalog,log=TRUE))
	if(length(s)>0)  {
		if(any(s<=0))  {
			s2<-NULL
			for(i in 1:length(s))  {
				if(s[i]>0) {s2<-c(s2,s[i])}
			}
			s<-s2
			if(length(s)>0)  {
				suscomp<-sum(plnorm(s,Mulog,Sigmalog,lower.tail=FALSE,log.p=TRUE))
			}
		}else{
			suscomp<-sum(plnorm(s,Mulog,Sigmalog,lower.tail=FALSE,log.p=TRUE))
		}
	}
	value<-failcomp+suscomp
	value
}

## This is the likelihood ratio test.  When used to establish acceptance of 3p model over its 2p counterpart on the same data
## The 2p model is taken as alternate (H1), 3p as null (Ho).  Degrees of freedom (df) is 1 for this constrained model test.
## Acceptance of the 3p model should require an LRT-P greater than 0.50.  ("If P is low Ho must go.")
## Comparison between Weibull and lognormal for failure data would place Weibull as alternate and lognormal as null.
## Degrees of freedom would be 2.
##LRT_P<-function(LLalt,LLnull,df)  {
##	lrtp<-pchisq(-2*(LLalt-LLnull),df)
##}
## Decided not to make this an exported function as it is easy enough to implement when needed.