## pivotal.rr
## This function name has been refactored to avoid confusion with pivotalMC in abremPivotals


pivotal.rr<-function(x, event=NULL, dist="weibull", npar=2, reg_method="XonY", R2, CI, unrel=NULL, P1=1.0, P2=1.0, S=10^4, seed=1234, ProgRpt=FALSE)  {		
			
	if(is.vector(x))  {			
		stop("use MRR functions for casual fitting, or pre-process with getPPP")		
	}else{			
		if(names(x)[1]=="time"&&names(x)[2]=="ppp")  {		
		## will handle the output from getPPP				
		}else{	
			if(length(x$ppp<3))  {
				stop("insufficient failure points")
			}else{
				stop("input format not recognized")	
			}
		}			
	}
		if(missing(event)){
		event<-c(rep(1,length(x[,1])))
	}else{
	## validate the event vector
		zeros<-length(event[sapply(event, function(x) x==0)])
		if(length(x[,1])!=(length(event)-zeros)) {
			stop("event vector has wrong length")
		}
	}
	
	if (R2 < 0|| R2>1) stop("Invalid R-squared value")
	if (CI < 0|| CI>1) stop("Invalid Confidence Interval")
		
	if(length(unrel)>0)  {
	dp<-unrel
	}else{
	## these descriptive percentiles match Minitab unchangeable defaults
	dp=c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))
	}

## It turns out that this code is general to all fitting methods:			
	if(tolower(dist) %in% c("weibull","weibull2p","weibull3p")){		
		dist_num <- 0	
	}else{		
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p", "lognormal3p")){	
			dist_num <-1
		}else{	
			if(!dist=="gumbel") {
		## Note: Original pivotalMC contained experimental support for "gumbel"	
			stop(paste0("dist argument ", dist, "is not recognized for distribution fitting"))
			}
			if(dist=="gumbel") {
			stop("experimental support for gumbel is currently suspended")
			}
		}	
	}		
	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){		
		npar<-3	
	}		
	if(!(npar %in% c(2,3))) stop("unsupported npar value")		
			
## trap mistakes that can be made with parameterization			
	if(P2 == 0) stop("neither weibull beta, nor lognormal sigmalog should be 0")		
	if(dist_num == 0 && !P1>0 ) stop("non-positive weibull eta")		
	if(dist_num == 1 && !P2>0 ) stop("non-positive lognormal sigmalog")		

	S = as.integer(S/10)*10	
	if (S < 10^3) {
## return the full vector or matrix output for special small sampled cases
	    if(R2>0) R2=1.0
		if(CI>0) CI=1.0
	   }	
	if(S>4*10^9)   {
		stop("Samples beyond MAX_INT")
	}
	
	reg_order=0
	if(reg_method=="YonX") reg_order=1
	
## a convergence limit is fixed here for 3rd parameter  convergence
## no longer an argument for the R function, but still an argument to C++ functions
	limit<-1e-5	
		
callargs<-list(ppp=x$ppp, event=event, R2=R2, CI=CI, P1=P1, P2=P2, S=S, seed=seed, dp=dp, 
				reg_order=reg_order, dist_num=dist_num, npar=npar, limit=limit,ProgRpt=ProgRpt)
				
##	result<-.Call("pivotalMC", x$ppp, event, c(R2,CI,P1,P2), S, seed, dq, ProgRpt, casenum , package="WeibullR")
##	result<-.Call(pivotalMC, x$ppp, event, c(R2,CI,P1,P2), S, seed, dq, ProgRpt, casenum)

	result<-.Call(pivotalMC, callargs)

return(result)				
}				
