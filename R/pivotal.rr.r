## pivotal.rr
## This function name has been refactored to avoid confusion with pivotalMC in abremPivotals


pivotal.rr<-function(x, event=NULL, dist="weibull", reg_method="XonY", R2, CI, unrel=NULL, P1=1.0, P2=1.0, S=10^4, seed=1234, ProgRpt=FALSE)  {		
				
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
	dq<-unrel
	}else{
	## these descriptive quantiles match Minitab unchangeable defaults
	dq=c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))
	}

		
		if(dist!="weibull" && P1==1.0) message("lognormal or gumbel sampled with P1=1.0")
		
	S = as.integer(S/10)*10	
	if (S < 10^3) {
## return the full vector or matrix output for special small sampled cases
	    if(R2>0) R2=1.0
		if(CI>0) CI=1.0
	   }	
	if(S>4*10^9)   {
		stop("Samples beyond MAX_INT")
	}
				
	casenum<-0			
	if(reg_method=="YonX") casenum=casenum+1						
	if(dist=="lnorm")casenum=casenum+2			
	if(dist=="gumbel") casenum=casenum+4			
				
##	result<-.Call("pivotalMC", x$ppp, event, c(R2,CI,P1,P2), S, seed, dq, ProgRpt, casenum , package="WeibullR")
	result<-.Call(pivotalMC, x$ppp, event, c(R2,CI,P1,P2), S, seed, dq, ProgRpt, casenum)

return(result)				
}				
