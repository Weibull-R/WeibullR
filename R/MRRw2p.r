## MRRw2p.r file
 ## A lightweight function using defaults to generate a median rank regression
 ## 
 ## Author: Jacob T. Ormerod
 ##   (c)2014-2017 OpenReliability.org
##

MRRw2p<-function(x, s=NULL, bounds=FALSE, CI=0.90, show=FALSE)  {			
	xName<-paste(deparse(substitute(x),500),collapse = "\n")		
	thisTitle=paste(xName,"Weibull fit")		
	gotppp<-getPPP(x, s,)		
	probability<-as.vector(gotppp[,2])		
	x<-as.vector(gotppp[,1])		
	fit<-lslr(gotppp)		
	fit<-c(fit,LL=LLw(x,s,fit[1],fit[2]))		
			
	if(bounds==TRUE)  {		
		## using  Jurgen Symynck's method for weibull	
		## the same complete failure adjustment as used for lognormal is shown here	
		## but it is believed to always be 1.0 for Weibull fits	
		## completePPP<-getPPP(x)	
		## qwPPP<-data.frame(data=qweibull(completePPP[,2],1,1), ppp=completePPP[,2])	
		P1<-lslr(getPPP(qweibull(probability,1,1)))[1]	
		P2<-lslr(getPPP(qweibull(probability,1,1)))[2]	
		## the adjustment divisor term here is always expected to be 1 for Weibull	
		## P2<-lslr(getPPP(qweibull(probability,1,1)))[2]/lslr(qwPPP)[2]	
			
		## descriptive quantiles for comparison with SuperSMITH (limit of 15 values)	
		dq<-c(.01, .02, .05, .10, .15, .20, .30, .40, .50,  .60, .70, .80, .90, .95, .99)	
		pivotals<-pivotal.rr(gotppp,R2=0, CI=CI,unrel=dq,P1=P1,P2=P2)	
		## use the slope of the median pivotals to get correction to 1.0	
		median_slope<-(log(log(1/(1-dq[15])))-log(log(1/(1-dq[1]))))/(pivotals[15,2]-pivotals[1,2])	
		median_intercept<-pivotals[10,2]-log(log(1/(1-dq[10])))/median_slope	
			
		adj_piv<-(pivotals-median_intercept)*median_slope	
		## intermediate readings of the slope were generated during development, now not used	
		## median_slope2<-(log(log(1/(1-dq[15])))-log(log(1/(1-dq[9]))))/(adj_piv[15,2]-adj_piv[9,2])	
		## median_intercept2<-adj_piv[10,2]-log(log(1/(1-dq[10])))/median_slope	
			
		## interpret the pivotals for the log plot	
		plot_piv<-(adj_piv)/fit[2]+log(fit[1])	
		## again intermediate readings of the slope and intercept were generated only for development	
		## median_slope3<-(log(log(1/(1-dq[15])))-log(log(1/(1-dq[9]))))/(plot_piv[15,2]-plot_piv[9,2])	
		## median_intercept3<-plot_piv[10,2]-log(log(1/(1-dq[10])))/median_slope3	
			
		## prepare the pivotals for plotted output 	
		LB<-as.vector(plot_piv[,1])	
		DATUM<-as.vector(plot_piv[,2])	
		HB<-as.vector(plot_piv[,3])	
			
		## prepare the pivotals for print output	
		print_piv<-exp(plot_piv)	
		outDF<-data.frame(unrel=dq,	
			LB=as.vector(print_piv[,1]),
			DATUM=as.vector(print_piv[,2]),
			HB=as.vector(print_piv[,3]))
	}		
	     if(show==TRUE)   {		
		 plot(log(x),log(log(1/(1-probability))),pch=19,col="red",	
		main=thisTitle)	
		Xintercept<--log(fit[1])*fit[2]	
		abline(Xintercept,fit[2],col="blue")	
		if(bounds==TRUE)  {	
			lines(LB,log(log(1/(1-dq))),col="blue")
			lines(HB,log(log(1/(1-dq))),col="blue")
			lines(DATUM,log(log(1/(1-dq))), col="magenta")
		}	
	     }		
			
	if(bounds==TRUE)  {		
		return(list(fit,outDF))	
	}else{		
		return(fit)	
	}		
}			
