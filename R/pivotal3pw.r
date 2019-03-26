##File pivotal3pw.r is placed in the WeibullR package temporarily during testing
## this content is expected to be incorporated into pivotal.rr.r in the future.


pivotal3pw<-function(x, s=NULL, CI=0.9, unrel=NULL, S=1000, listout=FALSE, show=FALSE)  {			
## x and s arguments will be interpreted in this getPPP call	
		da<-getPPP(x,s)
		fit3p<-lslr(da, npar=3)	
		## now that the fit return is a vector must use numbered items
		Beta_opt<-fit3p[2]	
		Eta_opt<-fit3p[1]	
		t0_opt<-fit3p[3]	

		## none of this applies at this time
#		if( !is.null(attr(fit3p,"unstable"))) stop("unstable 3rd parameter fit")	
			
#		if(!is.null(attr(fit3p, "rebound_pt")))  {	
#			maxtz<-attr(fit3p,"rebound_pt")
#		}else{	
#			maxtz<-min(x)
#		}	
## actually, I am not sure of this:			
		if(t0_opt==0) stop("t0 = 0, nothing to do")
		
	if(length(unrel)>0)  {
		dq<-unrel
	}else{				
## descriptive quantiles for comparison with SuperSMITH (limit of 15 values)
		dq<-c(.01, .02, .05, .10, .15, .20, .30, .40, .50,  .60, .70, .80, .90, .95, .99)
	}

	#Datum<-qweibull(dq,Beta_opt, Eta_opt) + t0_opt
	set.seed(1234)
# passing the dataframe passes the ppp column	
	sample_data<-da	
	boot.mat<-NULL	
	## generation of  parametric bootstrap	
	for(set in 1:S) {	
		repeat{
			sample_data$time<-sort(rweibull(nrow(da),Beta_opt,Eta_opt)+t0_opt)
			if(!any(sample_data<0)) break
		}
		fit<-lslr(sample_data, npar=3)
		xvals<- qweibull(dq,fit[2], fit[1])+fit[3]
		boot.mat<-rbind(boot.mat, xvals)
	}	

	boot.mat<-apply(boot.mat,2,sort)
	lo_row<-ceiling(S*(1-CI)/2)
	up_row<-floor(S*(1-(1-CI)/2))
	Lower<-boot.mat[lo_row,]
	Lower<-sapply(Lower,function(X) ifelse(X<=0,NA,X))
	Upper<-boot.mat[up_row,]
	Median<-boot.mat[floor(S/2),]

	if(show) {
	if(!exists("p2y")) {
	p2y <- function(p,log="x"){	
		# This is the inverse Cumulative Distribution function
		# used to transform a probability value to the
		# y-axis of the plot canvas.
		# Use of this transformation permits distributions
		# to appear as curves on unrelated canvas
		if(log =="x")ret <- log(qweibull(p,1,1))
		if(log =="xy") ret <- qlnorm(p,0,1)
		ret
	}	
	}
	
	obj1<-wblr.conf(wblr.fit(wblr(x-t0_opt,s-t0_opt, col="red"), col="grey"),lty=2, col="orange")
	plot(obj1, xlab="time - t0", main="Modified Data Plot")	
		
	lines(Median-t0_opt,p2y(dq), col="black")	
	lines(Lower-t0_opt,p2y(dq), col="blue")	
	lines(Upper-t0_opt,p2y(dq), col="blue")				

	## start a new graphics device
	#x11(xpos=-350, ypos=100)
	dev.new(xpos=-350, ypos=100)
	obj2<-wblr.fit(wblr(x,s, col="red"), npar=3, lty=0)	
	plot(obj2)	
		
	lines(Median,p2y(dq), col="black")	
	lines(Lower,p2y(dq), col="blue")	
	lines(Upper,p2y(dq), col="blue")
	
## The easiest way to get the bound lines for such a plot is to use the wblr object. 
## However placement of the conversion code should in the base technical code, not the wblr.conf function. 
## The return from the R function pivotal.rr should be the untransformed values, 
## without regard to graphic presentation.

	bounds2p<-obj1$fit[[1]]$conf[[1]]$bounds
	lines(bounds2p$Lower+t0_opt,p2y(bounds2p$unrel), lty=2, lwd=2, col="orange")
	lines(bounds2p$Upper+t0_opt,p2y(bounds2p$unrel), lty=2, lwd=2, col="orange")
		
	}


	bounds<-data.frame(unreliability=dq, Lower=Lower, Median=Median, Upper=Upper)

## listout could include the boot.mat as a return item
bounds
}			
