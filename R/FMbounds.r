FMbounds<-function(x, dist="weibull", CI=.90, unrel=NULL, debias="none", show=FALSE)  {					
					
##  x must be an lrq dataframe such as returned by mleframe					
	#if(class(x)!="data.frame") stop("FMbounds takes a structured dataframe input, use mleframe")	
	if(!is(x, "data.frame")) stop("FMbounds takes a structured dataframe input, use mleframe")	
	if(ncol(x)!=3)  {stop("FMbounds takes a structured dataframe input, use mleframe")}				
	xnames<-names(x)				
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {				
		 stop("FMbounds takes a structured dataframe input, use mleframe")  }			
					
	npar=2				
	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){				
		npar=3			
	}				
					
	if(tolower(dist) %in% c("weibull","weibull2p", "weibull3p")){				
		dist <- "weibull"			
	}else{				
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p","lognormal3p")){			
			dist <- "lognormal"		
		}else{			
			stop("distribution not identified in FMbounds")		
		}			
	}				
					
	if(length(unrel)>0)  {				
	dp<-unrel				
	}else{				
	## these descriptive percentiles match Minitab unchangeable defaults				
	dp=c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))				
	}				
					
	K<-qnorm(1-(1-CI)/2)				
	##  validity checking of arguments will also be performed in mlefit				
					
					
## get failure and suspension counts from the input lrq_frame for possible debias calculations					
	Qs<-sum(x$qty[x$right<0])				
	Qx<-sum(x$qty)-Qs				
					
	fit<-mlefit(x, dist, npar)				
					
	if(dist == "weibull") {				
## if bias adjustment is made, a modified hessian is determined					
## thus this becomes a "modified" Fisher Matrix bound calculation					
		if(debias=="rba") fit[2]<- fit[2]*rba(Qx,dist="weibull")			
		if(debias=="hrbu") fit[2]<- fit[2]*hrbu(Qx,Qs)			
		mu<-log(fit[1])             # log(Eta)			
		sigma<-1/fit[2]             # 1/Beta			
		q <- qweibull(dp, fit[2], fit[1])			
		z<- log(log(1/(1-dp)))			
	}else{				
## if bias adjustment is made, a modified hessian is determined					
## thus this becomes a "modified" Fisher Matrix bound calculation					
		if(debias=="rba") fit[2]<-fit[2]*rba(Qx, dist="lognormal")			
		mu<-fit[1]                  # Mulog			
		sigma<-fit[2]             #Ssigmalog			
		q <- qlnorm(dp, fit[1], fit[2])			
		z <- qnorm(dp, 0, 1)			
	}				
					
	if(npar==3) q=q+fit[3]				
					
	theta<- c(mu, sigma)				
	names(theta)<-c("mu", "sigma")				
	if(npar == 3) theta<- c(theta, gamma = fit[3])				
					
	hessian<-stats::optimHess(theta, LL, dist=dist, data=x)				
	varcov <- solve(hessian)				
					
	Q<-q				
	if(npar == 3) Q<-q-theta[3]				
	  dq.dmu <- Q   #with respect to mu				
	  dq.dsigma <- z*Q     #with respect to sigma				
	  dq.dgamma <- 1    #with respect to gamma				
					
	std.err<- NULL				
	for(sp in 1:length(dp)) {				
		dq.dtheta<-c(dq.dmu[sp], dq.dsigma[sp])			
		if(npar == 3) dq.dtheta<-c(dq.dtheta, dq.dgamma)			
		  std.err <- c(std.err, sqrt(t(dq.dtheta)%*%varcov%*%dq.dtheta))			
	}				
					
	  #compute confidence interval				
	w<-exp((K*std.err) /q)				
	lcl<-q/ w				
	ucl<- q * w				
					
	if(show == TRUE) {				
		plot(log(q),z, type="l")			
		lines(log(lcl),z, col="red")			
		lines(log(ucl),z, col="blue")			
	}				
					
	outDF<-data.frame(percentile=dp*100, lower=lcl, datum=q, upper=ucl)				
	return(outDF)				
}					

LL<-function(theta, dist, data) {		
	if(dist == "weibull") {	
		P2 <- exp(theta[1])
		P1 <- 1/theta[2]
	}	
	if(dist == "lognormal") {	
		P1 <- theta[1]
		P2 <- theta[2]
	}	
	tz<-theta[3]	
	if(is.na(tz)) tz<-0	
	return(wblrLoglike(c(P1,P2), x=data, dist=dist, sign= -1, tz=tz))	
}		
