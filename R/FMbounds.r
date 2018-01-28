FMbounds<-function(x, dist="weibull", CI=.90, unrel=NULL, debias="none", show=FALSE)  {

##  x must be an lrq dataframe such as returned by mleframe
	if(class(x)!="data.frame") {stop("FMbounds takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("FMbounds takes a structured dataframe input, use mleframe")}
	xnames<-names(x)
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {
		 stop("FMbounds takes a structured dataframe input, use mleframe")  }

	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){
		stop("3-parameter distributions not handled by FMbounds")
	}

	if(length(unrel)>0)  {
	dq<-unrel
	}else{
	## these descriptive quantiles match Minitab unchangeable defaults
	dq=c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))
	}

	K<-qnorm(1-(1-CI)/2)
	##  validity checking of arguments will also be performed in mlefit

	fit<-mlefit(x, dist=dist)

## get failure and suspension counts from the input lrq_frame for possible debias calculations
	Qs<-sum(x$qty[x$right<0])
	Qx<-sum(x$qty)-Qs

	if(tolower(dist) %in% c("weibull","weibull2p")){
		shape<-fit[2]
		scale<-fit[1]
## if bias adjustment is made, a modified hessian is determined
## thus this becomes a "modified" Fisher Matrix bound calculation
		if(debias=="rba") shape<- shape*rba(Qx,dist="weibull")
		if(debias=="hrbu") shape<- shape*hrbu(Qx,Qs)

		hessian<-optimHess(c(shape,scale),wblrLoglike, x=x, dist="weibull", sign=-1)
		V<-solve(hessian)
		yp<-log(log(1/(1-dq)))

		xp<-scale*(log(1/(1-dq)))^(1/shape)
		Vt<-V[2,2]/scale^2+yp^2*V[1,1]/shape^4-2*yp*V[1,2]/(shape^2*scale)
		Lb<-log(scale)+yp/shape-K*sqrt(Vt)
		Ub<-log(scale)+yp/shape+K*sqrt(Vt)


		if(show==TRUE)  {
			plot(log(xp),yp, type="l")
			lines(Lb,yp, col="red")
			lines(Ub,yp, col="blue")
		}
	}else{
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p")){
			meanlog<-fit[1]
			sdlog<-fit[2]
## if bias adjustment is made, a modified hessian is determined
## thus this becomes a "modified" Fisher Matrix bound calculation
		if(debias=="rba") sdlog<-sdlog*rba(Qx, dist="lognormal")

			hessian<-optimHess(c(meanlog,sdlog),wblrLoglike, x=x, dist="lognormal", sign=-1)
			V<-solve(hessian)
			yp<-qnorm(dq,0,1)
			Vt<-V[1,1] + yp^2*V[2,2] + 2*yp*V[1,2]


			lnxp<-yp*sdlog+meanlog
			Lb<-lnxp-K*sqrt(Vt)
			Ub<-lnxp+K*sqrt(Vt)
			xp<-exp(lnxp)

			if(show==TRUE)  {
				plot(lnxp,yp, type="l")
				lines(Lb,yp, col="red")
				lines(Ub,yp, col="blue")
			}
		}else{
			stop("distribution not identified in FMbounds")
		}
	}

	outDF<-data.frame(percentile=dq*100, lower=exp(Lb), datum=xp, upper=exp(Ub))
	return(outDF)
}
