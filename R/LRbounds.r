LRbounds<-function(x,  dist="weibull", CL=0.9, unrel=NULL,  contour=NULL, dof=1, ptDensity=120, debias="none", show=FALSE)  {

## check basic parameters of x
	if(class(x)!="data.frame") {stop("mlefit takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("mlefit takes a structured dataframe input, use mleframe")}
	xnames<-names(x)
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {
		 stop("mlefit takes a structured dataframe input, use mleframe")  }
		 
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
		 
		 
if(npar==2) {	
	ptDensity=120
	if(!is.null(control$ptDensity))  ptDensity<-control$ptDensity
		 
	if(missing(contour))  {
		contour<-MLEcontour(x, dist, CL, dof=dof, ptDensity=ptDensity, debias=debias)
	}

	if(dist=="weibull")  {
		ypts<-log(qweibull(dp,1,1))
	}else{
## Need lognormal p2y here
		ypts<-qnorm(dp,0,1)
	}
		
	j=1
## P1 was Eta
	P1<-contour[j,1]
## P2 was Beta
	P2<-contour[j,2]

	xvals=NULL
	for(k in 1:length(ypts) ) {
		if(dist=="weibull") {
			xval<-ypts[k]/P2+log(P1)
		}else{
## Need lognormal p2y here
			xval<-ypts[k]*P2+P1
		}

		xvals<-c(xvals,xval)
		names(xvals)<-NULL
	}

	outmat<-rbind(ypts,xlo=xvals, P1=rep(P1,length(ypts)),P2=rep(P2,length(ypts)),
		xhi=xvals, P1=rep(P1,length(ypts)),P2=rep(P2,length(ypts)))

			clen=length(contour[,1])
		for(j in 1:clen)  {
			P1<-contour[j,1]
			P2<-contour[j,2]
			xvals=NULL
			for(k in 1:length(ypts) ) {
				if(dist=="weibull") {
					xval<-ypts[k]/P2+log(P1)
				}else{
			## Need lognormal p2y here
					xval<-ypts[k]*P2+P1
				}
				if(xval<outmat[2,k])  {
					outmat[2,k]=xval
					outmat[3,k]=P1
					outmat[4,k]=P2
				}
				if(xval>outmat[5,k])  {
					outmat[5,k]=xval
					outmat[6,k]=P1
					outmat[7,k]=P2
				}
			}
		}

## calculate the Datum vector

	MLEfit<-mlefit(x, dist=dist)
	P1<-MLEfit[1]
	P2<-MLEfit[2]
## adjust P2 according to debias
## needed quantity of type information is provided as "data_types" attribute to MLEfit object

	dt<-unname(attributes(MLEfit)$data_types)
	Qx<-dt[1]-dt[2]
	Qs<-dt[2]
	if(!is.null(debias))  {
		if(debias=="rba")  {
			P2<- P2*rba(Qx,dist=dist)
		}else{
			if(debias=="hrbu" && dist=="weibull")  {
				P2<- P2*hrbu(Qx,Qs)
			}else{
				if(debias=="hrbu" && dist=="lnorm") P2<- P2*rba(Qx,dist=dist)
			}
		}
	}



	xvals=NULL
	for(k in 1:length(ypts) ) {

		if(dist=="weibull" )  {
			xval<-ypts[k]/P2+log(P1)
		}else{
			xval<-ypts[k]*P2+P1
		}
		xvals<-c(xvals,xval)
		names(xvals)<-NULL
	}



	if(show==TRUE)  {
		plot(xvals,ypts,type="l")
		lines(outmat[2,],outmat[1,],col="red")
		lines(outmat[5,],outmat[1,],col="blue")
	}




	outDF<-data.frame(percentile=dp*100, lower=exp(outmat[2,]), datum=exp(xvals), upper=exp(outmat[5,]))
	outList<-list(bounds=outDF,contour=contour)
## close 2p handling here.
}
	
	

outList
}

