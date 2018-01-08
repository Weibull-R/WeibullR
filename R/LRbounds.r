LRbounds<-function(x,  dist="weibull", CL=0.9, dq=(c(1,5,10,20,30,40,50,60,80,90,95,99)/100),  contour=NULL, dof=1, debias=FALSE,applyFF=FALSE, show=FALSE)  {
## check basic parameters of x
	if(class(x)!="data.frame") {stop("mlefit takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("mlefit takes a structured dataframe input, use mleframe")}
	xnames<-names(x)
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {
		 stop("mlefit takes a structured dataframe input, use mleframe")  }
	if(missing(contour))  {
		contour<-MLEcontour(x, dist, CL, DF=dof, debias=debias, applyFF=applyFF)
	}

	if(dist=="weibull")  {
		ypts<-log(qweibull(dq,1,1))
	}else{
## Need lognormal p2y here
		ypts<-qnorm(dq,0,1)
	}



##		i=1
		j=1

##	yval<-c(Blife=ypts)

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
##	for(i in 1:4)  {
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
##	}

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
				P2<- P2*hbru(Qx,Qs)
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




	outDF<-data.frame(percentile=dq*100, lower=exp(outmat[2,]), datum=exp(xvals), upper=exp(outmat[5,]))
##	outDF<-data.frame(ypts=ypts,Lower=outmat[2,],Datum=xvals, Upper=outmat[5,])
outDF
}
