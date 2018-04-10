





weibayes.mle<-function(x, beta=NULL, eta=NULL, incr=1e-7, listout=FALSE)  {
## check basic parameters of x
	if(class(x)!="data.frame") {stop("mlefit takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("mlefit takes a structured dataframe input, use mleframe")}
	xnames<-names(x)
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {
		 stop("mlefit takes a structured dataframe input, use mleframe")  }

## internal functions
	eta_bounds<-function(x, beta) {

	## find the contour that encompases this beta value
		for(cl in seq(.1,.9, by=.1)) {
			cntr<-MLEcontour(x, CL=cl)
			beta_range<-range(cntr$Beta)
			if(beta_range[1]<beta && beta_range[2]>beta)  break
			if(cl==.9) {
				stop("provided beta is beyond 90% confidence level")
			}
		}
	## get the bounding points around beta for interpolation of bounding eta values
		upper_pos<-which(cntr$Beta>beta)
		lower_pos<-which(cntr$Beta<beta)
		if(length(lower_pos)>length(upper_pos))  {
	## interpolate the eta bounds per upper_pos
			bratio1<-(cntr$Beta[min(upper_pos)]-beta)/(cntr$Beta[min(upper_pos)]-cntr$Beta[min(upper_pos)-1])
			eta1delta<-cntr$Eta[min(upper_pos)]-cntr$Eta[min(upper_pos)-1]
			eta1<-cntr$Eta[min(upper_pos)]-eta1delta *bratio1
			bratio2<-(cntr$Beta[max(upper_pos)]-beta)/(cntr$Beta[max(upper_pos)]-cntr$Beta[max(upper_pos)+1])
			eta2delta<-cntr$Eta[max(upper_pos)]-cntr$Eta[max(upper_pos)+1]
			eta2<-cntr$Eta[max(upper_pos)]-eta2delta *bratio2
		}else{
	## interpolate the eta bounds per lower_pos
			bratio1<-(cntr$Beta[min(lower_pos)]-beta)/(cntr$Beta[min(lower_pos)]-cntr$Beta[min(lower_pos)-1])
			eta1delta<-cntr$Eta[min(lower_pos)]-cntr$Eta[min(lower_pos)-1]
			eta1<-cntr$Eta[min(lower_pos)]-eta1delta *bratio1
			bratio2<-(cntr$Beta[max(lower_pos)]-beta)/(cntr$Beta[max(lower_pos)]-cntr$Beta[max(lower_pos)+1])
			eta2delta<-cntr$Eta[max(lower_pos)]-cntr$Eta[max(lower_pos)+1]
			eta2<-cntr$Eta[max(lower_pos)]-eta2delta *bratio2
		}

		return(c(eta1, eta2))
	}

	beta_bounds<-function(x, eta) {

	## find the contour that encompases this beta value
		for(cl in seq(.1,.9, by=.1)) {
			cntr<-MLEcontour(x, CL=cl)
			eta_range<-range(cntr$Eta)
			if(eta_range[1]<eta && eta_range[2]>eta)  break
			if(cl==.9) {
				stop("provided eta is beyond 90% confidence level")
			}
		}

	## get the bounding points around beta for interpolation of bounding eta values
		upper_pos<-which(cntr$Eta>eta)
		lower_pos<-which(cntr$Eta<eta)
		if(length(lower_pos)>length(upper_pos))  {
	## interpolate the beta bounds per upper_pos
			eratio1<-(cntr$Eta[min(upper_pos)]-eta)/(cntr$Eta[min(upper_pos)]-cntr$Eta[min(upper_pos)-1])
			beta1delta<-cntr$Eta[min(upper_pos)]-cntr$Eta[min(upper_pos)-1]
			beta1<-cntr$Beta[min(upper_pos)]-beta1delta *eratio1
			eratio2<-(cntr$Eta[max(upper_pos)]-eta)/(cntr$Eta[max(upper_pos)]-cntr$Eta[max(upper_pos)+1])
			beta2delta<-cntr$Beta[max(upper_pos)]-cntr$Beta[max(upper_pos)+1]
			beta2<-cntr$Beta[max(upper_pos)]-beta2delta *eratio2
		}else{
	## interpolate the eta bounds per lower_pos
			eratio1<-(cntr$Eta[min(lower_pos)]-eta)/(cntr$Eta[min(lower_pos)]-cntr$Eta[min(lower_pos)-1])
			beta1delta<-cntr$Beta[min(lower_pos)]-cntr$Beta[min(lower_pos)-1]
			beta1<-cntr$Beta[min(lower_pos)]-beta1delta *eratio1
			eratio2<-(cntr$Eta[max(lower_pos)]-eta)/(cntr$Eta[max(lower_pos)]-cntr$Eta[max(lower_pos)+1])
			beta2delta<-cntr$Beta[max(lower_pos)]-cntr$Beta[max(lower_pos)+1]
			beta2<-cntr$Beta[max(lower_pos)]-beta2delta *eratio2
		}

		return(c(beta1, beta2))
	}

## Main Function in two parts, known beta, or known eta
## first known beta (most likely to be called)
	if(length(beta)>0)  {
		ebounds<-eta_bounds(x,beta)
		lb<-ebounds[1]
		ub<-ebounds[2]
		de<-(ub-lb)*incr
		loop<-1
		lbgrad<- (wblrLoglike(c(beta,lb),x)-wblrLoglike(c(beta,lb-de),x))/de
		ubgrad<- (wblrLoglike(c(beta,ub),x)-wblrLoglike(c(beta,ub-de),x))/de
		outDF<-data.frame(lb,ub,lbgrad,ubgrad)

		while(abs( ubgrad)> 1e-10 && lbgrad > 1e-10 && loop < 100 ) {
			# simply a binary search now
			test<- ub-(ub-lb)*0.5
			testgrad<- (wblrLoglike(c(beta,test),x)-wblrLoglike(c(beta,test-de),x))/de
			if(testgrad>0) lb<-test else ub<-test

			if(testgrad==lbgrad || testgrad==ubgrad) break

			lbgrad<- (wblrLoglike(c(beta,lb),x)-wblrLoglike(c(beta,lb-de),x))/de
			ubgrad<- (wblrLoglike(c(beta,ub),x)-wblrLoglike(c(beta,ub-de),x))/de
			outDF<-rbind(outDF, data.frame(lb,ub,lbgrad, ubgrad))

			loop<-loop+1
		}


	}else{
## now the known eta
		bbounds<-beta_bounds(x,eta)

		lb<-bbounds[1]
		ub<-bbounds[2]
		de<-(ub-lb)*incr
		loop<-1
		lbgrad<- (wblrLoglike(c(lb, eta),x)-wblrLoglike(c(lb-de, eta),x))/de
		ubgrad<- (wblrLoglike(c(ub, eta),x)-wblrLoglike(c(ub-de, eta),x))/de
		outDF<-data.frame(lb,ub,lbgrad,ubgrad)

		while(abs( ubgrad)> 1e-10 && lbgrad > 1e-10 && loop < 100 ) {

			test<- ub-(ub-lb)*.5
			testgrad<- (wblrLoglike(c(test, eta),x)-wblrLoglike(c(test-de, eta),x))/de
			if(testgrad>0) lb<-test else ub<-test

			if(testgrad==lbgrad || testgrad==ubgrad) break

			lbgrad<- (wblrLoglike(c(lb, eta),x)-wblrLoglike(c(lb-de, eta),x))/de
			ubgrad<- (wblrLoglike(c(ub, eta),x)-wblrLoglike(c(ub-de, eta),x))/de
			outDF<-rbind(outDF, data.frame(lb,ub,lbgrad, ubgrad))

			loop<-loop+1
		}
	}

		if(lbgrad<abs(ubgrad)) outval<-lb else outval<-ub

		if(listout==FALSE) {
			return(outval)
		}else{
			return(list(outval, outDF))
		}

}
