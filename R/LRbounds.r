LRbounds<-function(x,  dist="weibull", CL=0.9, unrel=NULL,  contour=NULL, dof=1, control=NULL, debias="none", show=c(FALSE,FALSE)) {							
							
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
							
get2pbounds<-function()  {							
	ptDensity=120						
	if(!is.null(control$ptDensity))  ptDensity<-control$ptDensity[1]
	RadLimit<-1e-5
	if(!is.null(control$RadLimit))	RadLimit<-control$RadLimit
	# if(missing(contour))  {						
	# Error in missing(contour) : 'missing' can only be used for arguments						
	if(length(contour)==0) {						
		contour<-MLEcontour(x, dist, CL, dof=dof, ptDensity=ptDensity, RadLimit=RadLimit, debias=debias)					
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
							
	if(show[1]==TRUE)  {						
		plot(xvals,ypts,type="l")					
		lines(outmat[2,],outmat[1,],col="red")					
		lines(outmat[5,],outmat[1,],col="blue")					
	}						
							
							
							
							
	outDF<-data.frame(percentile=dp*100, lower=exp(outmat[2,]), datum=exp(xvals), upper=exp(outmat[5,]))						
	outList<-list(bounds=outDF,contour=contour)						
	outList						
}							
## close the get2pbounds function definition							
							
if(npar==2) {							
	outList<-get2pbounds()						
}							
							
# it will be cleaner to place 3p code in its own if block rather than just an else							
if(npar==3) {							
	modx<-function(x, tz) {						
		left<- ifelse(x$left<0, 0, x$left-tz)					
		right<- ifelse(x$right>0, x$right-tz, -1)					
		qty<-x$qty					
		return(data.frame(left=left, right=right, qty=qty))					
	}						
							
	ptDensity<-100						
	if(!is.null(control$ptDensity))  ptDensity<-control$ptDensity						
	tzpoints<-c(10,10,1)						
	if(!is.null(control$tzpoints))  tzpoints<-control$tzpoints						
	RadLimit<-1e-5						
	if(!is.null(control$RadLimit)) RadLimit<-control$RadLimit						
	maxit<-10						
	if(!is.null(control$maxit)) maxit<-control$maxit						
							
	#usage mlefit(x, dist="weibull", npar=2, debias="none", optcontrol=NULL)						
	MLEfit<-mlefit(x, dist, npar=3)						
	## must extract these values as vectors, numbered items return dataframes						
	P2_opt<-as.vector(MLEfit[2])						
	P1_opt<-as.vector(MLEfit[1])						
	t0_opt<-as.vector(MLEfit[3])						
	MLLx3p<-as.vector(MLEfit[4])						
							
	ratioLL  <-  MLLx3p- qchisq(CL,dof)/2						
	if( !is.null(attr(MLEfit,"message"))) {						
		if(attr(MLEfit, "message") == "t0 cutoff at minimal change") {					
		warning("t0 cutoff at minimal change")					
		}					
	## need to determine what to do when cutoff is LL at minimal change, characteristic of negative t0 instability						
		if(attr(MLEfit, "message") == "optimum not found, t0 cutoff at minimal gof change") {					
		warning("optimum not found, t0 cutoff at minimal gof change")					
		}					
	}						
							
	if(!is.null(attr(MLEfit, "rebound")))  {						
		maxtz<-attr(MLEfit,"rebound")					
	}else{						
		## Note discoveries continue to be discoveries until x$right-tz = zero					
		maxtz<-min(x$right[x$right != -1])					
	}						
							
	if(t0_opt==0) {						
	outList<-get2pbounds()						
	# control then drops to end of else block						
	# break would have not meaning here					
	}else{												
			tzpoints_arg<-tzpoints[1]				
			tzpoints_now<-tzpoints[1]				
			valid_tz<-NULL				
			invalid_tz<-NULL				
							
		repeat{					
			if(t0_opt> 0 )  {				
				if(!is.null(attr(MLEfit, "rebound")))  {			
					tzvec<-seq(0, maxtz, by=maxtz/tzpoints_now)		
				}else{			
					tzvec<-seq(0, maxtz-maxtz/tzpoints_now, by=maxtz/tzpoints_now)		
				}			
				if(!is.null(invalid_tz)) {			
					tzvec<-unlist(sapply(tzvec, function(X) if(X>max(invalid_tz)) X))		
				}			
			}else{				
			## just a guess for negative t0 will actually give one extra point beyond t0_opt				
			##	tzvec<-seq(0, t0_opt+t0_opt/tzpoints, by=t0_opt/tzpoints)			
				tzseq<-seq(maxtz-maxtz/tzpoints_now, t0_opt, by=(t0_opt-maxtz)/tzpoints_now)			
			## remove 0 if it is within this sequence, so it can be made first entry and appear only once				
				tzseq<-tzseq[!tzseq %in% 0]			
				tzvec<-c(0, seq(maxtz-maxtz/tzpoints_now, t0_opt, by=(t0_opt-maxtz)/tzpoints_now))			
				if(!is.null(invalid_tz)) {			
					tzvec<-unlist(sapply(tzvec, function(X) if(X<min(invalid_tz)) X))		
				}			
			}				
							
			for(tz in 1:length(tzvec)) {				
				## get the mle for 2-parameter fit at tz for validation testing			
				mle2p<-mlefit(modx(x, tzvec[tz]), dist)			
				if(mle2p[3]<ratioLL) {			
					invalid_tz<-unique(c(invalid_tz,tzvec[tz]))		
				}else{			
					valid_tz<-unique(c(valid_tz,tzvec[tz]))		
				}			
			}				
							
			if(length(valid_tz)<tzpoints_arg*.7) {				
				tzpoints_now<-tzpoints_now * 2			
				if(tzpoints_now>100) break			
			}else{				
				break			
			}				
		}					
							
		## add negative t0 candidates to valid_tz if zero is a valid_tz					
		if((t0_opt> 0 && 0 %in% valid_tz) )  {					
			tzpoints_neg<-seq(tzpoints[3],tzpoints[2], by=tzpoints[3])				
			step_size<-valid_tz[2]				
							
			for(try in tzpoints_neg) {				
				## get the mle for 2-parameter fit at tz for validation testing			
				neg_tz<- (-1)*step_size*try			
				mle2p<-mlefit(modx(x, tzvec[tz]), dist)			
				if(mle2p[3]<ratioLL) {			
					invalid_tz<-unique(c(invalid_tz,neg_tz))		
				}else{			
					valid_tz<-unique(c(valid_tz,neg_tz))		
				}			
			}				
		}					
							
		if(t0_opt<0) {					
			tzpoints_neg<-seq(tzpoints[3],tzpoints[2], by=tzpoints[3])				
			step_size<-valid_tz[length(valid_tz)]-valid_tz[length(valid_tz)-1]				
							
			for(try in tzpoints_neg) {				
				## get the mle for 2-parameter fit at tz for validation testing			
				neg_tz<- t0_opt+step_size*try			
				mle2p<-mlefit(modx(x,neg_tz), dist)			
				if(mle2p[3]<ratioLL) {			
					invalid_tz<-unique(c(invalid_tz,neg_tz))		
				}else{			
					valid_tz<-unique(c(valid_tz,neg_tz))		
				}			
			}				
		}					
							
		## Need to add t0_opt into valid_tz, but no need to sort?					
		valid_tz<-c(valid_tz,t0_opt)					
							
		# initialize objects to be constructed in the tz loop					
		contour_list<-list()					
		bounds_list<-list()					
		xlb_mat<-xub_mat<-NULL					
		# get seed parameters for contour range determination					
		fit2p<-mlefit(modx(x, t0_opt), dist)					
		minP1<-maxP1<-fit2p[1]					
		minP2<-maxP2<-fit2p[2]					
							
		list_item=0					
	## this is the super loop through all tz's						
	for(tz in valid_tz)  {						
			list_item=list_item+1				
		## get adjusted contour points for this modified.by.tz data					
		## can use MLEcontour in place of  test_contour3 with WeibullR version >= 1.0.10.3					
		##	contourpts<-test_contour3(x-tz, s-tz, MLLx=MLLx3p)				
			contourpts<-MLEcontour(modx(x, tz), dist, MLLx=MLLx3p, debias=debias)				
		## and then adjust for the modified.by.t0 basis					
			if(dist=="weibull") mod_P1<-contourpts[,1]+tz-t0_opt				
			if(dist=="lognormal") mod_P1<-log(exp(contourpts[,1])+tz-t0_opt)				
			contourpts<-data.frame(mod_P1,P2=contourpts[,2])				
		## while we are here let's get extents for a contour plot					
			maxP2<-max(c(maxP2,contourpts[,2]))				
			minP2<-min(c(minP2, contourpts[,2]))				
			minP1<-min(c(minP1, contourpts[,1]))				
			maxP1<-max(c(maxP1, contourpts[,1]))				
		## get points for the bounds for this tz case (on a modified.by.t0 basis)					
		##  using LRbounds directly from WeibullR					
			boundpts<-LRbounds(modx(x,t0_opt), dist, unrel=dp, contour=contourpts, CL=CL, dof=dof, debias=debias)				
		## assemble matrices for all the bounds, so we can find min and max values later					
			xlb_mat<-rbind(xlb_mat, boundpts$bounds$lower)				
			xub_mat<-rbind(xub_mat, boundpts$bounds$upper)				
		## finally assemble  lists for all tz bounds and contours					
			contour_list[[list_item]]<-boundpts$contour				
			bounds_list[[list_item]]<-boundpts$bounds				
	}						
							
							
	boundsDF<-data.frame(percentile=boundpts$bounds$percentile,						
		lower=apply(t(xlb_mat),1,min)+t0_opt,					
		datum=boundpts$bounds$datum+t0_opt,					
		upper=apply(t(xub_mat),1,max)+t0_opt 					
	)						
	
	if(debias!="none")  {	
		attr(boundsDF,"bias_adj")<-debias
	}	
	
	ylo<-floor(minP2)						
	if(maxP2<(floor(maxP2)+.5)) {						
		yhi<-floor(maxP2)+.5					
	}else{						
		yhi<-floor(maxP2)+1					
	}						
							
	P1Dec<-10^(floor(log(minP1)/log(10))-1)						
	xlo<-P1Dec*(floor(minP1/P1Dec)-1)						
	xhi<-P1Dec*(floor(maxP1/P1Dec)+1)						
							
	contour_range<-list(xlim=c(xlo,xhi),ylim=c(ylo,yhi))						
							
	outList<-list(bounds=boundsDF,bounds_list=list(xlb_mat,xub_mat),						
		 contour_list=contour_list,					
		 contour_range=contour_range					
	)						
							
	if(show[1] && show[2]) {						
## set the graphic device for double plot output							
## order of numbers in first arg to matrix determines order top to bottom of plots							
## additional values represent relative width and height of graphic panes							
		layout(matrix(c(2,1),1,2, byrow=TRUE))					
		layout.show(n=2)					
	}						
		if(show[1]) {					
			if(!exists("callWblr")) {				
				callWblr<-function(data, canvas) {			
					# data is in mleframe form		
					tefail<-data[which((data$right-data$left)==0 ),]		
					tesusp<-data[which(data$right==-1),]		
					interval_data<-data[which((data$right-data$left)>0),]		
					teq<-data.frame(time=NULL, event=NULL, qty=NULL)		
					if(dim(tefail)[1]>0) teq<-rbind(teq,data.frame(time=tefail$left, event=1, qty=tefail$qty))		
					if(dim(tesusp)[1]>0) teq<-rbind(teq,data.frame(time=tesusp$left, event=0, qty=tesusp$qty))		
					if(dim(teq)[1]==0) teq=NULL		
					if(dim(interval_data)[1]==0) interval_data=NULL		
					return(wblr(teq, interval=interval_data, canvas=canvas))		
				}			
			}				
							
		if(!exists("p2y")) {					
		p2y <- function(p,canvas="weibull"){					
			# This is the inverse Cumulative Distribution function				
			# used to transform a probability value to the				
			# y-axis of the plot canvas. 				
			# Use of this transformation permits distributions				
			# to appear as curves on unrelated canvas				
			if(canvas =="weibull")ret <- log(qweibull(p,1,1))				
			if(canvas =="lognormal") ret <- qlnorm(p,0,1)				
			ret				
		}					
		}					
							
			mod.obj<-callWblr(modx(x,t0_opt), canvas=dist)				
			mod.obj<-wblr.fit(mod.obj, method.fit="mle", dist=dist)				
			mod.obj<-wblr.conf(mod.obj,method.conf="lrb",lwd=2, lty=2,col="orange")				
			if(show[1] && show[2]) {				
				plot(mod.obj, xlab="time - t0",			
				canvas = dist,			
				in.legend=FALSE,			
				main="Modified Data Plot")			
			}else{				
				plot(mod.obj, xlab="time - t0", in.legend.blives=FALSE,			
				canvas = dist,			
				main="Modified Data Plot")			
			}				
			lines(boundsDF$lower-t0_opt,p2y(boundsDF$percentile/100,canvas = dist),lwd=2,col="red")				
			lines(boundsDF$upper-t0_opt,p2y(boundsDF$percentile/100,canvas = dist),lwd=2,col="red")				
		}					
		if(show[2]) {					
			C2P<-contour_list				
			plotargs <- list(x=NA,axes=TRUE)				
			plotargs$xlim <- contour_range$xlim				
			plotargs$ylim <- contour_range$ylim				
			plotargs$main <- "Contour Plots"				
			plotargs$sub  <- "data modified by potential t0 values"				
			plotargs$log <- ""				
			plotargs$xlab <- names(C2P[[1]])[1]				
			plotargs$ylab <- names(C2P[[1]])[2]				
							
			do.call("plot.default",plotargs)				
							
            abline(							
                h=pretty(contour_range$ylim,10),							
                v=pretty(contour_range$xlim,10),							
                col = "grey")							
							
			for(cntr in 1:length(C2P) )  {				
			# plot the contours				
				lwd<-1			
				lty<-1			
				col<-"black"			
				if(cntr==1&&!(0 %in% invalid_tz)) {col="darkgreen";lwd=2}			
				if(cntr==length(C2P)) {col="orange";lwd=2;lty=2}			
				points(C2P[[cntr]],type="l", lwd=lwd, lty=lty, col=col)			
							
			}				
		}					
							
	if(show[1] && show[2]) {						
	## reset the graphics device for single output						
		par(mfrow=c(1,1))					
	}						
							
	} # close else block for t0_opt!=0						
} # close the 3p code block							
							
							
							
							
outList							
}							
