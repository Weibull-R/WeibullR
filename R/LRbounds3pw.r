##File LRbounds3pw.r is placed in the WeibullR package temporarily during testing
## this content is expected to be incorporated into LRbounds.r in the future.


LRbounds3pw<-function(x, s=NULL, CL=0.9, DF=1 ,ptDensity=120, tzpoints=c(10,10,1), RadLimit=1e-5, listout=FALSE, show=c(FALSE,FALSE))  {			
		if(class(x)=="data.frame") {		
	##Note that any qty column would be ignored here. 
	##Actually there should be an expansion by qty if this code remains.			
		s=x$time[x$event==0]	
		x=x$time[x$event==1]	
		}		

		MLEfit<-mlefit(mleframe(x,s), npar=3)	
		## now that the fit return is a vector must use numbered items
		Beta_opt<-MLEfit[2]	
		Eta_opt<-MLEfit[1]	
		t0_opt<-MLEfit[3]	
		MLLx3p<-MLEfit[4]	
		ratioLL  <-  MLLx3p- qchisq(CL,DF)/2	
		if( !is.null(attr(MLEfit,"message"))) {	
			if(attr(MLEfit, "message") == "t0 cutoff at minimal change") {
			warning("t0_opt at optcontrol limit, 2p bounds on modified data applies")
			}
		## need to determine what to do when cutoff is LL at minimal change, characteristic of negative t0 instability	
			
			
		}	
	
		if(!is.null(attr(MLEfit, "rebound")))  {	
			maxtz<-attr(MLEfit,"rebound")
		}else{	
			maxtz<-min(x)
		}	
			
		if(t0_opt==0) stop("t0 = 0, nothing to do")	
		
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
			mle2p<-mlefit(mleframe(x-tzvec[tz], s-tzvec[tz]))		
			if(mle2p[3]<ratioLL) {		
				invalid_tz<-unique(c(invalid_tz,tzvec[tz]))	
			}else{		
				valid_tz<-unique(c(valid_tz,tzvec[tz]))	
			}		
		}			
										
		if(length(valid_tz)<tzpoints_arg*.7) {			
			tzpoints_now<-tzpoints_now * 2
			if(tzpoints_now<1000) break		
		}else{			
			break		
		}			
	}	

	## add negative t0 candidates to valid_tz				
	if((t0_opt> 0 && 0 %in% valid_tz) )  {				
		tzpoints_neg<-seq(tzpoints[3],tzpoints[2], by=tzpoints[3])			
		step_size<-valid_tz[2]			
					
		for(try in tzpoints_neg) {			
			## get the mle for 2-parameter fit at tz for validation testing		
			neg_tz<- (-1)*step_size*try		
			mle2p<-mlefit(mleframe(x-neg_tz, s-neg_tz))		
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
			mle2p<-mlefit(mleframe(x-neg_tz, s-neg_tz))		
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
		fit2p<-mlefit(mleframe(x-t0_opt,s-t0_opt))
		minEta<-maxEta<-fit2p[1]
		minBeta<-maxBeta<-fit2p[2]
		

	
		list_item=0	
		## this is the super loop through all tz's	
	for(tz in valid_tz)  {	
			list_item=list_item+1
		## get adjusted contour points for this modified.by.tz data	
		## can use MLEcontour in place of  test_contour3 with WeibullR version >= 1.0.10.3	
		##	contourpts<-test_contour3(x-tz, s-tz, MLLx=MLLx3p)
			contourpts<-MLEcontour(mleframe(x-tz, s-tz), MLLx=MLLx3p)
		## and then adjust for the modified.by.t0 basis	
			mod_eta<-contourpts[,1]+tz-t0_opt
			contourpts<-data.frame(mod_eta,beta=contourpts[,2])
		## while we are here let's get extents for a contour plot	
			maxBeta<-max(c(maxBeta,contourpts[,2]))
			minBeta<-min(c(minBeta, contourpts[,2]))
			minEta<-min(c(minEta, contourpts[,1]))
			maxEta<-max(c(maxEta, contourpts[,1]))
		## get points for the bounds for this tz case (on a modified.by.t0 basis)	
		##  using LRbounds directly from WeibullR	
			boundpts<-LRbounds(mleframe(x-t0_opt,s-t0_opt), contour=contourpts, dof=DF)
		## assemble matrices for all the bounds, so we can find min and max values later	
			xlb_mat<-rbind(xlb_mat, boundpts$bounds$lower)
			xub_mat<-rbind(xub_mat, boundpts$bounds$upper)
		## finally assemble  lists for all tz bounds and contours	
			contour_list[[list_item]]<-boundpts$contour
			bounds_list[[list_item]]<-boundpts$bounds
	}		
			
			
		boundsDF<-data.frame(percentile=boundpts$bounds$percentile,	
			lower=apply(t(xlb_mat),1,min),
			datum=boundpts$bounds$datum,
			upper=apply(t(xub_mat),1,max) 
		)	
			
		ylo<-floor(minBeta)	
		if(maxBeta<(floor(maxBeta)+.5)) {
			yhi<-floor(maxBeta)+.5
		}else{
			yhi<-floor(maxBeta)+1	
		}
			
		EtaDec<-10^(floor(log(minEta)/log(10))-1)	
		xlo<-EtaDec*(floor(minEta/EtaDec)-1)	
		xhi<-EtaDec*(floor(maxEta/EtaDec)+1)

		contour_range<-list(xlim=c(xlo,xhi),ylim=c(ylo,yhi))		
						
		outlist<-list(bounds=boundsDF,bounds_list=list(xlb_mat,xub_mat),	
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
			mod.obj<-wblr(x-t0_opt,s-t0_opt)
			mod.obj<-wblr.fit(mod.obj, method.fit="mle")
			mod.obj<-wblr.conf(mod.obj,method.conf="lrb",lwd=2, lty=2,col="orange")
			if(show[1] && show[2]) {
				plot(mod.obj, xlab="time - t0",
				in.legend=FALSE,
				main="Modified Data Plot")
			}else{
				plot(mod.obj, xlab="time - t0", in.legend.blives=FALSE,
				main="Modified Data Plot")
			}
			lines(boundsDF$lower,p2y(boundsDF$percentile/100),lwd=2,col="red")
			lines(boundsDF$upper,p2y(boundsDF$percentile/100),lwd=2,col="red")
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

	if(listout) {
		return(outlist)	
	}else{
		return(boundsDF)
	}
}			
