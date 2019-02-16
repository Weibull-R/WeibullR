##File LRbounds3pw.r is placed in the WeibullR package temporarily during testing
## this content is expected to be incorporated into LRbounds.r in the future.


LRbounds3pw<-function(x, s=NULL, CL=0.9, DF=1 ,ptDensity=100, tzpoints=10, RadLimit=1e-5, show=FALSE)  {			
		## require WeibullR version >= 1.0.10.2 for alloydata	
		##require(WeibullR)	
		MLEfit<-mlefit(mleframe(x,s), npar=3)	
		## when the fit return was a dataframe must these vlues as vectors, 
		## now that the fit return is a vector must use numbered items
		Beta_opt<-MLEfit[2]	
		Eta_opt<-MLEfit[1]	
		t0_opt<-MLEfit[3]	
		MLLx3p<-MLEfit[4]	
		##ratioLL  <-  MLLx3p- qchisq(CL,DF)/2	
		if( !is.null(attr(MLEfit,"unstable"))) stop("unstable 3rd parameter fit")	
			
		if(!is.null(attr(MLEfit, "rebound_pt")))  {	
			maxtz<-attr(MLEfit,"rebound_pt")
		}else{	
			maxtz<-min(x)
		}	
			
		if(t0_opt==0) stop("t0 = 0, nothing to do")	
			
		if(t0_opt> 0 )  {	
			
			tzvec<-seq(0, maxtz-maxtz/tzpoints, by=maxtz/tzpoints)
		}else{	
		## just a guess for negative t0 will actually give one extra point beyond t0_opt	
			tzvec<-seq(0, t0_opt+t0_opt/tzpoints, by=t0_opt/tzpoints)
		}	
			
		## Need to add t0_opt into tzvec, but no need to sort?	
		tzvec<-c(tzvec,t0_opt)	
			
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
	for(tz in tzvec)  {		
			list_item=list_item+1
		## get adjusted contour points for this modified.by.tz data	
		## can use MLEcontour in place of  test_contour3 with WeibullR version >= 1.0.10.3	
		##	contourpts<-test_contour3(x-tz, s-tz, MLLx=MLLx3p)
			contourpts<-MLEcontour(mleframe(x-tz, s-tz), MLLx=MLLx3p)
		## and then adjust for the modified.by.t0 basis	
			mod_eta<-contourpts[,1]+tz-t0_opt
			contourpts<-data.frame(mod_eta,contourpts[,2])
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
		yhi<-floor(maxBeta)+1	
			
		EtaDec<-10^(floor(log(minEta)/log(10))-1)	
		xlo<-EtaDec*(floor(minEta/EtaDec)-1)	
		xhi<-EtaDec*(floor(maxEta/EtaDec)+1)	
			
			
			
		outlist<-list(bounds=boundsDF,bounds_list=bounds_list,	
			 contour_list=contour_list,
			 contour_range=list(xlim=c(xlo,xhi),ylim=c(ylo,yhi))
		)	
			
## some response to show=TRUE to be developed here	
		if(show==TRUE) {
		mod.obj<-wblr(x-t0_opt,s-t0_opt)
		mod.obj<-wblr.fit(mod.obj, method.fit="mle")
		mod.obj<-wblr.conf(mod.obj,method.conf="lrb",lwd=1, lty=2,col="red")
		plot(mod.obj, xlab="time - t0", main="Modified Data Plot")
		lines(boundsDF$lower,p2y(boundsDF$percentile/100),lwd=2,col="red")
		lines(boundsDF$upper,p2y(boundsDF$percentile/100),lwd=2,col="red")
		}
		
return(outlist)			
}			
