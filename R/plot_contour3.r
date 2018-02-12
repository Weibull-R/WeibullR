## plot_contour2.r

## inspired by contour.abrem originally authored by Jurgen Symynck, April 2014
## Extensive re-write by David J. Silkworth removes the S3 functionality that never applied
## to this data format which specifies a single level of contour in only x,y coordinates.

## This rewrite is intended to add ability to calculate contours based on the data contained in
## a wblr object list, but permitting specification of distribution, multiple confidence levels,
## and degrees of freedom, the latter of which is important when contours are used for comparison
## of datasets.



plot_contour <- function(x ,CL=NULL ){
if(!is.null(CL)) {
if(class(CL)=="wblr") stop("multiple wblr objects must be entered as a list")
}
	    if(identical(class(x),"wblr")) x <- list(x)
	    if(!all(sapply(x,function(x)identical(class(x),"wblr")))){
	        stop("Argument \"x\" is not of class \"wblr\" or ",
	        "a list of \"wblr\" objects.")
	    }
	    # as of this point, x is always a list of one or more wblr objects

    # get options list from first object
## opa options appear to set main.contour and sub.contour, which could have better been added directly as arguments to plot_contour
## also used to set MLEpoint graphic parameters as cex.points and lwd.points 
## temporarily setting lty and lwd for contour lines 
	    opa <- x[[1]]$options




	if(is.null(CL))  {
		C2P<-ExtractContoursFromObjects(x)
	}else{
		C2P<-CalculateContours(x, CL)
	}

	## it could be desirable to exponentiate the log parameters of lognormal contours here.

	    # +--------------------------+
	    # |  create new plot canvas  |
	    # +--------------------------+
	    contourRanges <- findContourRanges(C2P)
	    if(!is.null(contourRanges)){
	        xlimits <- range(contourRanges[,1])
	        ylimits <- range(contourRanges[,2])
	 opanames <- names(opa)
	        plotargs <- c(list(x=NA,axes=TRUE),
	  opa[opanames %in% plot_default_args()])
	        plotargs$xlim <- xlimits
	        plotargs$ylim <- ylimits
	        plotargs$main <- opa$main.contour
	        plotargs$sub  <- opa$sub.contour
	        plotargs$log <- ""
	        plotargs$xlab <- names(C2P[[1]]$contour)[1]
	        plotargs$ylab <- names(C2P[[1]]$contour)[2]

	        do.call("plot.default",plotargs)
	        if(opa$is.plot.grid){
	            abline(
	                h=pretty(contourRanges[,2],10),
	                v=pretty(contourRanges[,1],10),
	                col = opa$col.grid)
	                # TODO: add userchoice in grid density here
	        }
	    }else message("plot.contour: No contours available in (list of) wblr objects.")


	    # +------------------+
	    # |  plot contours   |
	    # +------------------+

	for(cntr in 1:length(C2P) )  {

	# plot MLE points same options as data points in first object
		points(x=C2P[[cntr]]$MLEpt[1],y=C2P[[cntr]]$MLEpt[2],
#			pch=opa$options$pch,
			pch=opa$pch,
			#col=C2P[[cntr]]$color,
			col=opa$col,
			lwd=opa$lwd.points,
			cex=opa$cex.points)
#browser()
	# plot the contours
		points(C2P[[cntr]]$contour,type="l",
			lwd=opa$lwd,
			lty=opa$lty,
			col=C2P[[cntr]]$color)

	}
return(C2P)
}

getContoursFromSingleObject<-function(wblr) {

FOUND=FALSE

	if(!is.null(wblr$fit)) {
## a fit list exists

		for(fit_num in 1:length(wblr$fit)) {
			fit<-wblr$fit[[fit_num]]
			if(!is.null(fit$conf)) {
				for(conf_num in 1:length(fit$conf)) {
					conf<-fit$conf[[conf_num]]
					if(!is.null(conf$contour)) {
	## Yeah!, we found a contour
#						if(!exists("CP")) {
#use of exists("CP") became problematic when CP indeed existed outside of function environment
						if(!FOUND) {
							FOUND=TRUE
							CP<-list()
							j=1
							CP[[j]]<-list()

						}else{

				## here is the place to confirm fit types
							if(names(CP[[length(CP)]]$MLEpt[1])!= names(fit$MLEfit[1])) {
								warning("contours of mixed fit type found, mismatch ignored")
								break
							}

						j<-length(CP)+1
						CP[[j]]<-list
						}

						CP[[j]]$contour<-conf$contour
						CP[[j]]$MLEpt<-fit$MLEfit[-3]
						if(!is.null(conf$options$col)) {
							CP[[j]]$color<-conf$options$col
						}else{
							if(!is.null(fit$options$col)) {
								CP[[j]]$color<-fit$options$col
							}else{
								CP[[j]]$color<-wblr$options$col
							}
						}
					}
				}
			}
		}
	}
	CP
}

ExtractContoursFromObjects<-function(x)  {
	CP<-list()
	for(obj in 1:length(x))  {
		objCP<-getContoursFromSingleObject(x[[obj]])
		if(obj ==1) {
			for(li in 1:length(objCP))  {
				CP[[li]]<-objCP[[li]]
			}
		}else{
			if(names(CP[[length(CP)]]$contour[1])!=names(objCP[[1]]$contour[1])) {
				stop("dist mismatch in contours from objects")
			}
			cp_num<-length(CP)
			for(li in 1:length(objCP))  {
				cp_num<-cp_num+1
				CP[[cp_num]]<-objCP[[li]]
			}
		}
	}

	CP
}

contourRange <- function(contour){
	ra<-contour
	data.frame(range(ra[,1]),range(ra[,2]))
}


findContourRanges <- function(cplist) {
	findrange <- function(cp){
		 contourRange(cp$contour)
	}
	do.call("rbind",lapply(cplist,findrange))
}

ExtractContourParamsFromObject<-function(wblr) {
	str_eval=function(x) {return(eval(parse(text=x)))}
	getParam<-function(par_name) {
		val<-NULL
		if( !is.null(str_eval(paste0('conf$options$',par_name))) )  {
			val<-str_eval(paste0('conf$options$',par_name))
		}else{
			if( !is.null(str_eval(paste0('fit$options$',par_name))) )  {
				val<-str_eval(paste0('fit$options$',par_name))
			}else{
				if( !is.null(str_eval(paste0('wblr$options$',par_name))) )  {
					val<-str_eval(paste0('wblr$options$',par_name))
				}
			}
		}
		val
	}


	if(!is.null(wblr$fit)) {
## a fit list exists

		for(fit_num in 1:length(wblr$fit)) {
			fit<-wblr$fit[[fit_num]]
			if(!is.null(fit$conf)) {
				for(conf_num in 1:length(fit$conf)) {
					conf<-fit$conf[[conf_num]]
					if(!is.null(conf$contour)) {
	## Yeah!, we found a contour
						dist<-getParam("dist")
						dof<-getParam("dof")
						col<-getParam("col")
						lty<-getParam("lty")
						lwd<-getParam("lwd")
					}
				}
			}
		}
	}
	outlist<-NULL
	if(exists("dist")) {
		outlist<-list(dist=dist,dof=dof,col=col,lty=lty,lwd=lwd)
	}
	outlist
}

CalculateContours<-function(x, CL)  {
	c2p<-list()
	wblr_num<-0
	while(wblr_num < length(x))  {
		wblr_num<-wblr_num+1
		params<-ExtractContourParamsFromObject(x[[wblr_num]])
		fit<-unname(mlefit(x[[wblr_num]]$data$lrq_frame, dist=params$dist))
		
		for(cl_num in 1:length(CL))  {
## ptDensity could be a function of CL[cl_num]
## 360 for CL=.9, 40 for CL=.1
			dens<-ceiling(360*CL[cl_num]/.9)

			if(wblr_num==1 && cl_num==1)  {
				c2p_num<-1
			}else{
				c2p_num<-length(c2p)+1
			}
## usage MLEcontour(x,  dist="weibull", CL=0.9,dof=1,MLEfit=NULL, RadLimit=1e-5,
##		ptDensity=120, debias="none", show=FALSE)  {
			c2p[[c2p_num]]<-list()
			c2p[[c2p_num]]$contour<-MLEcontour(
					x[[wblr_num]]$data$lrq_frame, 
					dist=params$dist, 
					CL=CL[cl_num], 
					dof=params$dof, 
					MLEfit=fit,
					ptDensity=dens
					)
			c2p[[c2p_num]]$MLEpt<-fit
			c2p[[c2p_num]]$color<-params$col
# check implementation of ExtractParamsFromObject here
#			c2p[[c2p_num]]$lty<-params$lty
#			c2p[[c2p_num]]$lwd<-params$lwd
		}
	}
	c2p
}