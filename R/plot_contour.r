## plot_contour.r

## inspired by contour.abrem originally authored by Jurgen Symynck, April 2014
## Extensive re-write by David J. Silkworth removes the S3 functionality that never applied
## to this data format which specifies a single level of contour in only x,y coordinates.

## This rewrite is intended to add ability to calculate contours based on the data contained in
## a wblr object list, but permitting specification of multiple confidence levels.



plot_contour <- function(x ,CL=NULL, AL=TRUE,...){
if(!is.null(CL)) {
if(class(CL)=="wblr") stop("multiple wblr objects must be entered as a list")
}
	    if(identical(class(x),"wblr")) x <- list(x)
	    if(!all(sapply(x,function(x)identical(class(x),"wblr")))){
	        stop("Argument \"x\" is not of class \"wblr\" or ",
	        "a list of \"wblr\" objects.")
	    }
	    # as of this point, x is always a list of one or more wblr objects
		
	arg <- list(...)

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

## convert lognormal parameters to antilog (AL) if specified
##	antilog<-function(y) exp(y)
	if(names(C2P[[1]]$contour)[1]=="Mulog" && AL==TRUE)  {
		for(C2Pli in 1:length(C2P))  {
##			C2P[[C2Pli]]$contour<-antilog(C2P[[C2Pli]]$contour)
			C2P[[C2Pli]]$contour<-exp(C2P[[C2Pli]]$contour)
		}
## axis labels drawn from first contour names, but all have been converted
		names(C2P[[1]]$contour)<-c("MuAL", "SigAL")
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
			plotargs$main <- opa$main.contour
## overrides from the dots
			if(!is.null(arg$xlim)) plotargs$xlim<-arg$xlim
			if(!is.null(arg$ylim)) plotargs$ylim<-arg$ylim
			if(is.null(arg$main)) plotargs$main<-arg$main
			if(is.null(arg$sub)) plotargs$sub<-arg$sub
			
			
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

	# plot MLE points always a black 'x' symbol
		points(x=C2P[[cntr]]$MLEpt[1],y=C2P[[cntr]]$MLEpt[2],
#			pch=opa$options$pch,
			pch=4,
			#col=C2P[[cntr]]$color,
			#col=opa$col,
			col="black",
			lwd=opa$lwd.points,
			cex=opa$cex.points)
#browser()
			lwd=C2P[[cntr]]$lwd
			lty=C2P[[cntr]]$lty
			col=C2P[[cntr]]$color
## overrides from the dots
			if(!is.null(arg$lwd)) lwd<-arg$lwd
			if(!is.null(arg$lty)) lty<-arg$lty
			if(!is.null(arg$col)) col<-arg$col
#browser()
	# plot the contours
		points(C2P[[cntr]]$contour,type="l", lwd=lwd, lty=lty, col=col)

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
						CP[[j]]$dist<-wblr$options$dist
						CP[[j]]$MLEpt<-fit$MLEfit[-3]
						if(!is.null(conf$options$lty)) {
							CP[[j]]$lty<-conf$options$lty
						}else{
							if(!is.null(fit$options$lty)) {
								CP[[j]]$lty<-fit$options$lty
							}else{
								CP[[j]]$lty<-wblr$options$lty
							}
						}
						if(!is.null(conf$options$lwd)) {
							CP[[j]]$lwd<-conf$options$lwd
						}else{
							if(!is.null(fit$options$lwd)) {
								CP[[j]]$lwd<-fit$options$lwd
							}else{
								CP[[j]]$lwd<-wblr$options$lwd
							}
						}
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
	if(!exists("CP")) {
		stop("no contour found in listed object, try adding CL argument to calculate")
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

# don't really want to get parameter modifications from potentially multiple fits
# so this function is abandoned
#	getParam2<-function(par_name) {
#		val<-NULL
#		if( !is.null(str_eval(paste0('fit$options$',par_name))) )  {
#			val<-str_eval(paste0('fit$options$',par_name))
#		}else{
#			if( !is.null(str_eval(paste0('wblr$options$',par_name))) )  {
#			}
#		}
#		val
#	}

# no fit exists, so get param from base object options only, perhaps defaults
	getParam3<-function(par_name) {
		val<-NULL
		if( !is.null(str_eval(paste0('wblr$options$',par_name))) )  {
			val<-str_eval(paste0('wblr$options$',par_name))
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
## Yeah!, we found a contour, get parameters from here seeking back toward
## base object options list if necessary.
						dist<-getParam("dist")
						dof<-getParam("dof")
						col<-getParam("col")
						lty<-getParam("lty")
						lwd<-getParam("lwd")
					}else{
# a conf exists, but not a contour, so get the params from base object options
# for some reason extraction getParam2 did not work here, did not want to debug further
							dist<-getParam3("dist")
							dof<-getParam3("dof")
							col<-getParam3("col")
							lty<-getParam3("lty")
							lwd<-getParam3("lwd")
					}
				}
			}else{
# a fit exists, but get the params from base object options
# for some reason extraction getParam2 did not work here, did not want to debug further
					dist<-getParam3("dist")
					dof<-getParam3("dof")
					col<-getParam3("col")
					lty<-getParam3("lty")
					lwd<-getParam3("lwd")
			}
		}
	}else{
# no fit exists, get the params from base object options
			dist<-getParam3("dist")
			dof<-getParam3("dof")
			col<-getParam3("col")
			lty<-getParam3("lty")
			lwd<-getParam3("lwd")
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
# warn and drop any 3p suffix from params$dist
# this should never happen because 3p distribution specification should only appear in wblr.fit, not wblr
# there would be no conf in the 3p fit, so only base wblr dist option would be returned here
# but lets just avoid this strange case.
		if(substr(params$dist,nchar(params$dist)-1,nchar(params$dist))=="3p"){
		params$dist<-substr(params$dist,1,nchar(params$dist)-2)
		warning("3p dist modification specified in wblr has been ignored")
		}

# test for dist mismatch here
		if(wblr_num > 1) {
			if(c2p[[length(c2p)]]$dist!=params$dist) {
				stop("dist mismatch in entered objects")
			}
		}
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
			c2p[[c2p_num]]$dist<-params$dist
			c2p[[c2p_num]]$MLEpt<-fit
			c2p[[c2p_num]]$color<-params$col
# check implementation of ExtractParamsFromObject here
			c2p[[c2p_num]]$lty<-params$lty
			c2p[[c2p_num]]$lwd<-params$lwd
		}
	}
	c2p
}