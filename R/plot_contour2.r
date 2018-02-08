## plot_contour2.r

## inspired by contour.abrem originally authored by Jurgen Symynck, April 2014
## Extensive re-write by David J. Silkworth removes the S3 functionality that never applied
## to this data format which specifies a single level of contour in only x,y coordinates.

## This rewrite is intended to add ability to calculate contours based on the data contained in
## a wblr object list, but permitting specification of distribution, multiple confidence levels,
## and degrees of freedom, the latter of which is important when contours are used for comparison
## of datasets.



plot_contour <- function(x ,dist=NULL,CL=NULL, dof=NULL,...){
if(!is.null(dist)) {
if(class(dist)=="wblr") stop("multiple wblr objects must be entered as a list")
}
	    if(identical(class(x),"wblr")) x <- list(x)
	    if(!all(sapply(x,function(x)identical(class(x),"wblr")))){
	        stop("Argument \"x\" is not of class \"wblr\" or ",
	        "a list of \"wblr\" objects.")
	    }
	    # as of this point, x is always a list of one or more wblr objects

    # get options list from first object
	    opa <- x[[1]]$options
## handling of dots to be performed later. Dots arguments cannot override fixed args
## but want to permit revision of opa$main.contour and opa$sub.contour
## this might be the way to introduce debias options for calculated contours?



	if(length(dist)==0 && length(CL)==0 && length(dof)==0 )  {
		C2P<-ExtractContoursFromObjects(x)
	## it could be desirable to exponentiate the log parameters of lognormal contours here.
	}
## }else{
## 		C2P<-CalculateContours(x, dist, CL, dof)   # yet to be written
## }

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
			pch=opa$options$pch,
			#col=C2P[[cntr]]$color,
			col=opa$col,
			lwd=opa$lwd.points,
			cex=opa$cex.points)

	# plot the contours
		points(C2P[[cntr]]$contour,type="l",
			lwd=opa$lwd,
			lty=opa$lty,
			col=C2P[[cntr]]$color)

	}

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

