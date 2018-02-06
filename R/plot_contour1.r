#plot_contour1.r

# This file presents a plot_contour function to the WeibullR package based on a direct adaptation of 
# function abrem::contour.abrem originally written by Jurgen Syminck. 
# Adaptations have included containment of non-exported code as originally contained in findContourRanges.R in the abrem package.
# Within abrem plot.contour was specified as an S3 function, overriding default function graphics::contour.
# However, this contour plot is structured differently than the expectation for a graphics::contour argument set.
# For this reason the WeibullR::plot.contour will not be declared as an S3 function, merely a direct call function
# within the package.

# Adaptation by David Silkworth



plot_contour <- function(x,...){
	if(identical(class(x),"wblr")) x <- list(x)
	if(!all(sapply(x,function(x)identical(class(x),"wblr")))){
		stop("Argument \"x\" is not of class \"wblr\" or ",
		"a list of \"wblr\" objects.")
	}
	# as of this point, x is always a list of one or more abrem objects



	# get options list from first object
	opa <- x[[1]]$options


	opa <- modifyList(opa, list(...))


	# +--------------------------+
	# |	 create new plot canvas	 |
	# +--------------------------+
#	 contourRanges <- findContourRanges(x,opa$verbosity)
	contourRanges <- findContourRanges(x)
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
		plotargs$xlab <- "Eta"
		plotargs$ylab <- "Beta"

		do.call("plot.default",plotargs)
		if(opa$is.plot.grid){
			abline(
				h=pretty(contourRanges[,2],10),
				v=pretty(contourRanges[,1],10),
				col = opa$col.grid)
				# TODO: add userchoice in grid density here
		}
	}else message("plot.contour: No contours available in (list of) wblr objects.")
#	 r <- seq.log(opa$xlim[1]/10,opa$xlim[2]*10,c(1,5))

	# +------------------+
	# |	 plot contours	 |
	# +------------------+
	plotContours <- function(wblr){
		if(!is.null(wblr$fit)){
			plotContours2 <- function(fit){
				if(!is.null(fit$options)){
					opafit <- modifyList(wblr$options,fit$options)
				}else{opafit <- wblr$options}
				is_MLE <- any(c("mle","mle-rba") %in% tolower(fit$options$method.fit))
#				 if(!is.null(fit$conf$blives)){
				if(!is.null(fit$conf)){
					plotContours3 <- function(conf){
						if(!is.null(conf$options)){
							opaconf <- modifyList(opafit,conf$options)
						}else{opaconf <- opafit}
						 opaconf <- modifyList(opaconf,list(...))
 #						 if(!is.null(blicon$MLEXContour)){
						if(!is.null(conf$contour)){
## This only applies to weibull
							if(all(c(!is.null(fit$beta),!is.null(fit$eta)))){
								points(x=fit$eta,y=fit$beta,pch=wblr$options$pch,col=wblr$options$col,
								lwd=wblr$options$lwd.points,cex=wblr$options$cex.points)
								points(conf$contour,type="l",lwd=opaconf$lwd,lty=opaconf$lty,col=opaconf$col)
							}
						}
					}
					#mtrace(plotContours3)
					do.call("rbind",lapply(fit$conf,plotContours3))
						# combine the ranges from all MLEXContours
						# found in the list of blicons
				}
			}
			do.call("rbind",lapply(wblr$fit,plotContours2))
				# combine the ranges from all MLEXContours
				# found in the list of fits
		}
	}
	if(!is.null(contourRanges)) lapply(x,plotContours)
	invisible()
}

## Adapted from code formerly found in abrem package findContourRanges.R

contourRange <- function(contour){
## no idea why this worked in abrem
    #ra <- do.call("rbind",MLEXContour)
	ra<-contour
    data.frame(range(ra[,1]),range(ra[,2]))
}


findContourRanges <- function(x){
# +---------------------------------------------------+
# |  find absolute maximum and minimum contour range  |
# |     over the (list of) wblr objects               |
# +---------------------------------------------------+
# x is always a list of wblr object(s)

	findrange1 <- function(wblr){
		if(!is.null(wblr$fit)){
			findrange2 <- function(fit){
				if(!is.null(fit$conf)){
					 findrange3 <- function(conf){
						if(!is.null(conf$contour)){
							# a contour is available
							contourRange(conf$contour)
						}
					}
					do.call("rbind",lapply(fit$conf,findrange3))

				}
			}
	# combine the ranges from all contours
	# found in the list of fits
			do.call("rbind",lapply(wblr$fit,findrange2))
		}
	}
	do.call("rbind",lapply(x,findrange1))
}
