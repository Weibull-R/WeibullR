## plot.wblr.R
## refactored from code originally authored by Jurgen Symynck, April 2014
## Copyright 2014-2017 OpenReliability.org
#
# For more info, visit http://www.openreliability.org/
#
# For the latest version of this file, check the Subversion repository at
# http://r-forge.r-project.org/projects/Weibull-R/
##-------------------------------------------------------------------------------
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

plot.wblr <- function(x,...){
    # +------------------------------+
    # |  move wblr objects to list  |
    # |      of wblr objects        |
    # +------------------------------+
    if(identical(class(x),"wblr")) x <- list(x)
    if(!all(sapply(x,function(x)identical(class(x),"wblr")))){
        stop("Argument \"x\" is not of class \"wblr\" or ",
        "a list of \"wblr\" objects.")
    }
    # as of this point, x is always a list of one or more wblr objects
    
    # +------------------------------------+
    # |  create default options arguments  |
    # +------------------------------------+
    opa <- x[[1]]$options
    opa <- modifyList(opa, list(...))


    
    # +--------------------------+
    # |  create new plot canvas  |
    # +--------------------------+
    ra <- findMaxDataRange(x,opa$log)
## negative failure times should never reach this code anyway. Perhaps a stop would be better
        # NA values can be part of ra, when log scales are to be used
        # and there are negative failure times
    xlimits <- range(ra$xrange,na.rm=TRUE)
    ylimits <- range(ra$yrange,na.rm=TRUE)
    if(is.null(opa$xlim)){
        opa$xlim <- c(10^(log10(xlimits[1])-0.5),
            10^(log10(xlimits[2])+1))
    }
    if(is.null(opa$ylim)){
        if(ylimits[1] < 0.01) opa$ylim <- c(signif(ylimits[1],1),0.99)
        else opa$ylim <- c(0.01,0.99)
        # do not care about the upper limit
    }
    opanames <- names(opa)
    plotargs <- c(list(x=NA,axes=FALSE),
        opa[opanames %in% plot_default_args()])
    if(!is.null(plotargs$ylim)){
        plotargs$ylim <- F0inv(plotargs$ylim,opa$log)
    }
    plotargs$main <- NULL
        # do not plot "main" just yet...
    if(!is.null(opa$mar))par(mar=opa$mar)
    if(!is.null(opa$mai))par(mai=opa$mai)
    do.call(plot.default,plotargs)
    if(opa$is.plot.grid){
        abline(
            h=F0inv(seq.wb(opa$ylim[1]/10,1-(1-opa$ylim[2])/10),opa$log),
            v=seq.log(opa$xlim[1]/10,opa$xlim[2]*10,seq(0,10,1)),
            col = opa$col.grid)
    }
    r <- seq.log(opa$xlim[1]/10,opa$xlim[2]*10,c(1,5))
    #lin <- 0.0
    for(t in c(1,3)){
        axis(t,at=seq.log(opa$xlim[1]/10,opa$xlim[2]*10,seq(0,10,0.2)),
            labels=NA,tcl=-0.25)#,line=0.0
            # plot top and bottom axis tickmarks
        axis(t,at=r,labels=r,tcl=-0.75)#,line=0.0
            # plot top and bottom axis labels
    }
    r <- c(seq.wb(opa$ylim[1]/10,1-(1-opa$ylim[2])/10,c(1,2,5)),0.9)
    for(t in c(2,4)){
        # TODO: rewrite as do.call() or apply()
        axis(t,at=F0inv(seq.wb(opa$ylim[1]/10,1-(1-opa$ylim[2])/10),
            opa$log),labels=NA,tcl=-0.25)#,line=0.0
            # plot left and right axis tickmarks
        axis(t,at=F0inv(r,opa$log),
            labels=r*100,tcl=-0.75)#,line=0.0
            # plot left and right axis labels
    }
    abline(h=0,lty = 3,col = opa$col.grid)
    title(main=opa$main,line=3)
    # plot the 63.2 [%] rank line

    # +--------------------------+
    # |  plot confidence bounds  |
    # +--------------------------+
    plotConfs <- function(wblr){
        #opadata <- modifyList(x$options, arg)
        if(!is.null(wblr$fit)){
            ret <- lapply(wblr$fit,plotConfsInFit,opadata=wblr$options,...)
        }else{
            message("plotConfs: This wblr object contains no fits" )
        }
    }
    lapply(x,plotConfs)

    # +-------------+
    # |  plot fits  |
    # +-------------+
    plotFits <- function(wblr){
        opadata <- modifyList(wblr$options,list(opa$xlim,opa$ylim))
        if(!is.null(wblr$fit)){
            ret <- lapply(wblr$fit,plotSingleFit,
                opadata=opadata,...)
        }else{
            warning("plotFits: This wblr object contains no fits.")
        }
    }
    lapply(x,plotFits)

    # +----------------------------------+
    # |  plot points and interval lines  |
    # +----------------------------------+
    plotSingleDataSet <- function(x){
        if(opa$is.plot.pp){
		if(!is.null(x$data)) {
            # TODO: possibly, this does not allow much flexibility in plotting.
            opadata <- modifyList(x$options,list(...))

##	Removal of threshold option effects, t0 was not used here anyway		
	
## it is poor coding practice to place object assignments inside meaningless conditional test.			 
            if(!is.null(x$data$dpoints)){
## Note: x$data objects hava already been modified by t0 as may have been called for		
				ti <- x$data$dpoints$time
				pos <- x$data$dpoints$ppp
## ra was earlier used as the output of findMaxDataRange this is very poor/dangerous coding practice

## note how the log option is used to control F0inv here, defining lognormal versus weibull canvas
                #points(ti-t0,F0inv(pos,opadata$log),pch = opadata$pch,
## Note: x$data objects hava already been modified by t0 as may have been called for
                points(ti,F0inv(pos,opadata$log),pch = opadata$pch,
                    col = opadata$col,lwd = opadata$lwd.points,cex=opadata$cex.points)
                    # option "log" should only be set and read from either
                    # the arguments of plot.wblr
                    # Other instances should be ignored
            }
			
			if(!is.null(x$data$dlines)){
				dlines<-x$data$dlines
## Note chart x limits are  in opa$xlim[1]/10 and opa$xlim[2]*10
				for(intline in 1:dim(dlines)[1])  {
					if(dlines$t1[intline]==0) {
						x0=opa$xlim[1]/10
					}else{
						#x0=log(dlines$t1[intline])
						# log on x axis must have been specified elsewhere?
						x0=dlines$t1[intline]
					}
## this was specific for weibull canvas only				
					 	##y0=log(log(1/(1-dlines$ppp[intline])))
						y0=F0inv(dlines$ppp[intline],opadata$log)
						# log on x axis is specified in plotargs where log="xy" or log="x"
						#x1=log(dlines$t2[intline])
						x1=dlines$t2[intline]
						segments(x0, y0, x1, y0,
							col = opadata$interval.col,
							lty=opadata$interval.lty,
							lwd=opadata$interval.lwd
						)
				}			
			}
		}else{stop("no data to plot")}
        }
    }
    lapply(x,plotSingleDataSet)


    # +----------------+
    # |  plot legends  |
    # +----------------+
    lolegends <- NULL
    buildListOfLegends <- function(wblr){
        ret <- NULL

        if(!is.null(wblr$fit) && any(sapply(wblr$fit,function(fi)!is.null(fi)))){

            ret <- lapply(wblr$fit,buildSingleFitLegend,
                opadata=wblr$options,...)
        }else{
            if(wblr$options$is.plot.legend && opa$is.plot.legend){
                ret <- list(buildSingleDataLegend(wblr,opadata=wblr$options,...))
                if(!is.null(opa)) warning("buildListOfLegends: This wblr object contains no fits.")
            }
        }
        ret
    }
    lolegends <- unlist(lapply(x,buildListOfLegends),FALSE)
        # TODO: likely, unlist is NOT the best way to go here, investigate
    lolegends <- lolegends[sapply(lolegends,function(lol)!is.null(lol))]
        # omit any list entries that contain only NULL


    plotSingleLegend <- function(le,x,y){
        if(identical(label <- le$label,""))label <- NULL
        if(is.null(le$legend))le$legend <- ""
        legend(
            x=x,
            y=y,
            legend=le$legend,
            title=label,

            cex = le$legend.text.size,
            bg = "white",
            lty = unlist(le$lty),
            lwd = unlist(le$lwd),
            pch = unlist(le$pch),
            col = unlist(le$col),

            text.col = "black",
            xpd=TRUE

            )
            # TODO: Warning: unlist coerces numeric colors to character!
    }
    #if(!is.null(lolegends)){
    if(!is.null(lolegends) && any(sapply(lolegends,function(lol)!is.null(lol)))){
        lx <- rep(lolegends[[1]]$rect$left,length(lolegends))
        ly <- lolegends[[1]]$rect$top +
            c(0,cumsum(sapply(lolegends,function(le)le$rect$h)[-1]))
        if(opa$log %in% c("x","xy","yx")) lx <- 10^lx
        if(opa$log %in% c("y","xy","yx")) ly <- 10^ly
            # TODO: F0(ly): looks very suspicious that this works -> investigate!
        for(i in 1:length(lolegends)){
            plotSingleLegend(lolegends[[i]],lx[i],ly[i])
            # TODO: replace with lapply
        }
    }else{
        warning("plot.wblr: There is no legend to plot.")
    }

    invisible()
        # TODO: return the wblr object with updated graphical options
        # TODO: check if this makes sense when supplying a list
}

F0 <- function(q)
   1-exp(-exp(q))

F0inv <- function(p,log="x"){
    # transformation function to plot its argument
    # on the y-axis of the Weibull plot. This transformation function
    # lets the Weibull curve appear as a straight line on the weibull paper
    #
    # This is also the inverse Cumulative Distribution function of the
    # standardized Weibull plot with beta=eta=1
    # comparing  both implementationss of F0inv() with
    # system.time() does not show any significant difference
    #   log(log(1/(1-p)))}
    if(log %in% c("x",""))ret <- log(qweibull(p,1,1)) else ret <- qlnorm(p,0,1)
    ret
}

findMaxDataRange <- function(x,log=""){
    # +-------------------------------------------+
    # |  find absolute maximum and minimum range  |
    # |     over the (list of) wblr objects      |
    # +-------------------------------------------+
    # x is always a list of wblr object(s)
    findrange <- function(wblr){
        if(!is.null(wblr$data)){
			alltimes<-c(wblr$data$dpoints$time, wblr$data$dlines$t1[wblr$data$dlines$t1>0],wblr$data$dlines$t2)			
            if(!is.null(alltimes)){
				## wblr$data$dpoints$time constructed by getppp should never have na's
                ret <- data.frame(xrange=range(alltimes))
            }else{
                stop("no time data, cannot create plot canvas.")
            }
			allppp<-c(wblr$data$dpoints$ppp, wblr$data$dlines$ppp)
            if(!is.null(allppp)){
				## wblr$data$dpoints$ppp constructed by getppp should never have na's
                ret <- cbind(ret,yrange=range(allppp))
            }else{
                stop("$data$dpoints or $dlines contains no ppp -> ",
                    "cannot create plot canvas.")
            }
        }else{stop('Argument \"x\" contains no \"$data \" list object.')}
        ret
    }
#    if(identical(class(x),"wblr")){
#        if(v>= 2)message(match.call()[[1]],
#            ": Argument \"x\" is a single wblr object...")
#        ret <- findrange(x)
#    }else{
    if(all(sapply(x,function(x)identical(class(x),"wblr")))){
        ret <- do.call("rbind",lapply(x,findrange))
    }else{
        stop("Argument \"x\" is no list of \"wblr\" objects.")
    }
    # TODO: the above still needed? because x is always list of wblrs?
	
## this code should be elimintated, perhaps a stop if any zero or negative found	
    if(tolower(log) %in% c("x","xy","yx")){
        # if log scale is to be used then omit zero and negative time values 
        # from the range dataset
## first of all zero and negative values should never reach this point
## as they should have been elimintated from dpoints and are eliminated from alltimes above
## this code did not replace with a valid value anyway.
        ret[ret$xrange <=0,1] <- NA
    }
    ret
}


seq.wb <- function(from,to,base=seq(1,9,1)){
   # define gridline positions for 'Median Rank' axis (= y axis)
   r <- c(seq.log(from,.9,base),rev(1-seq.log(1-to,0.1,base))[-1])
   r[r >= from & r<=to]}

seq.log <- function(from,to,base=c(1,2,5)){
   r <- NULL
   for(i in floor(log10(from)):floor(log10(to)))r <- c(r,base*10^i)
   r[r >= from & r<=to]}
   
plotSingleConfBound <- function(blc,opafit,...){
    if(!is.null(blc$options)){
        opaconf <- modifyList(opafit,blc$options)
    }else{opaconf <- opafit}
    opaconf <- modifyList(opaconf,list(...))
    if(opaconf$is.plot.cb){
	
## Entry of t0 here is mearly an artifact from abrem code	
## Neither abrem nor wblr permit calculation of conf bounds on 3p models
## Validity of this would be questionable because uncertainty in t0 is not accounted for
        t0 <- 0
        

        if(!is.null(blc$bounds$Datum))
            lines(y=F0inv(blc$bounds$unrel,opaconf$log),
                x=blc$bounds$Datum-t0,
                col=opaconf$col,lwd=1,lty=2)
        if(!is.null(blc$bounds$Lower))
            lines(y=F0inv(blc$bounds$unrel,opaconf$log),
                x=blc$bounds$Lower-t0,col=opaconf$col,
                lwd=opaconf$lwd,lty=opaconf$lty)
        if(!is.null(blc$bounds$Upper))
            lines(y=F0inv(blc$bounds$unrel,opaconf$log),
                x=blc$bounds$Upper-t0,col=opaconf$col,
                lwd=opaconf$lwd,lty=opaconf$lty)
    }
}

plotConfsInFit <- function(fit,opadata,...){
    arg <- list(...)
    if(!is.null(fit$conf$blives)){
        if(!is.null(fit$options)){
            opafit <- modifyList(opadata,fit$options)
        }else{opafit <- opadata}
        lapply(fit$conf$blives,plotSingleConfBound,opafit=opafit,...)
    }
#    else{if(arg$v >= 1)message(match.call()[[1]],
#        ": This fit contains no confidence calculations for B-lives.")}
}
plotSingleFit <- function(fit,opadata,...){
    opafit <- opadata
    if(!is.null(fit$options)){
        opafit <- modifyList(opadata,fit$options)
    }
    opafit <- modifyList(opafit,list(...))
    if(opafit$is.plot.fit){
## removing $threshold influence	


		tz<-0
		if(!is.null(fit$t0) && !fit$modified) {
			tz<-fit$t0
		}
		
        if(!is.null(fit$beta) && !is.null(fit$eta)){
# One routine suits all . . .

                cret <- curve(F0inv(pweibull(x-tz,
                    fit$beta,fit$eta),opafit$log),
                    add=TRUE,n=1001,
                        # n=1001 is needed for displaying the extreme
                        # curvature towards -Inf with low Beta values
                        # like 0.1
                    col=opafit$col,lwd=opafit$lwd,lty=opafit$lty,
                    xlim=getPlotRangeX(opafit$log),
                    log=opafit$log)
                cret$y[is.infinite(cret$y)] <- NA
                    # works for weibull canvas
                cret$y[cret$y==0] <- NA
                    # replacing zero's is needed for lognormal canvas.
                imin <- which.min(cret$y)
                lines(rep(cret$x[imin],2),
                    y=c(cret$y[imin],getPlotRangeY(opafit$log)[1]),
                    col=opafit$col,lwd=opafit$lwd,lty=opafit$lty)
                    # plot vertical line towards -Inf
#            }
        }

        if(!is.null(fit$meanlog) && !is.null(fit$sdlog)){
            ### lognormal ###
##            if(opafit$verbosity >= 1)message(
##                "plotSingleFit: Adding Lognormal fit ...")
            curve(F0inv(plnorm(x-tz,fit$meanlog,fit$sdlog),opafit$log),
                add=TRUE,
                col=opafit$col,lwd=opafit$lwd,lty=opafit$lty,
                xlim=getPlotRangeX(opafit$log),
                log=opafit$log)
                # TODO: deal with Inf and -Inf values in the curve argument zo that the curce always extends to the edges of the plotting regions
        }
        if(!is.null(fit$rate)){
            ### exponential ###
##            if(opafit$verbosity >= 1)message(
##                "plotSingleFit: Adding Exponential fit ...")
            curve(F0inv(pexp(x+t0,fit$rate),opafit$log),add=TRUE,
                col=opafit$col,lwd=opafit$lwd,lty=opafit$lty,
                xlim=getPlotRangeX(opafit$log),
                log=opafit$log)
        }
    }
    invisible()
}

getPlotRangeX <- function(log){
    if(log %in% c("x","xy","yx")) 10^par("usr")[1:2]
    else par("usr")[1:2]
}

getPlotRangeY <- function(log){
    if(log %in% c("y","xy","yx")) 10^par("usr")[3:4]
    else par("usr")[3:4]
#    if(log %in% c("y","xy","yx"))) 10^par("usr")[1:2]
#    else par("usr")[1:2]
}
