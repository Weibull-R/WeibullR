## wblr.fit.R
## refactored from abrem.fit.R originally authored by Jurgen Symynck, April 2014
## Copyright 2014-2017 OpenReliability.org
#
# For more info, visit http://www.openreliability.org/
#
# For the latest version of this file, check the Subversion repository at
# http://r-forge.r-project.org/projects/Weibull-R/
#
#-------------------------------------------------------------------------------
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
#

wblr.fit <- function(x, modify.by.t0=FALSE,...){
    # x is a single wblr
	if(!class(x)=="wblr"){
		stop("\"x\" argument is not of class \"wblr\".")
	}


## okay lets do some realistic validation 
## It will probably be better to splitfitargs(...) for more comprehensive validations
	arg <- list(...)
	if(!is.null(c(arg$log,arg$canvas))) stop("cannot set log or canvas option in wblr.fit")
	
	if(!is.null(x$data$dlines)) {
	if(!is.null(arg$method.fit)) {
		if("rr" %in% tolower(arg$method.fit)){
## The user tried to set the method fit to rank regression with intervals in data
## so hard stop here, should be a very rare occurence
			stop("rank regression is not performed on interval data, use an mle fit")
		}
	}
	}
## silently correct a possible R user error	
	if(!is.null(arg$dist)) {
	if(arg$dist=="lnorm") arg$dist<-"lognormal"
	}

## actually it is not necessary to have a call to such an enclosed function now
## it is shown here for original code reference, only	
##		calculateSingleFit() 
#                        x <- lapply(x,calculateSingleFit,...)



##calculateSingleFit <- function(){
    # x is still the original single wblr object
	opadata <- x$options
	opafit <- modifyList(opadata,arg)
	
	if(!is.null(x$data$dlines)) {
		if("rr" %in% tolower(opafit$method.fit)){
## change from the default should have been called for in the (...), just warn and move on
			warning("rank regression is not performed on interval data method.fit has been set to mle")
			opafit$method.fit<-"mle"
		}
	}
		
## here are Jurgen's original validations, corrected for opafit
    supported_dist <- c(
        "weibull","weibull2p","weibull3p",
        "lognormal","lognormal2p","lognormal3p")
    supported_fit <-  c("rr","mle","mle-rba","mle-unbias")
	if(is.null(opafit$dist)){
			opafit$dist <- "weibull2p"
	}
	if(!any(tolower(opafit$dist) %in% supported_dist)){
		stop(paste0(opafit$dist," is not a supported fit distribution."))
	}
	if(!any(tolower(opafit$method.fit) %in% supported_fit)){
		stop(paste0(opafit$dist," is not a supported fit method."))
	}	
	
		
	if(modify.by.t0==TRUE) {
		if(tolower(opafit$dist) %in% c("weibull3p", "lognormal3p")){
# wipe the fit slate clean for this new wblr object
			x$fit<-NULL
		}else{
			modify.by.t0<-FALSE
			warning("modify.by.t0 ignored for non-3p fitting")
		}
	
	}
	########################
    #  main function body  #
    ########################

    atleastonefit <- FALSE
    if(is.null(x$fit)){
        ## Creating the first fit in the wblr object...")
        i <- 1
        x$fit <- list()
    }else{
        ## Appending a new fit to the existing wblr object...")
        i <- length(x$fit)+1
    }
    x$fit[[i]] <- list()
	
    op <- unique(c(names(x$options),names(opafit)))
        # this is needed to add options from opafit into li that
        # are NULL in x$options
        # TODO:tolower() needed?
    if(length(li <- opafit[sapply(op,function(y){
        !identical(x$options[[y]], opafit[[y]])})]) > 0){
        x$fit[[i]]$options <- li
        # the above enlists only options that are different from the wblrs
        # 'main' options. This excludes options$dist and options$method.fit
    }
	
	## the code above does not fill x$fit$options$dist or x$fit$options$method.fit
	## it turns out that these were filled within the regression sections of original calculateSingleFit
	## I will force this now, without these a legend will not be produced.
	if(is.null(x$fit[[i]]$options)) {
	x$fit[[i]]$options<-list()
	}
	
	x$fit[[i]]$options$dist<-opafit$dist
## need to assure interval data is not attempted by rank regression
	if(!is.null(x$data$dlines) && any(c("rr","rr2") %in% tolower(opafit$method.fit))){
		modifyList(opafit, list(method.fit="mle"))
	}
	x$fit[[i]]$options$method.fit<-opafit$method.fit
		

	
## why the heck this duplication?	
##(These are most likely used for the Legend which might be generated for each fit)
    x$fit[[i]]$n    <- x$n
    x$fit[[i]]$fail <- x$fail
    x$fit[[i]]$cens <- x$cens
	x$fit[[i]]$discovery <- x$discovery
	x$fit[[i]]$interval  <- x$interval

## It turns out that this code is general to all fitting methods:

	if(tolower(opafit$dist) %in% c("weibull","weibull2p","weibull3p")){	
		fit_dist<-"weibull"
	}else{

		if(tolower(opafit$dist) %in% c("lnorm","lognormal","lognormal2p", "lognormal3p")){
			fit_dist<-"lnorm"
		}else{
		## Note: lslr contains experimental support for "gumbel"
			stop(paste0("dist option ", opafit$dist, "is not recognized for distribution fitting"))
		}
	}
	
	npar<-2

	if(tolower(opafit$dist) %in% c("weibull3p", "lognormal3p")){
		npar<-3
	}	
	
    if(any(c("rr","rr2") %in% tolower(opafit$method.fit))){
        #  ____             _                                       _
        # |  _ \ __ _ _ __ | | __  _ __ ___  __ _ _ __ ___  ___ ___(_) ___  _ __
        # | |_) / _` | '_ \| |/ / | '__/ _ \/ _` | '__/ _ \/ __/ __| |/ _ \| '_ \
        # |  _ < (_| | | | |   <  | | |  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
        # |_| \_\__,_|_| |_|_|\_\ |_|  \___|\__, |_|  \___||___/___/_|\___/|_| |_|
        #                                   |___/

		regression_order<-"XonY" #default for anything other than "yonx" that might be at method.fit[2]
		if(length(opafit$method.fit)>1){
		if(tolower(opafit$method.fit[2])=="yonx") {
			regression_order<-"YonX"
		}
		}

		
		## content of fit_vec result depends on fit_dist, and npar
		try(fit_vec<-lslr(x$data$dpoints, dist=fit_dist, npar=npar, reg_method=regression_order))

		if(!is.null(fit_vec)){
## Jurgen's use of super assignment (<<-) appears to have been an error here		
			atleastonefit<-TRUE
			if(fit_dist=="weibull") {
				x$fit[[i]]$beta <- fit_vec[2]
				x$fit[[i]]$eta <- fit_vec[1]				
				}
			if(fit_dist=="lnorm") {
				x$fit[[i]]$meanlog<-fit_vec[1]
				x$fit[[i]]$sdlog<-fit_vec[2]	
			}
## gof is common to both fit_dist types			
			x$fit[[i]]$gof <- list()
			if(npar==2){
				x$fit[[i]]$gof$r2 <- fit_vec[[3]]
				x$fit[[i]]$gof$prr <- fit_vec[[4]]
			}else{  ## only other option is 3p
				x$fit[[i]]$t0 <- fit_vec[3]
				x$fit[[i]]$gof$r2 <- fit_vec[[4]]
## removing effects of threshold option				
##				if(!is.null(opafit$threshold)){
##					if(is.logical(opafit$threshold) && opafit$threshold)
##						 x$options$threshold <- fit_vec[[3]]
##				}
					# this overwrites any threshold setting at the data level with a number
					# TODO: Jurgen notes this is not the way to go when trying to implement support for
					# threshold with plot.wblr()
			}
		}else{
			##vm(0,"calculateSingleFit: Fitting failed.")
			x$fit[i] <<- list(NULL)
				# note that is.null(x$fit[[i]]) will exit with an error
				# TODO: replace with x$fit[[i]] <<- list(NULL)

		}		
	} ## close rank regression block
	
	

    if(any(c("mle","mle-rba","mle-unbias") %in% tolower(opafit$method.fit))){
        #  __  __ _     _____
        # |  \/  | |   | ____|
        # | |\/| | |   |  _|
        # | |  | | |___| |___
        # |_|  |_|_____|_____|	
		
		
## prepare arguments for mlefit		
## usage: mlefit(x, dist="weibull", npar=2, debias=NULL, optcontrol=NULL)
## x is the x$data$lrq_frame
## dist is fit_dist
## npar is npar
## no need for optcontrol changes from default

## interpret input for debias
		debias<-NULL
		if(tolower(opafit$method.fit) == "mle-rba")  debias <- "rba" 
		if(tolower(opafit$method.fit) == "mle-unbias") { 
			if(fit_dist == "weibull") {
				debias <- "hirose-ross" 
			}else{
#mle-unbias taken as mle-rba for lognormal
				debias <- "rba"
			}
		}
		
		try(fit_vec<-mlefit(x$data$lrq_frame, fit_dist, npar, debias))
		
        if(!is.null(fit_vec)){
			atleastonefit<-TRUE
			if(fit_dist=="weibull") {
				x$fit[[i]]$beta <- fit_vec[2]
				x$fit[[i]]$eta <- fit_vec[1]				
			}
			if(fit_dist=="lnorm") {
				x$fit[[i]]$meanlog<-fit_vec[1]
				x$fit[[i]]$sdlog<-fit_vec[2]				
			}
		
		   x$fit[[i]]$gof <- list()
			if(npar==2){
				x$fit[[i]]$gof$loglik <- fit_vec[[3]]

			}else{  ## only other option is 3p
				x$fit[[i]]$t0 <- fit_vec[3]
				x$fit[[i]]$gof$loglik <- fit_vec[[4]]
## Removing effects of threshold option
##				if(!is.null(opafit$threshold)){
##					if(is.logical(opafit$threshold) && opafit$threshold)
##						 x$options$threshold <- fit_vec[[3]]
##				}
					# this overwrites any threshold setting at the data level with a number
					# TODO: Jurgen notes this is not the way to go when trying to implement support for
					# threshold with plot.wblr()
			}		   

		}
				
	} ## close mle block
		
    if(!atleastonefit){
        warning("*** calculateSingleFit: Nothing has been fitted.  ***\n",
                '*** Does \"method.fit\" include sensible options?   ***')
        # x$fit[[i]] <- NULL
    }
##    if(is.numeric(opafit$threshold))
##        x$options$threshold <- opafit$threshold
        # overwrite any previously set data-level-t0 to the one specified as an argument to wblr.fit()
        # Don't know why - here - you MUST use <- in favor of <<- ...
		## correct - the x object AND the opafit object are both in current environment scope
		## indeed x$options$threshold is unknown to the global environment on first iteration here.
		
## calculateSingleFit has no return object as it modifies x$fit[[i]] directly
 
	#x$data$modified<-FALSE  ## turns out there is no use for this as notation in the fit is enough
	x$fit[[i]]$modified<-FALSE
	 if(modify.by.t0==TRUE) {
## subtract all input data with x$fit[[i]]$t0
		if(!is.null(x$fit[[i]]$t0) ){
			if(!x$fit[[i]]$t0 < min(x$data$lrq_frame$left)) {
				stop("t0 too large for data modification")
			}
			x$data$lrq_frame$left<-x$data$lrq_frame$left - x$fit[[i]]$t0
			if(!is.null(x$data$dpoints$time)) {
			x$data$dpoints$time<-x$data$dpoints$time - x$fit[[i]]$t0
			}
			if(!is.null(x$data$dlines$t1)) {
			x$data$dlines$t1<-x$data$dlines$t1e - x$fit[[i]]$t0
			}

## needed to control legend and fitted curve			
			x$fit[[i]]$modified<-TRUE
			## No, do not remove the "3p" suffix from the fit$dist to permit conf calculations
			## This is still marked as a 3p fit, so no conf calculations can be made.
			#SL<-nchar(x$fit[[i]]$options$dist)
			#x$fit[[i]]$options$dist<-substr(x$fit[[i]]$options$dist,1,SL-2)		
		}else{
			 warning("t0 not found for data modification")
		}



	}	  
    
x

}