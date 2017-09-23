## wblr.conf.R
## inspired by code originally authored by Jurgen Symynck, April 2014
# Extensive re-write by David J. Silkworth simplifying calls to independent, generalized
# functions to generate confidence interval bounds from which Blife estimates can also be determined.
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
wblr.conf <- function(x,...){
    # x is a single wblr object
	if(!(class(x)=="wblr")){
		stop('Argument \"x\" is not of class \"wblr\" ')
	}
## using findMaxDataRange from plot.wblr, which takes a list of wblr objects
## so simply convert this x to a single item list
	#dr <- findMaxDataRange(list(x))
	datarange<- findMaxDataRange(list(x))

	if(!is.null(x$fit)){	
## usage: calculateSingleConf(fit,xdata,opadata,datarange,...)	
## only acting on the last fit added to the object
	fit<-x$fit[[length(x$fit)]]
	xdata<-x$data
	opadata<-x$options
	
	#	x$fit[[length(x$fit)]]<- calculateSingleConf(
	#		x$fit[[length(x$fit)]],
	#		x$data, opadata=x$options,datarange=dr,...
	#	)


## calculateSingleConf is now inserted here
    # fit is a single fit
    arg <- list(...)

	if(is.null(fit$options$dist)) {
		stop("missing a fit distribution")
	}else{ 
	##perhaps test for valid fit distribution?? (How could this be wrong at this point?) 
	}
		
	if(!is.null(fit$options))  {
		opafit <- modifyList(opadata,fit$options)
	}
	opaconf <- modifyList(opafit,arg)

## prepare the descriptive quantiles  -  0.5 and F0(0) are pivot points for pivotal corrections			
	unrel <- c(F0(seq(F0inv(1e-3), F0inv(0.999),length.out=25)),
				opaconf$unrel, 0.5, F0(0))

	unrel <- unique(signif(unrel[order(unrel)]))
		# signif() has been used to eliminate
		# any identical looking unreliability
		# levels that differ only at place far
		# from the decimal point
		
## prepare the list objects
	if(is.null(fit$conf)){
		i <- 1
		fit$conf <- list()
	}else{
## Appending a new confidence calculation to the fit	
		i<-length(fit$conf)+1
		fit$conf[[i]]<- list()
	}
							
## I am not convinced that the second layer list fit$conf$blives is needed
## certainly for now, plot.wblr is looking for this though.							
							
## removing fit$conf$blives
#	if(is.null(fit$conf$blives)){
##Creating the first B-life confidence calculation in the fit
#		i <- 1
#		fit$conf$blives <- list()
#	}else{
## Appending a new B-life confidence calculation to the fit
#		i <- length(fit$conf$blives)+1
#	}
#	fit$conf$blives[[i]] <- list()
	
	atLeastOneBLifeConf <- FALSE
                            
                            
                            

	if("bbb" %in% tolower(opaconf$method.conf.blives)){
		#  ____  ____  ____
		# | __ )| __ )| __ )
		# |  _ \|  _ \|  _ \
		# | |_) | |_) | |_) |
		# |____/|____/|____/
		
## This is a modified version of the BBB confidence interval calculation from abrem.
## it warrants an independent, generalized function returning consistent output
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type <- "bbb"
		fit$conf[[i]]$cl <- opaconf$cl
		fit$conf[[i]]$sides <- opaconf$conf.blives.sides
## This is why data was redundantly added to the fit, but we have exposed the x$data list
##  as an xdata argument providing access to xdata$lrq_frame, xdata$dpoints and xdata$dlines

## soon it will be necessary to combine the xdata$dlines (modified to conform)
		sx<-xdata$dpoints
		
##		Beta Binomial "Z" factors are non-parametric
		Zlo<-qbeta((1-opaconf$cl)/2,sx$adj_rank,fit$n-sx$adj_rank+1)
		Zhi<-qbeta(1-(1-opaconf$cl)/2,sx$adj_rank,fit$n-sx$adj_rank+1)
		
		if(tolower(fit$options$dist) %in% c("weibull","weibull2p")){
			#Lower= bbb(sx$adj_rank,fit$n,(1-opaconf$cl)/2,fit$beta,fit$eta)
			#Upper= bbb(sx$adj_rank,fit$n,1-(1-opaconf$cl)/2,fit$beta,fit$eta)
			Lower<- qweibull(Zlo,fit$beta,fit$eta)
			Upper<- qweibull(Zhi,fit$beta,fit$eta)
		}else{
			if(tolower(fit$options$dist) %in% c("lnorm","lognormal","lognormal2p")){
				##stop("lognormal bbb not implemented yet")
				## build da for lognormal
			Lower<- qlnorm(Zlo,fit$meanlog,fit$sdlog)
			Upper<- qlnorm(Zhi,fit$meanlog,fit$sdlog)
			}else{
				stop(paste0("distribution ", fit$options$dist, " not supported for bbb."))
			}
		
		}

## using log transform on ppp for better linearization for interpolation
#		lo <- approxfun(log(sx$ppp),log(da$Lower))
#		up <- approxfun(log(sx$ppp),log(da$Upper))
		lo <- approxfun(log(sx$ppp),log(Lower))
		up <- approxfun(log(sx$ppp),log(Upper))
		bl <- log(unrel)
		da <- data.frame(unrel=unrel,Lower=exp(lo(bl)),Upper=exp(up(bl)))
## Add in only the two extreame data points for graphing		
		da <- rbind(da,data.frame(unrel=min(sx$ppp), Lower=min(Lower), Upper=min(Upper)))
		da <- rbind(da,data.frame(unrel=max(sx$ppp), Lower=max(Lower), Upper=max(Upper)))
		da <- da[order(da$unrel),] 
		da <- da[!duplicated(da$unrel),]
		fit$conf[[i]]$bounds <- da
		
		op <- unique(c(names(opafit),names(opaconf)))
			# this is needed to add options from opafit into li that
			# are NULL in opafit
			# TODO:tolower() not needed?
		if(length(li <- opaconf[sapply(op,function(y){
			!identical(opafit[[y]], opaconf[[y]])})]) > 0){
			fit$conf[[i]]$options <- li
		}
		
  		atLeastOneBLifeConf <- TRUE									
	}
								
								
	if("lrb" %in% tolower(opaconf$method.conf.blives)){
		#  Likelihood Ratio Bounds to be added back in when appropriate

	}
								
								
	if(any(c("mcpivotals","mcpivotal") %in% tolower(opaconf$method.conf.blives))){
		#                       _            _        _
		#  _ __ ___   ___ _ __ (_)_   _____ | |_ __ _| |___
		# | '_ ` _ \ / __| '_ \| \ \ / / _ \| __/ _` | / __|
		# | | | | | | (__| |_) | |\ V / (_) | || (_| | \__ \
		# |_| |_| |_|\___| .__/|_| \_/ \___/ \__\__,_|_|___/
		# 
	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {
		stop("3-parameter fits are not handled by mcpivotals")
	}			
		if(!is.integer(opaconf$seed))  {
			##warning(paste0("opaconf$seed: ",opaconf$seed,"is not an integer"))
			opaconf$seed<-1234
		}
		
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type   <- "mcpivotals"
		fit$conf[[i]]$S      <- 10^4
		fit$conf[[i]]$seed   <- opaconf$seed
		fit$conf[[i]]$rgen   <- opaconf$rgen
		fit$conf[[i]]$cl     <- opaconf$cl
		fit$conf[[i]]$sides  <- opaconf$conf.blives.sides
		fit$conf[[i]]$unrel <- opaconf$unrel
		ret <- NULL

## now using pivotalMC from abremPivotals . . .									
## \usage{
##pivotalMC(x, event=NULL, dist="weibull", reg_method="XonY", R2, CI, unrel,
##  P1=1.0, P2=1.0, S=10^4,seed=1234, ProgRpt=FALSE)

## data with intervals should already have been excluded from pivotalMC
		if(!is.null(xdata$dlines)) stop("mcpivotals not performed with interval data")
 
## must form an event vector from xdata$lrq_frame$right, if value is >0 event=1, else event=0 
## what about qty field? Okay, lets run a quick brute force loop:

 
		event_vec<-NULL
		for(di in 1:dim(xdata$lrq_frame)[1]) {
			if(xdata$lrq_frame$right[di]>0) {
				this_event=rep(1,xdata$lrq_frame$qty[di])
			}else{
				this_event=rep(0,xdata$lrq_frame$qty[di])
			}
			event_vec<-c(event_vec,this_event)
		}

		if(tolower(fit$options$dist) %in% c("weibull","weibull2p")) {
			fit_dist<-"weibull"
## prepare the weibull parameters for pivotal evaluation per Jurgen Symynck's	method	
			P1<-lslr(getPPP(qweibull(xdata$dpoints$ppp,1,1)))[1]	
			P2<-lslr(getPPP(qweibull(xdata$dpoints$ppp,1,1)))[2]			
		}else{
			if(tolower(fit$options$dist) %in% c("lognormal","lognormal2p")) {
				fit_dist<-"lnorm"
## prepare the lognormal parameters for pivotal evaluation per 
## David Silkworth's adaptation of Jurgen Symynck's weibull method				
				completePPP<-getPPP(xdata$dpoints$time)	
				comp_prob<- as.vector(completePPP[,2])	
				P1=mean(qnorm(xdata$dpoints$ppp,0,1))	
				P2=sd(qnorm(xdata$dpoints$ppp,0,1))/sd(qnorm(comp_prob,0,1))
			## division by sd(qnorm(comp_prob,0,1))	
			## corrects for resolution of sd for complete failures to 1.0	
			## critically needed for lognormal pivotal parameter determination			
				
			}else{
				stop("fit distribution not weibull or lognormal for pivotalMC")
			}
		}

		
		regression_order<-"XonY" #default for anything other than "yonx" that might be at method.fit[2]
		if(length(opafit$method.fit)>1){
		if(tolower(opafit$method.fit[2])=="yonx") {
			regression_order<-"YonX"
		}
		}
		
## R2 and CI control format of return object 
## given R2 value 0 suppresses prr and ccc2 output
## given CI value >0 <1 a "pivotals" dataframe of interval bounding values will be returned
		
		

		
		ret<-pivotal.rr(
			xdata$dpoints, 
			event=event_vec,
			dist=fit_dist, 
			reg_method=regression_order,
			R2=0,
			CI=opaconf$cl,
			unrel=unrel,
			S=10^4,
			P1=P1,
			P2=P2,
			seed=opaconf$seed,
			ProgRpt=FALSE
		)
				

			

				if(!is.null(ret)){
					atLeastOneBLifeConf <- TRUE
					if(fit_dist=="weibull") {
## David Silkworth's final adaptation to nail the bounds around the fitted line
		median_slope<-(log(log(1/(1-unrel[length(unrel)])))-log(log(1/(1-unrel[1]))))/(ret[length(unrel),2]-ret[1,2])	
		med_pos<-which(unrel > F0(0)-10^-5)[1]
		median_intercept<-ret[med_pos,2]-log(log(1/(1-unrel[med_pos])))/median_slope
		ret<-(ret-median_intercept)*median_slope					
										
						fit$conf[[i]]$bounds <- cbind(unrel,
							exp(log(fit$eta)+ ret/fit$beta))
					}
					if(fit_dist=="lnorm") {
## David Silkworth's final adaptation to nail the bounds around the fitted line
		## check the slope of the median pivotals to get correction to 1.0	
		median_slope<-(qnorm(unrel[length(unrel)], 0, 1)  - qnorm(unrel[1], 0, 1))/(ret[length(unrel),2]-ret[1,2])
		## don't know why I originally used 60% quantile in MRRln2p (same as MRRw2p)
		med_pos<-which(unrel > .5-10^-5)[1]		
		median_intercept<-ret[med_pos,2]-qnorm(unrel[med_pos], 0, 1)/median_slope	
		ret<-(ret-median_intercept)*median_slope

		## interpret the pivotals for the log plot	
		#plot_piv<-(adj_piv)*fit[2]+fit[1]		
## some confusion as to which way to rotate
#                                    exp(fit$meanlog + ret/fit$sdlog) ## just wrong
##                                    exp(fit$meanlog - ret/fit$sdlog)) ## appears to be flipped on y-axis at 50% intercept
						fit$conf[[i]]$bounds <- cbind(unrel,
							exp(ret*fit$sdlog + fit$meanlog))
					}
					
					names(fit$conf[[i]]$bounds) <- c("unrel","Lower","Datum", "Upper")
					op <- unique(c(names(opafit),names(opaconf)))
						# this is needed to add options from opafit into li that
						# are NULL in opafit
						# TODO:tolower() not needed?
					if(length(li <- opaconf[sapply(op,function(y){
						!identical(opafit[[y]], opaconf[[y]])})]) > 0){
						fit$conf[[i]]$options <- li
					}
				}else{
					message("calculateSingleConf: Confidence calculation failed.")
					fit$conf[[i]] <- NULL
				}
										
										
	}  ## end mcpivotals


	if("fmbounds" %in% tolower(opaconf$method.conf.blives)){
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type   <- "fmbounds"
		fit$conf[[i]]$cl     <- opaconf$cl
		fit$conf[[i]]$sides  <- opaconf$conf.blives.sides
		fit$conf[[i]]$unrel <- opaconf$unrel
		ret <- NULL
	
	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {
		stop("3-parameter fits are not handled by FMbounds")
	}
## assure valid dist names for FMbounds and the debias functions	
	if(tolower(fit$options$dist) %in% c("weibull","weibull2p")){
		fit_dist<-"weibull"
	}
	if(tolower(fit$options$dist) %in% c("lnorm","lognormal","lognormal2p")){
		fit_dist<-"lognormal"
	}			
		
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
if(!is.null(debias)) fit$conf[[i]]$debias <- debias	

## usage FMbounds(x, dist="weibull", CI=.95, unrel=NULL, debias=NULL, show=FALSE)		
		ret<-FMbounds(xdata$lrq_frame, dist=fit$options$dist, CI=opaconf$cl, unrel=unrel, debias=debias)

		if(!is.null(ret)){
			atLeastOneBLifeConf <- TRUE

			ret$percentile<-ret$percentile/100
			
			fit$conf[[i]]$bounds <- ret
			names(fit$conf[[i]]$bounds) <- c("unrel","Lower","Datum", "Upper")	
								
			op <- unique(c(names(opafit),names(opaconf)))
				# this is needed to add options from opafit into li that
				# are NULL in opafit
				# TODO:tolower() not needed?
			if(length(li <- opaconf[sapply(op,function(y){
				!identical(opafit[[y]], opaconf[[y]])})]) > 0){
				fit$conf[[i]]$options <- li
			}
		}else{
			message("calculateSingleConf: Confidence calculation failed.")
			fit$conf[[i]] <- NULL
		}
	
	}  # end FM bounds

	

					
		if(!atLeastOneBLifeConf)  {
				warning("calculateSingleConf: The fit argument is empty or contains no fits.")
		}
	
	x$fit[[length(x$fit)]]<-fit
	}	
	x	
		
	

}	


