## wblr.conf.R
## inspired by code originally authored by Jurgen Symynck, April 2014
# Extensive re-write by David J. Silkworth simplifying calls to independent, generalized
# functions to generate confidence interval bounds from which Blife estimates can also be determined.
## Copyright 2014-2021 OpenReliability.org
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
	if(!is(x, "wblr")) stop('Argument \"x\" is not of class \"wblr\" ')

## using findMaxDataRange from plot.wblr, which takes a list of wblr objects
## so simply convert this x to a single item list
	#dr <- findMaxDataRange(list(x))
	datarange<- findMaxDataRange(list(x))

	if(!is.null(x$fit)){
## only acting on the last fit added to the object
	fit<-x$fit[[length(x$fit)]]
	xdata<-x$data
	opadata<-x$options


## validation of the ... arguments
    arg <- list(...)
	if(!is.null(arg$method.conf.blives)) {
		warning("method.conf.blives has been depreciated in favor of method.conf")
		arg$method.conf<-arg$method.conf.blives
		arg<-modifyList(arg, list(method.conf.blives=NULL))
	}
	if(!is.null(arg$unrel.n)) modifyList(arg, list(num_dq=arg$unrel.n))
	## silent refactoring for undocumented abrem list option.
	##
	## tests for valid confidence calculations in args should be made here
	if(!is.null(c(arg$log,arg$canvas))) stop("cannot set log or canvas option in wblr.conf")
	if(!is.null(arg$dist)) stop("cannot set the fit distribution in wblr.conf")
	if(!is.null(arg$method.fit)) stop("cannot set the fit method in wblr.conf")

if(!is.null(arg$method.conf)) {
	if(substr(tolower(arg$method.conf),1,2) == "mc") {
		warning("'mcpivotal' has been depreciated for method.conf, use 'pivotal-rr'")
		arg$method.conf<-"pivotal-rr"
	}
}

#	if(!is.null(fit$options))  {
## it should be okay if fit$options is null anyway
		opafit <- modifyList(opadata,fit$options)
#	}
	opaconf <- modifyList(opafit,arg)


DescriptivePercentiles<-function(dplabel)  {		
## Descriptive percentiles are the percentile positions at which points on the curved bounds are		
## calculated, so that a smoothed curve can be plotted by linear interpolation.		
## This function will return a vector of percentiles based on a named set that can be stored		
##  in the options.wblr list. This is helpful for comparison with other software.		
## It is expected that this function will most often be called by its aleas, DP.		
		
	if(tolower(dplabel)=="minitab") {	
	## these descriptive percentiles match Minitab unchangeable defaults (27 values)	
		dp<-c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))
	}	
		
	if(tolower(dplabel)=="supersmith") {	
		## descriptive percentiles for comparison with SuperSMITH (limit of 15 values)
		dp<-c(.01, .02, .05, .10, .15, .20, .30, .40, .50,  .60, .70, .80, .90, .95, .99)
	}	
		
	if(tolower(dplabel)=="user") {	
		dp<-opaconf$user_dp
	}	
		
	if(tolower(dplabel)=="abrem")  {	
	## this is the original default by Jurgen Symynck for package abrem	
	## it produces evenly spaced points across the y limits of a weibull canvas	
	#F0(seq(F0inv(1e-3), F0inv(0.999),length.out=??)) Attempting to hold a constant number of points.	
	spec_pts<-c(opaconf$blife.pts, 0.5, 1-exp(-exp(0)))	
	len_out<-opaconf$num_dp-length(unique(spec_pts)) # in case any blife.points duplicate pivot points	
	mini <- min(c(opaconf$ylim[1]/10,datarange$yrange[1]/10),0.001)	
	maxi <- max(c((1-(1-opaconf$ylim[2])/10),(1-(1-datarange$yrange[2])/10),0.999))	
	dp<-1-exp(-exp(seq(log(qweibull((mini),1,1)), log(qweibull((maxi),1,1)),length.out=len_out)))	
	}	
		
	dp	
}		
		
DP<-DescriptivePercentiles		
		
## prepare the descriptive percentiles -  0.5 and F0(0) are pivot points for pivotal corrections		
	mini <- min(c(opaconf$ylim[1]/10,datarange$yrange[1]/10),0.001)	
	maxi <- max(c((1-(1-opaconf$ylim[2])/10),(1-(1-datarange$yrange[2])/10),0.999))	
		
	lodp<-1-exp(-exp(log(qweibull((mini),1,1))))	
	hidp<-1-exp(-exp(log(qweibull((maxi),1,1))))	
	unrel <- c(DP(opaconf$dp),opaconf$blife.pts, 0.5, 1-exp(-exp(0)),lodp,hidp)	
		
	unrel <- unique(signif(unrel[order(unrel)]))	
		# signif() has been used to eliminate any identical looking descriptive
		# percentiles that differ only at place far from the decimal point

## prepare the list objects
	if(is.null(fit$conf)){
		i <- 1
		fit$conf <- list()
	}else{
## Appending a new confidence calculation to the fit
		i<-length(fit$conf)+1
		fit$conf[[i]]<- list()
	}


	atLeastOneBLifeConf <- FALSE




	if("bbb" %in% tolower(opaconf$method.conf)){
		#  ____  ____  ____
		# | __ )| __ )| __ )
		# |  _ \|  _ \|  _ \
		# | |_) | |_) | |_) |
		# |____/|____/|____/

## This is a modified version of the BBB confidence interval calculation from abrem.
## it warrants an independent, generalized function returning consistent output
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type <- "bbb"
		fit$conf[[i]]$ci <- opaconf$ci
## all interval bounds are double sided having a Confidence Interval CI
## B-lives are reported with single side Confidence Level
##		fit$conf[[i]]$sides <- opaconf$conf.blives.sides
## This is why data was redundantly added to the fit, but we have exposed the x$data list
##  as an xdata argument providing access to xdata$lrq_frame, xdata$dpoints and xdata$dlines

## Need to combine ppp and adjusted ranks for points and lines.
		sx<-NULL
		if(!is.null(xdata$dpoints)) {
		sx<-xdata$dpoints[,2:4]
		}
		if(!is.null(xdata$dlines)) {
		sx<-rbind(sx,xdata$dlines[,3:5])
		sx<-sx[order(sx$adj_rank),]
		}

##		Beta Binomial "Z" factors are non-parametric
		Zlo<-qbeta((1-opaconf$ci)/2,sx$adj_rank,fit$n-sx$adj_rank+1)
		Zhi<-qbeta(1-(1-opaconf$ci)/2,sx$adj_rank,fit$n-sx$adj_rank+1)

		if(tolower(fit$options$dist) %in% c("weibull","weibull2p")){
			Lower<- qweibull(Zlo,fit$beta,fit$eta)
			Upper<- qweibull(Zhi,fit$beta,fit$eta)
		}else{
			if(tolower(fit$options$dist) %in% c("lnorm","lognormal","lognormal2p")){
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
## Add in only the two extreme data points for graphing
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



	if(substr(tolower(opaconf$method.conf),1,7)=="pivotal"){								
		#                       _            _        _							
		#  _ __ ___   ___ _ __ (_)_   _____ | |_ __ _| |___							
		# | '_ ` _ \ / __| '_ \| \ \ / / _ \| __/ _` | / __|							
		# | | | | | | (__| |_) | |\ V / (_) | || (_| | \__ \							
		# |_| |_| |_|\___| .__/|_| \_/ \___/ \__\__,_|_|___/							
		#							
	npar=2								
	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {								
		npar=3							
	}								
	if(substr(tolower(fit$options$method.fit),1,2)!= "rr") {								
		stop("pivotal bounds are only applied on rank regression fits")							
	}								
									
		if(!is.integer(opaconf$seed))  {							
			##warning(paste0("opaconf$seed: ",opaconf$seed,"is not an integer"))						
			opaconf$seed<-1234						
		}							
									
		fit$conf[[i]]        <- list()							
		fit$conf[[i]]$type   <- "pivotal-rr"							
		fit$conf[[i]]$S      <- 10^4							
		fit$conf[[i]]$seed   <- opaconf$seed							
		fit$conf[[i]]$ci     <- opaconf$ci							
		fit$conf[[i]]$blife.pts <- opaconf$blife.pts							
		ret <- NULL							
									
## data with intervals should already have been excluded from pivotal bounds									
		if(!is.null(xdata$dlines)) stop("pivotal bounds not performed with interval data")							
									
		if(tolower(fit$options$dist) %in% c("weibull","weibull2p", "weibull3p")) {							
			fit_dist<-"weibull"						
## prepare the weibull parameters for pivotal evaluation per Jurgen Symynck's method									
			piv_fit<-lslr(getPPP(qweibull(xdata$dpoints$ppp,1,1)), abpval=FALSE)						
		}else{							
			if(tolower(fit$options$dist) %in% c("lnorm", "lognormal","lognormal2p","lognormal3p")) {						
				fit_dist<-"lnorm"					
## prepare the lognormal parameters for pivotal evaluation per adaptation of Jurgen Symynck's weibull method									
				piv_fit<-lslr(getPPP(qlnorm(xdata$dpoints$ppp,0,1)), abpval=FALSE)					
			}else{						
				stop("fit distribution not weibull or lognormal for pivotal.rr")					
			}						
		}							
		P1<-piv_fit[1]							
		P2<-piv_fit[2]							
									
		regression_order<-"XonY" #default for anything other than "yonx" that might be an rr method.fit							
		if(nchar(opafit$method.fit)>2){							
		if(substr(tolower(opafit$method.fit),4,7)=="yonx") {							
			regression_order<-"YonX"						
		}							
		}							
									
## R2 and CI control format of return object									
## given R2 value 0 suppresses prr and ccc2 output									
## given CI value >0 <1 a "pivotals" dataframe of interval bounding values will be returned									
									
									
##pivotal.rr(x, event=NULL, dist="weibull", npar=2, reg_method="XonY", R2, CI, 									
##  unrel, P1=1.0, P2=1.0, S=10^4,seed=1234, ProgRpt=FALSE)									
		ret<-pivotal.rr(							
			xdata$dpoints,						
			dist=fit_dist,						
			npar = npar,						
			reg_method=regression_order,						
			R2=0,						
			CI=opaconf$ci,						
			unrel=unrel,						
			S=opaconf$S,						
			P1=P1,						
			P2=P2,						
			seed=opaconf$seed,						
			ProgRpt=FALSE						
		)							
									
									
									
									
		if(!is.null(ret)){							
			atLeastOneBLifeConf <- TRUE						
			if(fit_dist=="weibull") {						
## David Silkworth's final adaptation to nail the bounds around the fitted line									
## F0(0) was coded as (1-exp(-exp(0)))									
				rotation_pos<-which(unrel > (1-exp(-exp(0)))-10^-5)[1]					
				rise<-(log(log(1/(1-unrel[length(unrel)])))-log(log(1/(1-unrel[rotation_pos]))))					
				run<-(ret[length(unrel),2]-ret[rotation_pos,2])					
				rotation_slope<-rise/run					
				rotation_intercept<-ret[rotation_pos,2]-log(log(1/(1-unrel[rotation_pos])))/rotation_slope					
				ret<-(ret-rotation_intercept)*rotation_slope					
## 3p bounds translate and rotate according to the 2p fit line			
				if(npar == 3) {		
					fit2p<-lslr(xdata$dpoints, dist=fit_dist,npar=2, abpval=FALSE)	
					fit$conf[[i]]$bounds <- cbind(unrel,	
						exp(log(fit2p[1])+ ret/fit2p[2]))
				}else{		
					fit$conf[[i]]$bounds <- cbind(unrel,	
						exp(log(fit$eta)+ ret/fit$beta))
				}		
			}						
									
			if(fit_dist=="lnorm") {						
## David Silkworth's final adaptation to nail the bounds around the fitted line									
				## don't know why I originally used 60% quantile in MRRln2p (same as MRRw2p)					
				rotation_pos<-which(unrel > .5-10^-5)[1]					
				## check the slope of the rotation pivotals to get correction to 1.0					
				rise<-(qnorm(unrel[length(unrel)], 0, 1)  - qnorm(unrel[rotation_pos], 0, 1))					
				run<-(ret[length(unrel),2]-ret[rotation_pos,2])					
				rotation_slope<-rise/run					
				rotation_intercept<-ret[rotation_pos,2]-qnorm(unrel[rotation_pos], 0, 1)/rotation_slope					
				ret<-(ret-rotation_intercept)*rotation_slope					
									
									
## some confusion as to which way to rotate									
#                                    exp(fit$meanlog + ret/fit$sdlog) ## just wrong									
##                                    exp(fit$meanlog - ret/fit$sdlog)) ## appears to be flipped on y-axis at 50% intercept									
									
## 3p bounds translate and rotate according to the 2p fit line			
				if(npar == 3) {										
					fit2p<-lslr(xdata$dpoints, dist=fit_dist,npar=2, abpval=FALSE)	
					fit$conf[[i]]$bounds <- cbind(unrel,	
						exp(ret*fit2p[2]+ fit2p[1]))									
				}else{					
					fit$conf[[i]]$bounds <- cbind(unrel,					
						exp(ret*fit$sdlog + fit$meanlog))
				}
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
									
									
	}  ## end pivotal-rr


#############################################################################
############   FM bounds ####################################################
#############################################################################

	if(any(c("fm","fmbounds") %in% tolower(opaconf$method.conf))) {
	if(substr(tolower(fit$options$method.fit),1,3)!= "mle") {								
		stop("fm bounds are only applicable on mle fits")							
	}	
	
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type   <- "fm"
		fit$conf[[i]]$ci     <- opaconf$ci
##		fit$conf[[i]]$sides  <- opaconf$conf.blives.sides
		fit$conf[[i]]$blife.pts <- opaconf$blife.pts
		ret <- NULL
	npar=2								
	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {								
		npar=3							
	}
##	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {
##		stop("confidence bounds are not prepared on 3-parameter fits")
##	}
## assure valid dist names for FMbounds and the debias functions
	if(tolower(fit$options$dist) %in% c("weibull","weibull2p","weibull3p")){
		fit_dist<-"weibull"
	}else{
		if(tolower(fit$options$dist) %in% c("lnorm","lognormal","lognormal2p","lognormal3p")){
			fit_dist<-"lognormal"
		}else{
			stop(paste0("dist ",fit$options$dist, " not recognized"))
		}
	}
		debias<-"none"
		if(tolower(opafit$method.fit) == "mle-rba")  debias <- "rba"
		if(tolower(opafit$method.fit) == "mle-unbias") {
			if(fit_dist == "weibull") {
				debias <- "hrbu"
			}else{
#mle-unbias taken as mle-rba for lognormal
				debias <- "rba"
			}
		}
		fit$conf[[i]]$debias <- debias
		
## re-append the 3p when applicable
		if(npar==3) fit_dist <- paste0(fit_dist, "3p")

## usage FMbounds(x, dist="weibull", CI=.90, unrel=NULL, debias="none", show=FALSE)
		ret<-FMbounds(xdata$lrq_frame, dist=fit_dist, CI=opaconf$ci, unrel=unrel, debias=debias)


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


#############################################################################
############   Likelihood Ratio bounds ######################################
#############################################################################

	if(any(c("lrb","lrbounds") %in% tolower(opaconf$method.conf))) {
	
	if(substr(tolower(fit$options$method.fit),1,3)!= "mle") {								
		stop("likelinood ratio bounds are only applicable on mle fits")							
	}		
		fit$conf[[i]]        <- list()
		fit$conf[[i]]$type   <- "lrb"
		fit$conf[[i]]$ci     <- opaconf$ci
	## It is likely desirable to list input characteristics of the contour
		fit$conf[[i]]$dof <- opaconf$dof
##		fit$conf[[i]]$applyFF     <- opaconf$applyFF

		fit$conf[[i]]$blife.pts <- opaconf$blife.pts
		ret <- NULL
	## this qualifier needs to apply to all bounds
	##	if(any(c("weibull3p", "lognormal3p") %in% tolower(fit$options$dist))) {
	##		stop("confidence bounds are not prepared on 3-parameter fits")
	##	}

## just get the distribution from the fit object, no crossfire
	fit_dist<-fit$options$dist

		debias<-"none"
		if(tolower(opafit$method.fit) == "mle-rba")  debias <- "rba"
		if(tolower(opafit$method.fit) == "mle-unbias") {
			if(fit_dist == "weibull") {
				debias <- "hrbu"
			}else{
#mle-unbias taken as mle-rba for lognormal
				debias <- "rba"
			}
		}
if(!is.null(debias)) fit$conf[[i]]$debias <- debias


## usage LRbounds(x,  dist="weibull", CL=0.9, unrel=NULL,  contour=NULL, dof=1, control=NULL, debias="none", show=FALSE)
		ret<-LRbounds(xdata$lrq_frame,
			dist=fit$options$dist,
			CL=opaconf$ci,
			 unrel=unrel,
# this specific setting broke the code, just depend on default
#			contour=NULL,
			dof=fit$conf[[i]]$dof,
			control=list(ptDensity=opaconf$ptDensity, RadLimit=opaconf$RadLimit, tzpoints=opaconf$tzpoints),
			debias=debias
		)

		if(!is.null(ret)){
			atLeastOneBLifeConf <- TRUE

			ret$bounds$percentile<-ret$bounds$percentile/100

			fit$conf[[i]]$bounds <- ret$bounds
			names(fit$conf[[i]]$bounds) <- c("unrel","Lower","Datum", "Upper")

			fit$conf[[i]]$contour<-ret$contour

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

	}  # end LRB bounds



		if(!atLeastOneBLifeConf)  {
				warning("calculateSingleConf: The fit argument is empty or contains no fits.")
		}

	x$fit[[length(x$fit)]]<-fit
	}
	x



}


