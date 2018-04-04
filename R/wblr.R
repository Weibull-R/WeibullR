# wblr.R
# Inspired by Abrem.R originally written by Jurgen Symynck
# Extensive re-write by David J. Silkworth includes:
#	interval input argument
#	new [object]$data now as list containing objects lrq_frame, dpoints, and dlines
#   the [object]$data$dpoints corresponds to previous [object]$data
# copyright (c) OpenReliability.org 2011-2017
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

wblr<-function(x, s=NULL, interval=NULL,...) {
arg <- splitargs(...)
## arg$opa includes options.wblr list items for update
## arg$rem holds any remaining argument names
### need to split out opa$dat, modified by arg$dat to set data instead of using all modified arg$opa
## depreciate use of log option and validate log or canvas entry
	if(!is.null(arg$opa$log)) {
		warning("log option is to be depreciated in favor of canvas")
		#stop("log option is not to be set directly, use canvas") #future hard stop validation
		if(arg$opa$log=="xy" || arg$opa$log=="yx") {
			arg$opa$canvas<-"lognormal"
			arg$opa$log<-"xy"  #just in case "yx" was entered
		}else{
			if(arg$opa$log=="x") {
				arg$opa$canvas<-"weibull"
			}else{
				stop("if used, log argument must be \"x\", \"xy\", or \"yx\" ")
			}
		}
	}else{
		if(!is.null(arg$opa$canvas)) {
			if(tolower(arg$opa$canvas)=="lognormal") {
				arg$opa$log<-"xy"
				arg$opa$canvas<-"lognormal"
			}else{
				arg$opa$log<-"x"
				if(tolower(arg$opa$canvas)!="weibull") {
					warning("canvas option not recognized, default \"weibull\" is assumed")
					arg$opa$canvas<-"weibull"
				}
			}
		}
	}

# silently fix use of ties as an option
		if(!is.null(arg$rem$ties) && is.null(arg$opa$ties.handler)) arg$opa$ties.handler<-arg$rem$ties

##############################################################
## it is necessary to process the input data using arg$opa$pp
## to permit method.fit and method.conf validations in arg$opa
##############################################################


opa <- modifyList(options.wblr(), arg$opa)


## arg$rem holds any remaining argument names
# for now, if x exists it will be assumed that it is either a vector of exact failure times,
# or a time/event dataframe. The time/event df will be further validated in mleframe.

if(!missing(x)){
	ti <- c(arg$rem$time,arg$rem$fail)
	if(!is.null(ti)) warning("fail times are taken from first argument, time or fail arguments are ignored")
	if (class(x) == "data.frame") {
		su<-c(s,arg$rem$susp)
		if(!is.null(su)) warning("suspension times are taken from first argument dataframe, s or susp arguments are ignored")
	}
}else{
	if(!is.null(arg$rem$fail)) {
		x<-arg$rem$fail
		if(!is.null(arg$rem$time)) {
			warning("fail times are taken from fail argument, time argument has been ignored")
		}
	}else{
		if(!is.null(arg$rem$time)) {
			x<-arg$rem$time
		}
	}
}
if (!class(x) == "data.frame") {
if(!missing(s)){
	if(!is.null(arg$rem$susp)) warning("suspension times are taken from second argument s, argument susp is ignored")
}else{
	if(!is.null(arg$rem$susp)) s<-arg$rem$susp
}
}

lrq_frame<-mleframe(x,s,interval)

plot_data<-getPlotData(lrq_frame,opa)


    obj <- list()
	obj$data <- list()
	obj$data$lrq_frame<-lrq_frame
	obj$data$dpoints<-plot_data[[1]]
	obj$data$dlines<-plot_data[[2]]

## here was source of bug found in Legend upon ties handling
## previous use of the length of dpoints and dlines columns was not accurate
## now using the sum of weight values this is correct for all cases

	obj$fail <- sum(obj$data$dpoints$weight)
	obj$interval <- sum(obj$data$dlines$weight)
	obj$n    <- sum(lrq_frame$qty)
	obj$cens <- obj$n - obj$fail - obj$interval
	obj$discovery <-sum(lrq_frame$qty[lrq_frame$left==0])
	obj$interval <- obj$interval-obj$discovery

    obj$options <- opa
        # always store a full copy of the options.wblr structure here
    class(obj) <- "wblr"

obj
}





getPlotData<-function(x,opa) {

## adjustmets to ppp are made by using options
## it is also possible to alter the rank adjustment from "Johnson" to"KMestimator"
## it is also possible to adjust the handling of ties as covered by getPercentilePlottingPositions
	opa$pp<-tolower(opa$pp)
	if(opa$pp=="median") ppos<-"beta"
	if(opa$pp=="benard") ppos<-"Benard"
	if(opa$pp=="hazen") ppos<-"Hazen"
	if(opa$pp=="mean") ppos<-"mean"
	if(opa$pp=="kaplan-meier") ppos="Kaplan-Meier"
	if(opa$pp=="blom")  ppos<-"Blom"

	if(!(any(c("beta", "Benard", "Hazen", "mean","Kaplan-Meier","Blom") %in% ppos))) {
		stop("pp option not recognized")
	}

	if(tolower(opa$rank.adj)=="johnson")  {
		aranks<-"Johnson"
	}else{
		if(opa$rank.adj==tolower("kmestimator")) {
			if(ppos!="Kaplan-Meier") warning("KMestimator applied without Kaplan-Meier plotting positions")
			aranks<-"KMestimator"
		}else{
			stop("rank.adj option not recognized")
		}
	}

	if(!(any(c("none", "highest", "lowest", "mean","sequential") %in% opa$ties.handler))) {
		stop("ties.handler option not recognized")
	}


	exponent10<-(2+log10(max(x[,1:2])))
	mult<-10^exponent10
## the original reason for forming the mod1x was to enable the transformed (mean) time calculation
## it now also removes suspensions ready for use in forming dataDF
	mod1x<-cbind(x,mean=(x$left+x$right)/2)[x$right!=-1,]
	trans_time<-mod1x$mean*mult+mod1x$left
	dataDF<-data.frame(time=trans_time, event=1, qty=mod1x$qty)
	suspDF<-NULL
	if(length(x[x$right==-1,][,1])>0) {
	sx<-x[x$right==-1,]
	suspDF<-data.frame(time=sx$left*mult, event=0, qty=sx$qty)
	}
	p_argx<-rbind(dataDF,suspDF)

	mod2x<-cbind(mod1x[,-4],tmean=dataDF$time)

	p<-getPPP(p_argx, ppos=ppos, aranks=aranks, ties=opa$ties.handler)

if( nrow(mod2x)==nrow(p) )  {
## this is the R-efficient method when ties have been handled, 
## or there were no duplicates at all
	NDX<-order(mod2x[,4])
	mod3x<-mod2x[NDX,]
	lrdiv <-mod3x$left/mod3x$right
	mod3x<-cbind(mod3x, p$ppp, p$adj_rank, lrdiv)
	dpoints<-mod3x[mod3x$lrdiv==1, c(2,5,6,3)]
	if(nrow(dpoints)==0) {
		dpoints<-NULL
	}else{
		names(dpoints)<-c("time", "ppp", "adj_rank", "weight")
	}
	dlines<-mod3x[mod3x$lrdiv!=1, c(1,2,5,6,3)]
	if(nrow(dlines)==0)  {
		dlines<-NULL
	}else{
		names(dlines)<-c("t1", "t2", "ppp", "adj_rank", "weight")
	}

}else{
## This is the complex loop with un-handled duplicates
	dpoints<-NULL
	dlines<-NULL
#if(opa$ties.handler=="none") {
	for(frow in 1: nrow(mod2x))  {
		if(mod2x$left[frow]==mod2x$right[frow]) {
			for(q in 1:mod2x$qty[frow])  {
				prow<-match(mod2x$tmean[frow], p$time)
				point_row<-data.frame(
					time=mod2x$left[frow], 
					ppp=p$ppp[prow], 
					adj_rank=p$adj_rank[prow],
					weight=1)
				dpoints<-rbind(dpoints, point_row)
				p<-p[-prow,]
			}
		}else{
			for(q in 1:mod2x$qty[frow])  {
				prow<-match(mod2x$tmean[frow], p$time)
				line_row<-data.frame(
					t1=mod2x$left[frow], 
					t2=mod2x$right[frow], 
					ppp=p$ppp[prow], 
					adj_rank=p$adj_rank[prow],
					weight=1)
				dlines<-rbind(dlines, line_row)
				p<-p[-prow,]
			}
		}
	}
}
## This is the loop for the case where tie handling has taken place.
## Since the p vector should be same length and order of mod2x, there is 
## likely an efficient R method for handling this, no prow match required
#}else{
#	for(frow in 1: nrow(mod2x))  {
#		if(mod2x$left[frow]==mod2x$right[frow]) {
#			#for(q in 1:mod2x$qty[frow])  {
#				prow<-match(mod2x$tmean[frow], p$time)
#				point_row<-data.frame(
#					time=mod2x$left[frow], 
#					ppp=p$ppp[prow], 
#					adj_rank=p$adj_rank[prow],
#					weight=mod2x$qty[frow])
#				dpoints<-rbind(dpoints, point_row)
#				p<-p[-prow,]
#			#}
#		}else{
#			#for(q in 1:mod2x$qty[frow])  {
#				prow<-match(mod2x$tmean[frow], p$time)
#				line_row<-data.frame(
#					t1=mod2x$left[frow], 
#					t2=mod2x$right[frow], 
#					ppp=p$ppp[prow], 
#					adj_rank=p$adj_rank[prow],
#					weight=mod2x$qty[frow])
#				dlines<-rbind(dlines, line_row)
#				p<-p[-prow,]
#			#}
#		}
#	}
#}

outlist<-list(dpoints,dlines)
outlist
}


splitargs <- function(...){
    arg         <- list(...)
    argnames    <- names(arg)
    parplot     <- plot_default_args()
    ret         <- list()
    opanames    <- names(options.wblr())
    ret$opa     <- arg[tolower(argnames) %in% tolower(opanames)]
        # ret$opa can be an emply list, which is compatible with modifyList()
    ret$rem     <- arg[!(tolower(argnames) %in% tolower(opanames))]
    ret
}
