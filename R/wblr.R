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

##############################################################
## it is necessary to process the input data using arg$opa$pp
## to permit method.fit and method.conf validations in arg$opa
##############################################################
## additional features yet to be implemented could include
## handling of ties and rank adjustment method (i.e. "johnson" or "kmestimator")
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


	obj$fail <- length(obj$data$dpoints$time)
##    obj$cens <- length(s)
	obj$interval <- length(obj$data$dlines$t1)
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
	row.names(x)<-seq(1, dim(x)[1])
	f<-NULL
	s<-NULL
	m<-NULL
	for(frow in 1: dim(x)[1])  {

		if(x$left[frow]==x$right[frow]) {
			for(q in 1:x$qty[frow])  {
				f<-c(f,x$left[frow])
			}
			m<-c(m,x$left[frow])
		}else{

			if(x$right[frow]<0) {
				for(q in 1:x$qty[frow])  {
					s<-c(s,x$left[frow])
				}
				m<-c(m, NA)
			}else{
				for(q in 1:x$qty[frow])  {
					f<-c(f,mean(c(x$left[frow],x$right[frow])))
				}
				m<-c(m,mean(c(x$left[frow],x$right[frow])))
			}
		}

	}

	x<-cbind(x,mean=m)
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
#		if(!(any(c("none", "highest", "lowest", "mean") %in% opa$ties.handler))) {
			stop("ties.handler option not recognized")
		}

	p<-getPPP(f,s, ppos=ppos, aranks=aranks, ties=opa$ties.handler)
## set quantity of duplicated points to 1, so that dpoints and dlines loops below work
	if(opa$ties.handler!="none") {
		if(nrow(p)!=nrow(x)) {
			stop("ties.handler did not function as expected")
		}
		x$qty<-rep(1,nrow(x))
	}

	dpoints<-NULL
	dlines<-NULL
	##   dpoints<-data.frame(time=NULL, ppp=NULL)
	##   dlines<-data.frame(t1=NULL, t2=NULL,, ppp=NULL)

	for(frow in 1: nrow(x))  {
		if(x$left[frow]==x$right[frow]) {
			for(q in 1:x$qty[frow])  {
				prow<-match(x$mean[frow], p$time)
				point_row<-data.frame(time=x$mean[frow], ppp=p$ppp[prow], adj_rank=p$adj_rank[prow])
				dpoints<-rbind(dpoints, point_row)
				p<-p[-prow,]
			}

		}else{
			if(!is.na(x$mean[frow]))  {
				for(q in 1:x$qty[frow])  {
					prow<-match(x$mean[frow], p$time)
					line_row<-data.frame(t1=x$left[frow], t2=x$right[frow], ppp=p$ppp[prow], adj_rank=p$adj_rank[prow])
					dlines<-rbind(dlines, line_row)
					p<-p[-prow,]
				}
			}
		}
	}


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
