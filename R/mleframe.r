# mleframe.R
# loaded from preparatory package abremDebias originally written by David J Silkworth
# copyright (c) OpenReliability.org 2014-2017
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

mleframe<-function(x, s=NULL, interval=NULL)  {
## interval dataframe validation
	colname_error<-FALSE
	if(class(interval)=="data.frame")  {
## test names in first two columns
		test_names<-names(interval)
			if(test_names[1] !="left") {
				stop("'left' column name error in interval dataframe object")
			}
			if(test_names[2] !="right") {
				stop("'right' column name error in interval dataframe object")
			}
## additional validations on interval argument, such as positive numeric checking
## removal of potential na's, etc. could take place here
		if(anyNA(interval))  {
		stop("NA not handled in interval data")
		}
		if(any(c(interval$left,interval$right)<0)) {
		stop("negative values in interval data")
		}
		if(any((interval$right-interval$left)<=0))  {
		stop("non-positive interval")
		}
		
## add qty column if not provided
		if(ncol(interval)<3)  {

			ivalchar<- apply(interval,2,as.character)
			ivalstr<-paste0(ivalchar[,1],"_",ivalchar[,2])
			ivaldf<-as.data.frame(table(ivalstr))
			ivalstr2<-as.character(levels(ivaldf[,1]))
## much done here, but this returns the tabled left and right columns
## in a dataframe with rows corresponding to the tabled quantities
			lrdf<-data.frame(
				matrix(
					as.numeric(
						unlist(
							strsplit(ivalstr2,"_")
						)
					)
				,ncol=2, byrow=T
				)
			)
## now just complete the consolidation of duplicates in the interval dataframe
			interval<-cbind(lrdf,ivaldf[,2])
			names(interval)<-c("left","right","qty")

			# interval<- cbind(interval, qty=c(rep(1,nrow(interval))))
		}else{
## assure that a "qty" column exists (and is only extra column used)
			if(is.null(x$qty)) {
				stop("'qty' column name error in interval dataframe object")
			}
## strip any extraneous columns
			interval<-cbind(interval$time, interval$event, interval$qty)
		}

## sort to facilitate consolidation of any duplicated entries, may not be required
		NDX<-order(interval$left,interval$right)
		interval<-interval[NDX,]

## finally, reject any other object type but NULL
	}else{
		if(length(interval)>0)  {
			stop("error in interval argument type")
		}
	}

## now build dataframes for failures and suspensions
## could x be a dataframe with time and event columns??
	suspensions<-NULL
	if(is.vector(x))  {
		if(anyNA(x))  {
		stop("NA in failure data")
		}
		if(any(x<=0))  {
		stop("non-positive values in failure/occurrence data")
		}

		#x<-sort(x)
## I'm not convinced this needs to be sorted here, but it doesn't hurt
		fail_vec<-sort(x)

		# failures<-data.frame(left=x,right=x,qty=rep(1,length(x)))

		if(length(s)>0)  {
		if(anyNA(s))  {
		stop("NA  in suspension data")
		}
		if(any(s<=0))  {
		stop("non-positive values in suspension data")
		}
		susp_vec<-sort(s)

		# suspensions<-data.frame(left=s,right=-1,qty=rep(1,length(s)))
		}
## end pure vector argument processing
	}else{
	## here a time-event dataframe can be evaluated, if provided as x
	## This is the support for a time-event dataframe
		if (class(x) == "data.frame") {

		## this test is drawn from Abrem.R
			if(is.null(x$time) || is.null(x$event)){
				stop(': Argument \"x\" is missing $time and/or ","$event columns...')
			}

	## verify positive time values
			if (anyNA(x$time)) {
				stop("NA in failure or suspension data")
			}
			if (any(x$time<= 0)) {
				stop("non-positive values in failure or suspension data")
			}
	## verify 1's and 0's only in event
	## using Jurgen's validation code
			ev_info <- levels(factor(x$event))
			if(identical(ev_info,c("0","1")) || identical(ev_info,"1")){
			# okay x is holding event indicators
			}else{
			stop("event column not '1' or '0' ")
			}

			if(length(s)>0)  {
			warning("argument 's' ignored when time-event dataframe provided")
			}

			if(is.null(x$qty)) {
				fail_vec<-x$time[x$event==1]
				# failures <- data.frame(left = f, right = f, qty = rep(1, length(f)))
			}else{
	## The assumption is that data input with a qty field is appropriately  consolidated
	## But let's be sure the qty field is all integer, else future havoc could ensue
				if(any(!is.integer(x$qty))) x$qty<-ceiling(x$qty)
				failures <- data.frame(left = f, right = f, qty = x$qty[x$event==1])
			}

			if(identical(ev_info, c("0","1"))) {
			s<-x$time[x$event==0]
				if(is.null(x$qty)) {

					susp_vec<-s
						# suspensions <- data.frame(left = s, right = -1, qty = rep(1, length(s)))
					}else{
	## The assumption is that data input with a qty field is appropriately  consolidated
						suspensions <- data.frame(left = s, right = -1, qty = x$qty[x$event==0])
					}
			}
		}else {
			if (length(x) > 0) {
				stop("error in x argument type")
			}
		}
## end the time_event dataframe  evaluation and close argument processing
	}

## consolidate duplicates in any pure failure or suspension vectors

	if(exists("fail_vec")) {
		fdf<-as.data.frame(table(fail_vec))
		ft<-as.numeric(levels(fdf[,1]))
		fq<-fdf[,2]
		failures<-data.frame(left=ft, right=ft, qty=fq)
	}

	if(exists("susp_vec")) {
		sdf<-as.data.frame(table(susp_vec))
		st<-as.numeric(levels(sdf[,1]))
		sq<-sdf[,2]
		suspensions<-data.frame(left=st, right=-1, qty=sq)
	}

	outDF<-rbind(failures,suspensions,interval)

	outDF
}
