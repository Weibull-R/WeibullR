mlefit<-function(x, dist="weibull", npar=2, debias="none", optcontrol=NULL)  {
## these warnings only apply to 3p fitting, but the object must be created
#secant_warning<-FALSE
#stable<-TRUE
## tz is required for MLEloglike and MLEsimplex calls now
		default_tz=0
## sign is now required for MLEloglike call
		default_sign=1

## check basic parameters of x
	if(class(x)!="data.frame") {stop("mlefit takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("mlefit takes a structured dataframe input, use mleframe")}
	xnames<-names(x)
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {
		 stop("mlefit takes a structured dataframe input, use mleframe")  }
## test for any na's and stop, else testint below will be wrong

## It turns out that this code is general to all fitting methods:
	if(tolower(dist) %in% c("weibull","weibull2p","weibull3p")){
		fit_dist<-"weibull"
	}else{
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p", "lognormal3p")){
			fit_dist<-"lnorm"
		}else{
#			if(!dist=="gumbel") {
		## Note:  only lslr contains experimental support for "gumbel"
			stop(paste0("dist argument ", dist, "is not recognized for mle fitting"))
#			}
		}
	}

##	npar<-2 ## introducing 3p in dist argument will override any npar (or its default)
	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){
		npar<-3
	}

## initialize counts at zero, to be filled as found
	Nf=0
	Ns=0
	Nd=0
	Ni=0

## need this length information regardless of input object formation

	failNDX<-which(x$right==x$left)
	suspNDX<-which(x$right<0)
	Nf_rows<-length(failNDX)
	if(Nf_rows>0) {
		Nf<-sum(x[failNDX,3])
	}
	Ns_rows<-length(suspNDX)
	if(Ns_rows>0) {
		Ns<-sum(x[suspNDX,3])
	}
	discoveryNDX<-which(x$left==0)
	Nd_rows<-length(discoveryNDX)
	if(Nd_rows>0) {
		Nd<-sum(x[discoveryNDX,3])
	}
	testint<-x$right-x$left
	intervalNDX<-which(testint>0)
	interval<-x[intervalNDX,]
	intervalsNDX<-which(interval$left>0)
	Ni_rows<-length(intervalsNDX)
	if(Ni_rows>0) {
		Ni<-sum(interval[intervalsNDX,3])
	}


## rebuild input vector from components, because this order is critical
	fsiq<-rbind(x[failNDX,], x[suspNDX,], x[discoveryNDX,], interval[intervalsNDX,])

## Not sure what to place as restriction for C++ call
##	if((Nf+Ni)<3)  {stop("insufficient failure data")}


## now form the arguments for C++ call
## fsdi is the time vector to pass into C++
	fsd<-NULL
	if((Nf+Ns)>0)  {
		fsd<-fsiq$left[1:(Nf_rows + Ns_rows)]
	}
	if(Nd>0) {
		fsd<-c(fsd,fsiq$right[(Nf_rows + Ns_rows + 1):(Nf_rows +  Ns_rows + Nd_rows)])
	}
	if(Ni>0)  {
		fsdi<-c(fsd, fsiq$left[(Nf_rows + Ns_rows + Nd_rows + 1):nrow(fsiq)],
		fsiq$right[(Nf_rows + Ns_rows + Nd_rows + 1):nrow(fsiq)])
	}else{
		fsdi<-fsd
	}

	q<-fsiq$qty
## third argument will be c(Nf,Ns,Nd,Ni)
	N<-c(Nf_rows,Ns_rows,Nd_rows,Ni_rows)

	mrr_fail_data<- c(rep(x[failNDX,1],x[failNDX,3]),
		rep( x[discoveryNDX,2]/2, x[discoveryNDX,3]),
		rep((interval[intervalsNDX,1]+(interval[intervalsNDX,2]-interval[intervalsNDX,1])/2), interval[intervalsNDX,3])
		)
	mrr_susp_data<-rep(x[suspNDX,1], x[suspNDX,3])

## establish distribution number and start parameters
	if(fit_dist=="weibull"){
		dist_num=1
## single failure data set with suspensions (only) uses simplistic weibayes for vstart
		if(Nf==1 && Nd+Ni==0) {
		weibayes_scale <-x[failNDX,1]+sum(x[suspNDX,1])
		vstart<- c(1, weibayes_scale)
		warning("single failure data set may be candidate for weibayes fitting")
		}else{
# use of quick fit could have been circular here
#		mrr_fit<-MRRw2p(mrr_fail_data, mrr_susp_data)
			mrr_fit<-lslr(getPPP(mrr_fail_data, mrr_susp_data))
			shape<-mrr_fit[2]
			scale<- mrr_fit[1]
			vstart <- c(shape, scale)

		}

	}else{
		if(fit_dist=="lnorm"){
			dist_num=2
# use of quick fit could have been circular here
#			mrr_fit<-MRRln2p(mrr_fail_data, mrr_susp_data)
			mrr_fit<-lslr(getPPP(mrr_fail_data, mrr_susp_data), dist="lognormal")
			ml<- mrr_fit[1]
			sdl<- mrr_fit[2]
#			ml <- mean(log(data_est))
#			sdl<- sd(log(data_est))
			vstart<-c(ml,sdl)
		}else{
			stop("distribution not resolved for mle fitting")
		}
	}

## Optional optimization control list to be handled here
		## vstart will be as estimated
		limit<-1e-5
		maxit<-100
		listout<-FALSE

	if(length(optcontrol)>0)  {
		if(length(optcontrol$vstart>0))  {
			vstart<-optcontrol$vstart
		}
		if(length(optcontrol$limit)>0)  {
			limit<-optcontrol$limit
		}
		if(length(optcontrol$maxit)>0)  {
			maxit<-optcontrol$maxit
		}
		if(length(optcontrol$listout)>0)  {
			listout<-optcontrol$listout
		}
	}

	pos<-1
	Q<-sum(q)
	for(j in seq(1,4))  {
		if(N[j]>0) {
			Q<-c(Q, sum(q[pos:(pos+N[j]-1)]))
			pos<-pos+N[j]
		}else{
			Q<-c(Q, 0)
		}
	}
	names(Q)<-c("n","fo", "s", "d", "i")

	MLEclassList<-list(fsdi=fsdi,q=q,N=N)
## Test for successful log-likelihood calculation with given vstart
## tz is required for MLEloglike call now
##		LLtest<-.Call("MLEloglike",MLEclassList,vstart,dist_num, default_sign, default_tz, package="WeibullR")
		LLtest<-.Call(MLEloglike,MLEclassList,vstart,dist_num, default_sign, default_tz)
## This should have failed as left with abremDebias call.
		if(!is.finite(LLtest))  {
			stop("Cannot start mle optimization with given parameters")
		}

	ControlList<-list(dist_num=dist_num,limit=limit,maxit=maxit)

## here is a good place to validate any debias argument (before more calculations begin)
	if(debias!="none" && dist_num==1)  {
		if(tolower(debias)!="rba"&&tolower(debias)!="mean"&&tolower(debias)!="hrbu")  {
			stop("debias method not resolved")
		}
	}




## Handle the original 2 parameter case first
##	if(tolower(dist)=="weibull" || tolower(dist)=="lognormal" ||tolower(dist)=="weibull2p" || tolower(dist)=="lognormal2p"  )  {
	if(npar==2) {

## listout control is passed as an integer to C++, this enables temporary change of status without losing input argument value

			if(listout==TRUE)  {
				listout_int<-1
			}else{
				listout_int<-0
			}
##  tz  inserted here with a default of zero
##		result_of_simplex_call<-.Call("MLEsimplex",MLEclassList, ControlList, vstart, default_tz, listout_int, package="WeibullR")
		result_of_simplex_call<-.Call(MLEsimplex,MLEclassList, ControlList, vstart, default_tz, listout_int)
## extract fit vector from result of call to enable finishing treatment of the outvec
		if(listout==FALSE)  {
			resultvec<-result_of_simplex_call
		}else{
			resultvec<-result_of_simplex_call[[1]]
		}
		outvec<-resultvec[1:3]
		if(resultvec[4]>0)  {
			warn<-"likelihood optimization did not converge"
			attr(outvec,"warning")<-warn
		}

		if(dist_num == 1)  {
			names(outvec)<-c("Eta","Beta","LL")
			if(debias!="none")  {
				if(debias!="rba"&&debias!="mean"&&debias!="hrbu")  {
					stop("debias method not resolved")
				}
				if(debias=="rba")  {
					outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="weibull",basis="median")
				}
				if(debias=="mean")  {
					outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="weibull",basis="mean")
				}
				if(debias=="hrbu")  {
					outvec[2]<-outvec[2]*hrbu(Q[1]-Q[3], Q[3])
				}
##			outvec[3]<-.Call("MLEloglike",MLEclassList,c(outvec[2],outvec[1]),dist_num, default_sign, default_tz, package="WeibullR")
			outvec[3]<-.Call(MLEloglike,MLEclassList,c(outvec[2],outvec[1]),dist_num, default_sign, default_tz)
			attr(outvec,"bias_adj")<-debias
			}
		}

		if(dist_num == 2)  {
			names(outvec)<-c("Mulog","Sigmalog","LL")
			if(debias!="none")  {
				outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="lognormal")
				if(debias!="rba")  {
					warning("rba has been applied to adjust lognormal")
					debias="rba"
				}
##			outvec[3]<-.Call("MLEloglike",MLEclassList,c(outvec[1],outvec[2]),dist_num, default_sign, default_tz, package="WeibullR")
			outvec[3]<-.Call(MLEloglike,MLEclassList,c(outvec[1],outvec[2]),dist_num, default_sign, default_tz)
			attr(outvec,"bias_adj")<-debias
			}
		}


		if(listout==TRUE) {
			optDF<-as.data.frame(result_of_simplex_call[[2]])
			if(dist_num == 1)  {
				names(optDF)<-c("beta_est", "eta_est", "negLL", "error")
			}
			if(dist_num == 2)  {
				names(optDF)<-c("mulog_est", "sigmalog_est", "negLL", "error")
			}
		}

## end of 2p code
	}



##  this section of code is specifically addressing 3p models
##	if(tolower(dist)=="weibull3p" || tolower(dist)=="lognormal3p"  )  {
	if(npar==3) {

## For now, listout is passed as integer
## listout argument has different meaning for 3p models
		listout_int<-0

## for now enter a default tz=0
##			result_of_simplex_call<-.Call("MLEsimplex",MLEclassList, ControlList, vstart, default_tz, listout_int, package="WeibullR")
			result_of_simplex_call<-.Call(MLEsimplex,MLEclassList, ControlList, vstart, default_tz, listout_int)
			if(result_of_simplex_call[4]>0)  {
				stop("2p model does not converge")
			}
## restore the meaning of listout
			if(listout==TRUE)  {
				listout_int<-1
			}else{
				listout_int<-0
			}

			
			
	## n has been removed as an argument and placed in the control list						
	seek_control<-list(						
		num_points=10,					
		err_t0_limit= 1e-6,					
		err_gof_limit= 1e-5)
	optcontrol4seek<-list(optcontrol$num_points, optcontrol$err_t0_limit, optcontrol$err_gof_limit)
	if(length(optcontrol4seek)>0) seek_control<-modifyList(seek_control, optcontrol4seek)						
	## restore meaning of n for rest of code						
	n<-seek_control$num_points						
	if(n<5) {						
		n<-5					
		warning("n specified too small, n=5 used")					
	}						
							
## This is the point to go to C++						
## Will need to pass in MLEclassList, fit_dist and seek_control							
							
							
	# This function is called in context on each trial to establish the sequence of tz points						
	#  no arguments, start, end and maxtz labels are drawn from the environment at the time of call						
	tzvec<-function() {						
	# we either protect the maxtz end point						
		if(end == maxtz)  {					
			spacing<-(end-start)/n  ## protecting maxtz				
			return(seq(start, end-spacing, by=spacing))				
		}else{					
	# or we get a sequence including the end point						
	# spacing will always take the sign of the end point						
			spacing<-(end-start)/(n-1)				
			return(seq(start, end, by=spacing))				
		}					
	}						
							
							
							
	# establish the maximum limit for t0						
	## MLEmodel will treat convert any negative x$left-tz as zero						
	## Note discoveries continue to be discoveries until x$right-tz = zero						
	maxtz<-min(x$right[x$right != -1])						
							
## this may not be necessary as MLEmodel will treat any negative x$left-tz as zero							
##	## short circuit if any discoveries are in intervals						
##	t0_found<-FALSE						
##	if(maxtz==0)  {						
##			t0_found<-TRUE				
##			t0<-0				
##			fit<-mlefit(data_mod,dist="weibull")				
##			DF<-data.frame(P1=fit[1],P2=fit[2],tz,gof=fit[3])				
##	}						
							
	# set up for first trial						
	start<-0						
	end<-maxtz						
							
							
	t0_found<-FALSE						
	rebound<-FALSE						
	trial<-1						
	try_list<-list()						
## the simplex optcontrol is set to the defaults in mlefit							
	optcontrol<-list(limit=1e-5, maxit=100)						
							
while(!t0_found)  {							
	DF<-NULL						
							
## This is the view generation loop that establishes each trial							
	for(tz in tzvec()) {						
		right_mod<-sapply(x$right, function(X) if(X>0) X-tz else X)					
		data_mod<-data.frame(left=x$left-tz, right=right_mod, qty=x$qty)					
## note MLEmodel::Loglike ignores negative suspensions and will convert early intervals to discoveries as x$left-tz becomes =< 0							
		fit<-mlefit(data_mod,dist=fit_dist,optcontrol=optcontrol)					
		thisDFrow<-data.frame(P1=fit[1], P2=fit[2],tz,gof=fit[3])					
		DF<-rbind(DF,thisDFrow)					
	}						
							
	try_list[[trial]]<-DF						
## The rest of the loop is sorting out whether an optimum has been identified							
## or then setting up the next trial							
							
	max_ind<-which(DF$gof==max(DF$gof))	
## in case the max can't be defined by a single element of DF, we must be done	
if(length(max_ind)>1) {
	max_ind<-max_ind[1]
	t0_found<-TRUE
}else{
	if(max_ind !=1)  {						
## get error measures for tz and gof by comparison with max_ind-1, for use in various locations as well as output							
		err_t0<- abs((DF$tz[max_ind]-DF$tz[max_ind-1])/(DF$tz[max_ind]))					
		err_gof<-abs((DF$gof[max_ind]-DF$gof[max_ind-1])/(DF$gof[max_ind]))					
	}						
							
	if(max_ind != 1 && max_ind !=n) {						
		if( err_t0 > seek_control$err_t0_limit) {					
							
			## establish optcontrol for next trial if necessary				
			if(err_gof/n < 1e-5) optcontrol$limit=err_gof/n				
				start<-DF$tz[max_ind-1]			
				end<-DF$tz[max_ind+1]			
		}else{					
			t0_found<-TRUE				
		}					
	}else{						
		if(trial==1 && max_ind==1)  {					
		## reverse the seek for negative t0					
			start<- DF$tz[2]		# make sure to cover all early positve values		
			end<-(-1)*maxtz				
		}else{					
			if(DF$tz[n] >0)  {				
				if(max_ind == n)  {			
					## get err_t0, then check if less than default limit		
					if(err_t0 < seek_control$err_t0_limit) {		
						t0_found<-TRUE
						positive_runnout<-TRUE
					}else{		
						## establish optcontrol for next trial if necessary	
						if(err_gof/n < 1e-5) optcontrol$limit=err_gof/n	
					## check for rebound case		
						shift_gof<-c(DF$gof[1],DF$gof[-n])	
						progression<-DF$gof-shift_gof	
						rebound_ind<-which(progression==min(progression))	
						if(rebound_ind!=1)  {	
						## if so, next trial is a rework of this trial with end set to rebound point	
							## start is unchanged
							end<-DF$tz[rebound_ind]
							## decrement the trial so it will replace last and flag this unusual event (for further study?).
							trial<-trial-1
							rebound<-TRUE
						}else{	
							start<-DF$tz[n]
							end<-maxtz
						}	
					}		
				}else{			
				## max_ind must be 1, the seek is still positive t0			
					start <- try_list[[trial-1]]$tz[n-1]		
					end <- DF$tz[2]		
				}			
			}else{				
			## get err_gof, then check if less than default limit				
				if(err_gof < seek_control$err_gof_limit) {			
					t0_found<-TRUE
					negative_runnout<-TRUE
									}else{			
				## establish next mlefit limit if necessary			
				if(err_gof/n < 1e-5) optcontrol$limit=err_gof/n			
							
	## Code modification required here if negative rebound is a concern						
							
							
					start<-end		
					end<-end*10		
				}			
			}				
		}					
	}						
	if(!t0_found)  {						
		trial<- trial+1					
	}
#close the new block for more than one max_ind found	
}
#close the main loop finding t0							
}	

## must collect outvec and try_list from return of C++ call						
	#outvec<-DF[max_ind,]	
	# returning a single line dataframe causes later problems, needs to be a named vector	
	outvec<-c(DF$P1[max_ind], DF$P2[max_ind], DF$tz[max_ind],DF$gof[max_ind])	
	##names(outvec)<-""	
			
	if(dist_num==1)  {
		names(outvec)<-c("Eta","Beta", "t0", "LL")
		if(debias!="none")  {
			if(debias=="rba")  {
				outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="weibull",basis="median")
			}
			if(debias=="mean")  {
				outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="weibull",basis="mean")
			}
			if(debias=="hrbu")  {
				outvec[2]<-outvec[2]*hrbu(Q[1]-Q[3], Q[3])
			}
			outvec[3]<-.Call(MLEloglike,MLEclassList,c(outvec[2],outvec[1]),dist_num, default_sign, X0)
			attr(outvec,"bias_adj")<-debias
		}
	}
	if(dist_num == 2)  {
		names(outvec)<-c("Mulog","Sigmalog", "t0", "LL")
		if(debias!="none")  {
			outvec[2]<-outvec[2]*rba(Q[1]-Q[3], dist="lognormal")
			if(debias!="rba")  {
				warning("rba has been applied to adjust lognormal")
				debias="rba"
			}
			outvec[3]<-.Call(MLEloglike,MLEclassList,c(outvec[1],outvec[2]),dist_num, default_sign, X0)
			attr(outvec,"bias_adj")<-debias
		}
	}

if(exists("positive_runnout")) {
	attr(outvec, "message")<-"t0 cutoff at minimal change"
}

if(exists("negative_runnout")) {
	attr(outvec, "message")<-"optimum not found, t0 cutoff at minimal gof change"
}
## end of 3p code
}
## the following applies to both 2p and 3p results
## it is used by LRbounds to simplify data_type determination for debias adjustment
## but often removed for normal use by attributes(fit_vec)$data_types<-NULL
	attr(outvec,"data_types")<-Q[-2]

	if(listout==FALSE) {
		out_object<-outvec
	}else{
		if(npar==2)  out_object<-list(fit=outvec, opt=optDF)
		if(npar==3)  out_object<-list(fit=outvec, opt=try_list)
	}

	out_object

## end function
}
