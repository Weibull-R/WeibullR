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
###################################################################
##                      this would be the position to go to C++
###################################################################

				secant_warning=FALSE
##  Tao Pang's original variable labels from FORTRAN are used where possible
				DL<-ControlList$limit
## Introduce constraints for tz
## need to create the vector of potential minimums
#				fdr<-NULL
#				if(Nf>0) {fdr<-fsdi[1:Nf]}
#				if(Nd>0) {fdr<-c(fdr,fsdi[(Nf+Ns+1):(Nf+Ns+Nd)])}
#				if(Ni>0)  {fdr<-c(fdr, fsdi[(Nf+Ns+Nd+Ni+1):(Nf+Ns+Nd+2*Ni)])}
				fdr<-NULL
				if(Nf_rows>0) {fdr<-fsdi[1:Nf_rows]}
				if(Nd_rows>0) {fdr<-c(fdr,fsdi[(Nf_rows+Ns_rows+1):(Nf_rows+Ns_rows+Nd_rows)])}
				if(Ni_rows>0)  {fdr<-c(fdr, fsdi[(Nf_rows+Ns_rows+Nd_rows+Ni_rows+1):(Nf_rows+Ns_rows+Nd_rows+2*Ni_rows)])}
				C1<-min(fdr)
				maxit<-100

## initial step is based on min(x)*.1
				DX<-C1*0.1
				X0<-0.0
				istep<-0
				X1<-X0+DX
				if(X1>C1) {X1<-X0+0.9*(C1-X0)}
				tz=0
##				FX0vec<-.Call("MLEdMaxLLdx", MLEclassList, ControlList, vstart, tz, package="WeibullR")
				FX0vec<-.Call(MLEdMaxLLdx, MLEclassList, ControlList, vstart, tz)
				FX0<-FX0vec[1]
## new start estimate from last fit (without any modification)
				vstart<-FX0vec[-1]
## X1 is next proposed tz
				tz=X1
##				FX1vec<-.Call("MLEdMaxLLdx", MLEclassList, ControlList, vstart, tz, package="WeibullR")
				FX1vec<-.Call(MLEdMaxLLdx, MLEclassList, ControlList, vstart, tz)
				FX1<-FX1vec[1]
## new start estimate from last fit (without any modification)
				vstart<-FX1vec[-1]
## FX1 will contain slope sign information to be used only one time to find X2
				D<- abs(FX1-FX0)
				X2<-X1+abs(X1-X0)*FX1/D
				if(X2>C1) {X2<-X1+0.9*(C1-X1)}
				X0<-X1
				X1<-X2
				DX<-X1-X0
				istep<-istep+1
##  Detail output to be available with listout==TRUE
		DF<-data.frame(steps=istep,root=X0,error=DX,deriv=FX1)


		while(abs(DX)>DL&&istep<maxit)  {
				FX0<-FX1
## save last successful call with useful vstart
				FX0vec<-FX1vec
## X1 is next proposed tz
				tz=X1
##				FX1vec<-.Call("MLEdMaxLLdx", MLEclassList, ControlList, vstart, tz, package="WeibullR")
				FX1vec<-.Call(MLEdMaxLLdx, MLEclassList, ControlList, vstart, tz)
				FX1<-FX1vec[1]
				if(is.nan(FX1))  {
				FX1<-FX0
				secant_warning=TRUE
				break
				}
## new start estimate from last fit (without any modification)
				vstart<-FX1vec[-1]

## FX1 will contain slope information only one time
				D<- abs(FX1-FX0)
				X2<-X1+abs(X1-X0)*FX1/D
				if(X2>=C1) {X2<-X1+0.9*(C1-X1)}

				X0<-X1
				X1<-X2
## setting up to test for instability
last_abs_DX<-abs(DX)
				DX<-X1-X0

## Don't understand this instability, but it happened on example testing with interval data
## The multiplier of 2 here might be too tight, but worked well on example
if(abs(DX) > last_abs_DX*2) {
stable<-FALSE
##browser()
break
}
				istep<-istep+1

				DFline<-data.frame(steps=istep,root=X0,error=DX,deriv=FX1)
				DF<-rbind(DF,DFline)
## return to while loop
		}
## provide a last good vstart in case FX1vec indicated nan
		vstart<-FX0vec[-1]
		
## Can X0 be first trial, but ultimately subject to convergence problems??
		listout_int<-0
##		result_of_simplex_call<-.Call("MLEsimplex",MLEclassList, ControlList, vstart, X0, listout_int, package="WeibullR")
		result_of_simplex_call<-.Call(MLEsimplex,MLEclassList, ControlList, vstart, X0, listout_int)

## extract fit vector from result of call to enable finishing treatment of the outvec


		outvec<-c(result_of_simplex_call[1:2], X0, result_of_simplex_call[3])

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
##				outvec[3]<-.Call("MLEloglike",MLEclassList,c(outvec[2],outvec[1]),dist_num, default_sign, X0, package="WeibullR")
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
##				outvec[3]<-.Call("MLEloglike",MLEclassList,c(outvec[1],outvec[2]),dist_num, default_sign, X0, package="WeibullR")
				outvec[3]<-.Call(MLEloglike,MLEclassList,c(outvec[1],outvec[2]),dist_num, default_sign, X0)
				attr(outvec,"bias_adj")<-debias
			}
		}

		optDF<-DF[-1]

## end of 3p code
}

## restored the original secant_warning and added the new stabile warning
## perhaps using if(exists("secant_warning"))  etc. would be better.
#if(secant_warning==TRUE) {
if(exists("secant_warning")) {
	attr(outvec, "warn")<-"unstable fit 1"
}
#if(stable==FALSE) {
if(exists("stable")) {
	attr(outvec, "warn")<-"unstable fit 2"
}

## the following applies to both 2p and 3p results
## it is used by LRbounds to simplify data_type determination for debias adjustment
## but often removed for normal use by attributes(fit_vec)$data_types<-NULL
	attr(outvec,"data_types")<-Q[-2]

	if(listout==FALSE) {
		out_object<-outvec
	}else{
		out_object<-list(fit=outvec, opt=optDF)
	}

	out_object

## end function
}
