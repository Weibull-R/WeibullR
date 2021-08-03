simplex<-function(x, dist="weibull", tz=0, debias="none", optcontrol=NULL)  {
## tz is required for MLEloglike and MLEsimplex calls now
		default_tz=0
## sign is now required for MLEloglike call
		default_sign=1

## check basic parameters of x
	if(class(x)!="data.frame") {stop("mlefit takes a structured dataframe input, use mleframe")}
	if(ncol(x)!=3)  {stop("mlefit takes a structured dataframe input, use mleframe")}
	#if(x$right[1] != min(x$right[x$right != -1])) {stop("use mleframe to sort data")
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
			mrr_fit<-lslr(getPPP(mrr_fail_data, mrr_susp_data), abpval=FALSE)
			shape<-mrr_fit[2]
			scale<- mrr_fit[1]
			vstart <- c(shape, scale)

		}

	}else{
		if(fit_dist=="lnorm"){
			dist_num=2
# use of quick fit could have been circular here
#			mrr_fit<-MRRln2p(mrr_fail_data, mrr_susp_data)
			mrr_fit<-lslr(getPPP(mrr_fail_data, mrr_susp_data), dist="lognormal", abpval=FALSE)
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
		limit<-1e-6
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

	MLEclassList<-list(fsdi=fsdi,q=q,N=N,dist_num=dist_num)
## Test for successful log-likelihood calculation with given vstart
## tz is required for MLEloglike call now
		LLtest<-.Call("MLEloglike",MLEclassList,vstart,default_sign, tz, package="WeibullR")
##		LLtest<-.Call(MLEloglike,MLEclassList,vstart,dist_num, default_sign, default_tz)
## This should have failed as left with abremDebias call.
		if(!is.finite(LLtest))  {
			stop("Cannot start mle optimization with given parameters")
		}

	#ControlList<-list(dist_num=dist_num,limit=limit,maxit=maxit)
	ControlList<-list(limit=limit,maxit=maxit)

## here is a good place to validate any debias argument (before more calculations begin)
	if(debias!="none" && dist_num==1)  {
		if(tolower(debias)!="rba"&&tolower(debias)!="mean"&&tolower(debias)!="hrbu")  {
			stop("debias method not resolved")
		}
	}

#LLtest

## For now, listout is passed as integer
## listout argument has different meaning for 3p models
		listout_int<-0

## for now enter a default tz=0
			result_of_simplex_call<-.Call("MLEsimplex",MLEclassList, ControlList, vstart, tz, listout_int, package="WeibullR")
##			result_of_simplex_call<-.Call(MLEsimplex,MLEclassList, ControlList, vstart, default_tz, listout_int)
			if(result_of_simplex_call[4]>0)  {
				warning("simplex does not converge")
			}
result_of_simplex_call
}
