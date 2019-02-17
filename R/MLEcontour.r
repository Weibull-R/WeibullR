MLEcontour<-function(x,  dist="weibull", CL=0.9,dof=1,MLLx=NULL,MLEfit=NULL, RadLimit=1e-5,
		ptDensity=120, debias="none", show=FALSE)  {

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
		## Note:  only lslr contains experimental support for "gumbel"
			stop(paste0("dist argument ", dist, "is not recognized for mle fitting"))
		}
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
		Ni<-sum(x[intervalsNDX,3])
	}

## rebuild input vector from components, because this order is critical
	fsiq<-rbind(x[failNDX,], x[suspNDX,], x[discoveryNDX,], interval[intervalsNDX,])

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

## establish distribution number
	if(fit_dist=="weibull"){
		dist_num=1
	}else{
		if(fit_dist=="lnorm"){
			dist_num=2
		}else{
			stop("distribution not resolved for mle fitting")
		}
	}

	MLEclassList<-list(fsdi=fsdi,q=q,N=N)

## start of main contour procedures
	if(is.null(MLEfit))  {
##		require(WeibullR)
## in this case the x argument is already an mleframe
		if(dist_num==1) {
			MLEfit<-unname(suppressWarnings(mlefit(x)))
		}else{
			MLEfit<-unname(suppressWarnings(mlefit(x, dist="lognormal")))
		}
	}else{
		unname(MLEfit)
	}

if(dist_num==1)  {
## par is provided as a vector c(shape, scale) for weibull
	par_hat <- c(MLEfit[2], MLEfit[1])
}else{
## par is provided as a vector c(meanlog, sdlog) for lognormal
	par_hat <- c(MLEfit[1], MLEfit[2])
}

## MLLx argument is intended for use in establishing confidence intervals for 3-parameter models.
if(is.null(MLLx)){MLLx<-MLEfit[3]}
##browser()
	ratioLL  <-  MLLx- qchisq(CL,dof)/2
## assure ptDensity is an integer
	ptDensity<-ceiling(ptDensity)

## Call the C++ code to deliver a matrix of contour points
##	resultMat<- .Call("getContour", MLEclassList, par_hat, dist_num, MLLx, ratioLL, RadLimit, ptDensity, package="WeibullR")
	resultMat<- .Call(getContour, MLEclassList, par_hat, dist_num, MLLx, ratioLL, RadLimit, ptDensity)
	if(sum(resultMat[,3])>0) {
		warning("instability detected")
	}
	contourpts<-data.frame(resultMat[,1:2])



	if(dist_num==1) {
		names(contourpts)<- c("Eta", "Beta")
		if(debias!="none")  {

			if(!(debias=="rba" || debias=="mean" || debias=="hrbu"))  {
				stop("debias method not resolved")
			}
			if(debias=="rba")  {
				bias_adj<-rba(sum(N)-Ns, dist="weibull",basis="median")
			}
			if(debias=="mean")  {
				bias_adj<-rba(sum(N)-Ns, dist="weibull",basis="mean")
			}
			if(debias=="hrbu")  {
				bias_adj<-hrbu(sum(N)-Ns, Ns)
			}

			minbeta<-min(contourpts[,2])
			mindelta<-minbeta-minbeta*bias_adj
			contourpts$Beta<-contourpts$Beta-mindelta

			attr(contourpts,"bias_adj")<-debias
		}
	}else{
		names(contourpts)<- c("Mulog", "Sigmalog")
			if(debias!="none")  {
				bias_adj<-rba(sum(N)-Ns, dist="lognormal")
#				minsig<-min(contourpts[,2])
#				maxsig<-max(contourpts[,2])

				contourpts$Sigmalog<-contourpts$Sigmalog*bias_adj
#contourpts$Sigmalog<-sapply(contourpts$Sigmalog,
#				function(x) {
				#x*bias_adj
## eventhough this attempted adjustment created unexpectedly large alteration
## contours basically agreed with SuperSMITH, but lognormal canvas plot was unchanged
#				x*bias_adj*(1+((x-minsig)/(maxsig-minsig))*.05)
#				}
#				)
				if(debias!="rba")  {
					warning("rba has been applied to adjust lognormal")
					debias="rba"
				}
			attr(contourpts,"bias_adj")<-debias
			}
	}

	if(show==TRUE)  {

		maxYpar<-max(contourpts[,2])
		minYpar<-min(contourpts[,2])
		minXpar<-min(contourpts[,1])
		maxXpar<-max(contourpts[,1])

		ylo<-floor(minYpar)
		yhi<-floor(maxYpar)+1

		XparDec<-10^(floor(log(minXpar)/log(10))-1)
		xlo<-XparDec*(floor(minXpar/XparDec)-1)
		xhi<-XparDec*(floor(maxXpar/XparDec)+1)
		if(dist_num==1)  {
			Beta<-MLEfit[2]
			Eta<-MLEfit[1]
			plot(Eta,Beta,xlim=c(xlo,xhi),ylim=c(ylo,yhi), col="red")
		}else{
			Sigmalog<-MLEfit[2]
			Mulog<-MLEfit[1]
			plot(Mulog,Sigmalog,xlim=c(xlo,xhi),ylim=c(ylo,yhi), col="red")
		}
	
		lines(contourpts)

	}

	contourpts
}

