MLEcontour<-function(x,  dist="weibull", CL=0.9,DF=1,MLEfit=NULL, RadLimit=1e-5,
		ptDensity=120, debias=NULL, applyFF=FALSE, show=FALSE)  {
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
		require(WeibullR)
## in this case the x argument is already an mleframe
		MLEfit<-unname(mlefit(x))
	}else{
		unname(MLEfit)
	}
## Eta_hat and Beta_hat are plotting coordinates for show=TRUE
		Beta_hat<-MLEfit[2]
		Eta_hat<-MLEfit[1]

## par is provided as a vector c(shape, scale)
	par_hat <- c(MLEfit[2], MLEfit[1])
	MLLx<-MLEfit[3]

	FF<-1
	if(applyFF==TRUE) {
## The 'Fulton Factor' is a non-achademic component discussed in Abernethy's book.
## It was a further adjustement, beyond RBA, to add a "pleasing" shape to contour bounds
## when comparing them to pivotal rank regression interval bounds.
## This factor has never been demonstrated to have any statistical basis.
		if(debias!="rba")  {
			stop('FF is only applicable when debias is set to "rba"')
		}else{
			Nf<-length(x)
			FF<-(Nf-1)/(Nf+0.618)
		}
	}
	ratioLL  <-  MLLx- qchisq(CL,DF)/(2*FF)

## assure ptDensity is an integer
	ptDensity<-ceiling(ptDensity)

## Test for successful identification of a contour point at angle theta
	resultMat<- .Call("getContour", MLEclassList, par_hat, dist_num, MLLx, ratioLL, RadLimit, ptDensity, package="WeibullR")

##	resultMat
#}

	if(sum(resultMat[,3])>0) {
		warning("instability detected")
	}
	contourpts<-data.frame(resultMat[,1:2])
	names(contourpts)<- c("Eta", "Beta")

	if(show==TRUE)  {

		maxBeta<-max(contourpts[,2])
		minBeta<-min(contourpts[,2])
		minEta<-min(contourpts[,1])
		maxEta<-max(contourpts[,1])

		ylo<-floor(minBeta)
		yhi<-floor(maxBeta)+1

		EtaDec<-10^(floor(log(minEta)/log(10))-1)
		xlo<-EtaDec*(floor(minEta/EtaDec)-1)
		xhi<-EtaDec*(floor(maxEta/EtaDec)+1)

		if(!exists("ba")) ba<-1
		plot(Eta_hat,Beta_hat*ba,xlim=c(xlo,xhi),ylim=c(ylo,yhi), col="red")
		lines(contourpts)

	}

	contourpts
}

