wblrLoglike<-function(par, x, dist="weibull", sign=1, tz=0 )  {				
## check basic format of x				
				
	if(class(x)!="data.frame") {stop("wblrLoglike takes a structured dataframe input, use mleframe")}			
	if(ncol(x)!=3)  {stop("wblrLoglike takes a structured dataframe input, use mleframe")}			
	xnames<-names(x)			
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {			
		 stop("wblrLoglike takes a structured dataframe input, use mleframe")  }		
## test for any na's and stop, else testint below will be wrong				
				
				
## need this length information regardless of input object formation				
	testint<-x$right-x$left			
	failNDX<-which(testint==0)			
	suspNDX<-which(testint<0)			
	Nf<-length(failNDX)			
	Ns<-length(suspNDX)			
	discoveryNDX<-which(x$left==0)			
	Nd<-length(discoveryNDX)			
	intervalNDX<-which(testint>0)			
	interval<-x[intervalNDX,]			
	intervalsNDX<-which(interval$left>0)			
	Ni<-length(intervalsNDX)			
				
				
## need to stop if Nf<1?				
## or Nf+Ndi-Nd <3?				
				
## further validate the input arguments for non-frame.fsiq object				
	if(length(attributes(x)$fsiq)!=1)  {			
				
				
				
## stop if Nf+Ns+Ndi != nrow(x)				
	if( (Nf+Ns+Nd+Ni) != nrow(x))  {			
		stop("invalid input dataframe")		
	}			
				
## rebuild input vector from components, just to be sure				
	fsiq<-rbind(x[failNDX,], x[suspNDX,], x[discoveryNDX,], interval[intervalsNDX,])			
## end input validation code				
	}else{			
		fsiq<-x		
	}			
				
## now form the arguments for C++ call				
## no data limitation applies to getting a Loglikelihood value	
##	if((Nf+Ni)<3)  {stop("insufficient failure data")}	

	fsdi<-NULL
	if( (Nf+Ns)>0 )  {
		fsdi<-fsiq$left[1:(Nf + Ns)]
	}	
	if(Nd>0) {		
		fsdi<-c(fsdi,fsiq$right[(Nf + Ns + 1):(Nf +  Ns + Nd)])	
	}		
	if(Ni>0)  {		
		fsdi<-c(fsdi, fsiq$left[(Nf + Ns + Nd + 1):nrow(fsiq)], 	
			  fsiq$right[(Nf + Ns + Nd + 1):nrow(fsiq)])	
	}
	
## qualify the tz argument		
	if(tz>0)  {		
		fdr<-NULL	
		if(Nf>0) {fdr<-fsdi[1:Nf]}	
		if(Nd>0) {fdr<-c(fdr,fsdi[(Nf+Ns+1):(Nf+Ns+Nd)])}	
		if(Ni>0)  {fdr<-c(fdr, fsdi[(Nf+Ns+Nd+Ni+1):(Nf+Ns+Nd+2*Ni)])}	
			
		if(tz>min(fdr))  {	
			stop("tz is greater than data permits")
		}	
	}		
	
	q<-fsiq$qty			
## third argument will be c(Nf,Ns,Nd,Ni)				
	N<-c(Nf,Ns,Nd,Ni)	

## establish distribution number, not sure the 3p qualification would make any sense here.
	if(tolower(dist) %in% c("weibull","weibull2p","weibull3p")){
	dist_num=1
	}else{
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p", "lognormal3p")){
			dist_num=2
		}else{
			stop("distribution not resolved for loglike calculation")
		}
	}

	MLEclassList<-list(fsdi=fsdi,q=q,N=N,dist_num=dist_num)
	
	if(sign^2!=1)  {	
		stop("sign must be 1 or -1")
	}	

							
##	outval<-.Call("MLEloglike",MLEclassList, par, sign, tz, package="WeibullR")
 	outval<-.Call(MLEloglike,MLEclassList,par, sign, tz)				
			
outval
			
}