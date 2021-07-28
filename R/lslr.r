lslr<-function(x, dist="weibull", npar=2, reg_method="XonY", abpval=TRUE)  {
	## a convergence limit is fixed here for 3rd parameter  convergence
	## no longer an argument for the R function, but still an argument to C++ functions
	limit<-1e-5

	if(is.vector(x))  {
		stop("use MRR functions for casual fitting, or pre-process with getPPP")
	}else{
		if(names(x)[1]=="time"&&names(x)[2]=="ppp")  {
		## will handle the output from getPPP
		}else{
			if(length(x$ppp)<3)  {
				stop("insufficient failure points")
			}else{
				stop("input format not recognized")
			}
		}
	}

## It turns out that this code is general to all fitting methods:
	if(tolower(dist) %in% c("weibull","weibull2p","weibull3p")){
		fit_dist<-"weibull"
	}else{
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p", "lognormal3p")){
			fit_dist<-"lnorm"
		}else{
			if(!dist=="gumbel") {
		## Note: Original LSLR contained experimental support for "gumbel"
			stop(paste0("dist argument ", dist, "is not recognized for distribution fitting"))
			}
			if(dist=="gumbel") {
		## Note: experimental support for "gumbel" is currently suspended
			stop("experimental support for gumbel is currently suspended")
			}			
		}
	}

##	npar<-2 ## introducing 3p in dist argument will override any npar (or its default)
	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){
		npar<-3
	}
	if(!(npar %in% c(2,3))) stop("unsupported npar value")

	##casenum<-0
	reg_order=0
	if(reg_method=="YonX") reg_order=1
	##if(npar==3) casenum=casenum+2
	dist_num=0
	if(fit_dist=="lnorm")dist_num=1
	if(dist=="gumbel") dist_num=3

	parlist<-list(fail=x$time, ppp=x$ppp, reg_order=reg_order, dist_num=dist_num, npar=npar, limit=limit)
resultVec<-.Call("LSLR", parlist, package="WeibullR")
##resultVec<-.Call(LSLR, parlist)					
				
					
casenum<-0					
if(reg_method=="YonX") casenum=casenum+1					
if(npar==3) casenum=casenum+2					
if(fit_dist=="lnorm")casenum=casenum+4					
if(dist=="gumbel") casenum=casenum+8					
					
	if(casenum < 4) {				
		if(npar == 2)  {			
			if(abpval == TRUE) {		
				prr<-AbPval(dim(x)[1], resultVec[3])	
				outVec<-c(Eta=resultVec[1],Beta=resultVec[2],Rsqr=resultVec[3], AbPval=prr[1])	
				names(outVec)<-c("Eta", "Beta", "Rsqr", "AbPval")	
			}else{		
				## calling with abpval=FALSE returns an unnamed vector	
				outVec<-c(resultVec[1],resultVec[2],resultVec[3])	
			}   		
					
					
		}else{			
			outVec<-c(Eta=resultVec[1],Beta=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])		
			if(resultVec[5]==1)  {		
				warn="3p optimization did not converge"	
				attr(outVec,"warning")<-warn	
			}		
		}			
	}else{				
		if(casenum < 8) {			
			if(npar == 2)  {		
				if(abpval == TRUE) {	
					prr<-AbPval(length(x[,1]), resultVec[3],"lnorm")
					outVec<-c(Mulog=resultVec[1],Sigmalog=resultVec[2],Rsqr=resultVec[3], AbPval=prr[1])
					names(outVec)<-c("Mulog", "Sigmalog", "Rsqr", "AbPval")
				}else{	
					## calling with abpval=FALSE returns an unnamed vector
					outVec<-c(resultVec[1],resultVec[2],resultVec[3])
				}   	
			}else{		
				outVec<-c(Mulog=resultVec[1],Sigmalog=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])	
				if(resultVec[5]==1)  {	
					warn="3p optimization did not converge"
					attr(outVec,"warning")<-warn
				}	
			}		
##                gumbel support is currently suspended					
##		}else{			
##			if(length(resultVec)==3)  {		
##				outVec<-c(Etalog=resultVec[1],Betalog=resultVec[2],Rsqr=resultVec[3])	
##			}else{		
##				outVec<-c(Etalog=resultVec[1],Betalog=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])	
##				if(resultVec[5]==1)  {	
##					warn="3p optimization did not converge"
##					attr(outVec,"warning")<-warn
##				}	
##			}		
		}			
	}				
					
return(outVec)					
}					
