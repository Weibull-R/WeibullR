## debias.r
## This file includes a selection of functions for the adjustment of bias in the shape or spread parameters
## of Weibull or normal/log-normal distributions when performed on small data samples by the MLE method.  
## 
## (C) OpenReliability.org 2013-2017
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
## GNU General Public License for more details.


## RBA
## The RBA (Reduced Bias Adjustment) has been promoted by Dr Abernethy in "The New Weibull Handbook" for several decades.
## This concept is not new as Dr. Abernethy cites work on the underlying 'C4' term by Gosset on the normal distribution
## in the early 20th century. The primary difference noted by Dr. Abernethy is on the choice of the median bias reduction
## for the Weibull shape parameter, rather than the mean reduction, which is an alternative choice.  This difference does
## not exist for adjustment of the standard deviation of the normal, or log standard deviation of log-normal.
## As Dr. Abernethy explains, for the symmetrical distributions the mean and median bias are the same.

rba<-function(Qx, dist="weibull", basis="median")  {
	if(Qx<3)  {
		stop("insufficient data points")
	}
	if(Qx>343) return(1.0)

## factorial(x) is simply gamma(1+x)
	num<-gamma(1+(Qx-2)/2)	
	den<-gamma(1+(Qx-3)/2)	
	C4<-sqrt(2/(Qx-1))*num/den
	
	if(dist=="weibull")  {
		if(basis=="median")  {
			return(C4^3.5)	
		}
		if(basis=="mean")  {
			return(C4^6)
		}
	}
	if(tolower(dist) %in% c("lnorm","lognormal","lognormal2p"))  {
		return(sqrt(Qx/(Qx-1))/C4)
	}		
}	

## HRBU
## This is a combined implementation of the Hirose and Ross Beta-Unbias adjustments for Weibull MLE on small samples.
## These functions are used by ReliaSoft's Weibull++ product as presented in ReliaSoft Corporation,
## Life Data Analysis Reference, Tucson, AZ: ReliaSoft Publishing, 2005.  It can be observed that the Hirose
## correlation is nearly identical to the "mean" adjustment characterized by C4^6.  The Ross correlation adds a modest
## increase to the adjustment with increasing suspension, or right-censored, data.

hrbu<-function(Qx,Qs=NULL)  {
	r<-Qx
	n<-r+Qs
    if (length(Qs)==0||Qs==0) {
## This is the Hirose Beta Unbias factor for complete failure samples	
	BU<-1/( 1.00115+1.278/r+2.001/r^2+20.35/r^3-46.98/r^4 )
	}else{
## This is the Ross Beta Unbias factor for samples with suspensions	
	BU<- 1/(1+1.37/(r-1.92)*sqrt(n/r))
	}
BU	
}