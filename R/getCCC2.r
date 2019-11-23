getCCC2<-function(F, model="weibull")  {

# getCCC2_r.r
# Correlation developed by David Silkworth in Oct 2013
# Back-converted in Nov 2019 from ansi c code in abpv.c from abremPivotals latest modification (May 2014)

getCCC2_r<-function(nF, dist="weibull")  {							
	T1w<-c(0.7924059,0.7990601,0.8076817,0.8205148,0.8324683,0.842629,0.8515389,0.8593485,0.8662833,						
		0.8724438,0.8779514,0.8828717,0.8873306,0.8913891,0.89512,0.8986051,0.9017277,0.9045952,					
		0.90739,0.9098512,0.9122063,0.9144563,0.9165037)					
	T2w<-c(2.482953,2.592705,2.687305,2.772759,2.84817,2.916081,2.978664,3.033237,3.085852,3.136547,						
		3.18183,3.226624,3.267314,3.305018,3.3424,3.377739,3.412428,3.444103,3.474612,3.504789,					
		3.531535,3.559953,3.587629,3.61134,3.637438,3.661817)					
	T3w<-c(3.661817,3.870633,4.035835,4.172313,4.290143,4.392362,4.484968,4.565962,4.641202,4.710528,						
		4.774699,4.83303,4.888764,4.940351,4.989649,5.034619,5.079213,5.120318,5.160426,5.200657,					
		5.235989,5.270168,5.304064,5.336199,5.3666,5.396415)					
							
	T1l<-c(0.7938923,0.7992166,0.8143357,0.8286594,0.8416131,0.8531055,0.863076,0.8717764,0.8794219,						
		 0.8862083,0.8921895,0.8975986,0.9024265,0.9068011,0.9107908,0.9144347,0.9177708,0.9208458,					
		 0.9236726,0.9262948,0.9287454,0.931017,0.9331573)					
	T2l<-c(2.705413,2.847212,2.969813,3.077389,3.173427,3.260117,3.339296,3.412094,3.479407,3.542081,						
		3.600777,3.655789,3.707801,3.756996,3.803559,3.847988,3.890183,3.93063,3.969467,4.006691,					
		4.042462,4.076862,4.109992,4.142034,4.173005,4.202877)					
	T3l<-c(4.202877,4.458735,4.659025,4.823861,4.963904,5.085855,5.193705,5.290567,5.378335,5.45858,						
		5.532547,5.601127,5.6651,5.725016,5.781398,5.834598,5.884945,5.932743,5.978272,6.0218,					
		6.06339,6.103111,6.141497,6.178083,6.213509,6.2477,6.280481,6.312331)					
							
						
	if(dist=="weibull")  {						
		T1<-T1w					
		T2<-T2w					
		T3<-T3w					
	}else{						
		T1<-T1l					
		T2<-T2l					
		T3<-T3l					
	}						
	if(nF<26)  {						
		CCC2<-T1[nF-2]					
		}else{					
			if(nF<151)  {				
				i=5			
				if(nF%%i==0)  {			
	## The qweibull of CCC2 can be taken directly from T2						
	## Offset value is 25/i-1 (for R) Will be 25/i for C++						
					CCC2<-1-1/exp(T2[nF/i-4])		
				}else{			
	## The qweibull of CCC2 will have to be interpolated from T2						
	## establish  nF and qweibull bounds						
					nFbl<-i*as.integer(nF/i)		
					nFbu<-nFbl+i		
					qwl<-T2[nFbl/i-4]		
					qwu<-T2[nFbu/i-4]		
	## Then interpolate using log(F) and log(Fbounds)						
					qwccc2<-qwl+((log(nF)-log(nFbl))/(log(nFbu)-log(nFbl))*(qwu-qwl))		
					CCC2<-1-1/exp(qwccc2)		
				}			
			}else{				
				if(nF<1401)   {			
					i=50		
					if(nF%%i==0)  {		
	## The qweibull of CCC2 can be taken directly from T3						
	## Note there is a difference in the F/i offset for element selection!!!						
	## In this case the offset = 150/i-1 (for R), will be 150/i for C++						
						CCC2<-1-1/exp(T3[nF/i-2])	
					}else{		
	## The qweibull of CCC2 will have to be interpolated from T3						
							
	## establish  nF and qweibull bounds						
						nFbl<-i*as.integer(nF/i)	
						nFbu<-nFbl+i	
						qwl<-T3[nFbl/i-2]	
						qwu<-T3[nFbu/i-2]	
	## Then interpolate using log(nF) and log(nFbounds)						
						qwccc2<-qwl+((log(nF)-log(nFbl))/(log(nFbu)-log(nFbl))*(qwu-qwl))	
						CCC2<-1-1/exp(qwccc2)	
					}							
				}else{			
					warning(paste0("Quantity ",nF," failures has not been correlated to CCC2"))	
					CCC2<-NA
				}			
			}				
		}					
	CCC2						
}							
	

return(getCCC2_r(F, model))
}								
