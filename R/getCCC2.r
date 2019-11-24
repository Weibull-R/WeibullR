getCCC2<-function(F, model="weibull")  {

# getCCC2_r.r
# Correlation developed by David Silkworth in Oct 2013
# Back-converted in Nov 2019 from ansi c code in abpv.c from abremPivotals latest modification (May 2014)

getCCC2_r<-function(nF, dist="weibull")  {							
	T1w<-c(	0.792235,0.7990604,0.8076126,0.8204102,0.832331,0.8425375,0.8514909,0.8593213,0.8662665,
		 0.8724075,0.8779149,0.8828711,0.887337,0.8914107,0.8951199,0.8985272,0.9016913,0.9045873,     
		 0.9073225,0.9098415,0.9121919,0.9144172,0.9164957)    
	T2w<-c(	2.482857,2.593721,2.689116,2.773029,2.847969,2.915728,2.977454,3.034384,3.087205,3.136568,
		 3.182761,3.226223,3.267092,3.306225,3.343156,3.378561,3.412309,3.444714,3.475747,3.505642,     
		 3.534433,3.562236,3.588777,3.614618,3.639701,3.663851)   
	T3w<-c(	3.663851,3.872307,4.037513,4.17472,4.292113,4.394901,4.486445,4.568626,4.643578,4.712588,
		 4.776053,4.83515,4.890495,4.942337,4.991377,5.037647,5.081707,5.12344,5.163087,5.2011,     
		 5.237467,5.272338,5.305946,5.338348,5.369513,5.399487,5.428722,5.456727)					
							
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
