## getPercentilePlottingPositions.r file
 ##
 ## Author: David J. Silkworth
 ##   (c)2014-2017 OpenReliability.org
##

getPercentilePlottingPositions<-function(x, s=NULL, interval=NULL, ppos="Benard", aranks="Johnson", ties=NULL)  {							
	F<-length(x)						
	N<-F+length(s)						
	##  create the event vector						
	if(!missing(s)) {						
	## if(length(s)>0)  {						
		## suspension data has been provided					
		    data<-c(x,s)					
		    event<-c(rep(1,F),rep(0,N-F))					
		    prep_df<-data.frame(data=data,event=event)					
		## now sort the dataframe on data values					
		    NDX<-order(prep_df[,1])					
		    prep_df<-prep_df[NDX,]					
	  }else{						
		## this is simply a complete failure set					
		    data<-sort(x)					
		    event<-rep(1,F)					
		    prep_df<-data.frame(data=data,event=event)					
	  }						
							
	if(aranks=="Johnson")  {						
		## adjust ranks using Drew Auth's simplification of Leonard Johnson's method					
		## start with extra element to reference zero as previous rank to first					
		adj_rank<-0					
		for(k in 1:N)  {					
			rr=N-k+1				
			if(prep_df$event[k]>0)  {				
			this_rank<-(rr*adj_rank[k]+N+1)/(rr+1)				
			adj_rank<-c(adj_rank, this_rank)				
			}else{				
			adj_rank<-c(adj_rank,adj_rank[k])				
			}				
		}					
		prep_df<-cbind(prep_df,adj_rank=adj_rank[-1])					
		## now eliminate the suspension data					
		prep_df<-prep_df[sapply(prep_df$event, function(x) x>0),c(1,3)]					
	}else{						
		if(aranks=="KMestimator")  {					
		## adjust ranks using David Silkworth's adaptation of the modified					
		## Kaplan-Meier estimator used by Minitab as "nonparametric"					
			## start with extra element to reference zero as previous rank to first				
			adj_rank<-0				
			for(k in 1:N)  {				
				if(prep_df$event[k]>0)  {			
					this_rank<-1-((1-adj_rank[k])*(N-k)/(N-k+1))		
				adj_rank<-c(adj_rank, this_rank)			
				}else{			
				adj_rank<-c(adj_rank,adj_rank[k])			
				}			
			}				
			prep_df<-cbind(prep_df,adj_rank=adj_rank[-1])				
			## now eliminate the suspension data				
			prep_df<-prep_df[sapply(prep_df$event, function(x) x>0),c(1,3)]				
			## Now provide a modification for the final element if it was a failure
			## This adjustment is used by Minitab
			if(prep_df$adj_rank[F]==1)  {				
				prep_df$adj_rank[F]=1-((1-prep_df$adj_rank[F-1])*1/10)			
			}				
			## Finally reverse the Kaplan-Meier plotting position to reveal useable
			## adj_rank	for any plotting position method.			
			prep_df$adj_rank<-prep_df$adj_rank*N
		}else{	
			stop("aranks argument not recognized")		
		}					
	}						
							
	## now to handle ties, if called for						
	if(!missing(ties))  {						
	## if(length(ties)>0) {						
		test_hi<-prep_df$data-c(prep_df$data[-1],0)					
		test_lo<-prep_df$data-c(1,prep_df$data[1:(F-1)])					
		highest<-prep_df[sapply(test_hi, function(y) y!=0),]					
		lowest<-prep_df[sapply(test_lo, function(y) y!=0),]					
		if(ties=="highest")  {					
			prep_df<-highest				
		}else{					
			if(ties=="lowest")  {				
				prep_df<-lowest			
			}else{				
				if(ties=="mean")  {			
					prep_df<-data.frame(data=highest$data,adj_rank=(highest$adj_rank+lowest$adj_rank)/2)		
				}else{			
					if(ties=="sequential")  {		
						seq_adj<-cumsum(highest$adj_rank-lowest$adj_rank)	
						prep_df<-data.frame(data=highest$data, adj_rank=highest$adj_rank-seq_adj)	
					}else{		
						stop("ties argument not recognized")	
					}		
				}			
			}				
		}					
	}						
							
	## special note:  the number of data entries N remains as originally identified						
	## if CANNOT be changed due to removal of suspensions or ties						
							
	## prep_df now contains the data to be plotted with their final adjusted ranks						
	## finally we get the probability plotting positions						
							
	if(ppos=="Benard"||ppos=="beta"||ppos=="mean")  {						
		if(ppos=="Benard")  {					
			ppp<-(prep_df$adj_rank-0.3)/(N+0.4)				
		}else{					
			## the incomplete beta function is what				
			## Benard was approximating				
			## we can simply use this directly				
			if(ppos=="beta")  {				
				ppp<-qbeta(0.5, prep_df$adj_rank,N-prep_df$adj_rank+1)			
			}else{				
				if(ppos=="mean")  {			
					## mean ranks (aka Herd-Johnson) give even spacing 		
					## mean ranks were originally proposed by Weibull, but later		
					## he agreed that the Benard approximation gave		
					## more pleasing results when used with fatigue failure data		
					ppp<-prep_df$adj_rank/(N+1)		
				}			
			}				
		}					
	}else{						
		if(ppos=="Hazen"||ppos=="Kaplan-Meier"||ppos=="Blom")  {					
			if(ppos=="Hazen")  {				
				## aka modified Kaplan-Meier			
				ppp<-(prep_df$adj_rank-0.5)/N			
			}else{				
				if(ppos=="Kaplan-Meier")  {			
				## This is method used by SuperSMITH with Johnson aranks			
					ppp<-(prep_df$adj_rank)/N		
					K<-length(prep_df[,1])		
					if(prep_df$adj_rank[K]==N)  {		
						ppp[K]<-K/(N+0.001)	
					}		
				}else{			
					if(ppos=="Blom")  {		
						ppp<-(prep_df$adj_rank-0.375)/(N+0.25)	
					}		
				}			
			}				
		}else{
			stop("ppos argument not recognized")
		}		
	}						
							
	outDF<-cbind(prep_df$data,data.frame(ppp),prep_df$adj_rank)						
	colnames(outDF)<-c("time","ppp","adj_rank")						
return(outDF)							
}							
## assign the alias							
getPPP<-getPercentilePlottingPositions							
