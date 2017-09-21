Aitken<-function(X,Y,x)  {			
		n=length(X)	
		i=1	
	while(i <n+1)  {			
			j=1
		while( (j+i)< n+1 )  {
				Y[j]= (x-X[j])/(X[j+i]-X[j])*Y[j+1]+(x-X[j+i])/(X[j]-X[j+i])*Y[j]			
			j=j+1
		}	
		i=i+1	
	}			
	Y[1]			
}			
