# Extract a Fit Summary from a 'wblr' Object

xfit<-function(obj, fit_num=1)  {			
	if(!is(obj, "wblr")) stop('Argument \"obj\" is not of class \"wblr\" ')
	if(is.null(obj$fit)) {
		stop("no fit found in wblr object")
	}else{
		if(length(obj$fit)<fit_num) {
		stop(paste0("fit_num [[",fit_num,"]] not found in wblr object"))
		}
	}
	
	dist=obj$fit[[fit_num]]$options$dist		
	fit=obj$fit[[fit_num]]$fit_vec		
	if(is.null(attr(fit,"data_types")) ) {		
		data_types<-c(	
		obj$fit[[fit_num]]$n,	
		obj$fit[[fit_num]]$cens,	
		obj$fit[[fit_num]]$discovery,	
		obj$fit[[fit_num]]$interval)
		names(data_types)<-c("n","s","d","i")	
		attr(fit,"data_types")<-data_types	
	}		
	return(list(dist=dist,fit=fit))		
}