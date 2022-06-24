xbounds<-function(obj, fit_num=1, conf_num=1)  {		
	if(!is(x, "wblr")) stop('Argument \"x\" is not of class \"wblr\" ')	
	if(is.null(obj$fit)) {	
		stop("no fit found in wblr object")
	}else{	
		if(length(obj$fit)<fit_num) {
		stop(paste0("fit_num [[",fit_num,"]] not found in wblr object"))
		}
	}	
		
		
	if(is.null(obj$fit[[fit_num]]$conf)) {	
		stop("no bounds found in wblr object")
	}else{	
		if(length(obj$fit[[fit_num]]$conf)<conf_num) {
		stop(paste0("conf_num [[",conf_num,"]] not found in wblr object"))
		}
	}	
		
	out_obj<-obj$fit[[fit_num]]$conf[[conf_num]]$bounds	
		
	if(obj$fit[[fit_num]]$modified == TRUE)  {	
		return(list(bounds=out_obj, modified.by.t0=obj$fit[[fit_num]]$t0))
		
	}else{	
		return(list(bounds=out_obj))
	}	
}		
