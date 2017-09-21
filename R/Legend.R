legendConf <- function(fit,conftype,opadata,...){

    if(!is.null(fit)){

        if(!is.null(fit$options))
            opafit <- modifyList(opadata,fit$options)
        opafit <- modifyList(opafit,list(...))
        if(identical(tolower(conftype),"blives")){
            if(!is.null(fit$conf$blives)){
                for.each.blicon <- function(blicon){
                    if(!is.null(blicon$options)){
                        opaconf <- modifyList(opafit,blicon$options)
                    }else{opaconf <- opafit}
                    if(opaconf$in.legend){
                            # TODO: correct usage of this logical value?
                        li <- list()
                        li[[1]] <-  bsll(legend=paste0("CI bounds, type = ",
                            ifelse(is.null(blicon$type),"NA",
                            paste0("\"",blicon$type,"\""))),
                            col=opaconf$col,lwd=opaconf$lwd,lty=opaconf$lty)
                        li[[2]] <- bsll(legend=paste0("  CI = ",
                            ifelse(is.null(blicon$cl),"NA",
                                paste0(signif(blicon$cl*100,4)," [%]")),
                            ifelse(is.null(blicon$S),"",
                                paste0(", S = ",blicon$S))))
                        if(opaconf$in.legend.blives){
                            params <- unlist(list(beta=fit$beta,eta=fit$eta,t0=fit$t0,
                                meanlog=fit$meanlog,sdlog=fit$sdlog,rate=fit$rate))
                            if(is.null(bl <- blicon$unrel))bl <- opaconf$unrel
                            fu <- function(bl){
                                bsll(legend=Blifestring(bl,blicon,opafit$signif,params))
                            }
                            c(li,lapply(bl,fu))
                        }else(li)
                    }else NULL
                }
                #mtrace(for.each.blicon)
                unlist(lapply(fit$conf$blives,for.each.blicon),FALSE)
                    # TODO: replace by do.call ?
            }else{NULL}
        }
    }else{NULL}
}

buildSingleDataLegend <- function(x,opadata,...){
    arg <- list(...)
    si <- function(number)signif(number,opadata$signif)
    li <- list()
    if(opadata$in.legend){

        li[[10]]    <- bsll(legend=paste0("ranks = ",opadata$pp[1]),
            col=opadata$col,pch=opadata$pch,lwd=opadata$lwd.points)
		if(is.null(x$data$dlines)) {
			li[[15]]    <- bsll(legend=paste0("n (fail | cens.) = ",x$n,
				" (",x$fail," | ",x$cens,")"))
		}else{
			li[[15]]    <- bsll(legend=paste0("n (f | s | d | i) = ",x$n,
				" (",x$fail," | ",x$cens," | ",x$discovery," | ",x$interval,")"))		
		}
    }
    removeBadLegendEntries <- function(e){
        if(!is.null(e))!is.na(e$legend) else FALSE
    }
    if(length(li)>0)li <- li[sapply(li,removeBadLegendEntries)]
    else li <- ""
        # remove list items where the legend text = NA
    fu  <- function(x,i){if(i %in% names(x))x[[i]]}
    fu2 <- function(i,x){lapply(x,fu,i=i)}
    items <- c("legend","lty","lwd","pch","col")
    le  <- lapply(items,fu2,li)
    names(le) <- items
    if(identical(label <- opadata$label,""))label <- NULL
    le$rect <- legend(
        "bottomright",
        legend=le$legend,
        title=label,
        cex = opadata$legend.text.size,
        plot=FALSE)$rect
    le$label <- opadata$label
    le$legend.text.size <- opadata$legend.text.size
    le
}

buildSingleFitLegend <- function(fit,opadata,...){
    arg <- list(...)
    if(!is.null(fit$options)){
        opafit <- modifyList(opadata,fit$options)
    }else{opafit <- opadata}
    opafit <- modifyList(opafit,list(...))
    t0 <- NULL
    le <- NULL
    
    if(opafit$is.plot.legend){
	
## removing threshold effects, when needed t0 is fit$t0	 fit$modified indicates data modification by t0
        #if(is.logical(opafit$threshold))if(opafit$threshold){
        #   if(is.logical(opadata$threshold)){if(opadata$threshold)
        #       warning("opafit$threshold and opadata$threshold are logical values but numeric values were expected. Proceeding...")
        #    }else{
                # reuse the t0 value from the data level
        #       t0 <- opadata$threshold
        #    }
        #}
        #if(is.numeric(opafit$threshold))t0 <- opafit$threshold
		
		
        si <- function(number)signif(number,opafit$signif)
            # shorter writing form for signif()
        li <- list()
        if(opadata$in.legend){
    		modstr<-""
		if(fit$modified) modstr<- "*t0 mod* "
            li[[10]]    <- bsll(legend=paste0(modstr,"ranks = ",opafit$pp[1]),
                col=opadata$col,pch=opadata$pch,lwd=opadata$lwd.points)
		if((fit$discovery+fit$interval)==0) {
			li[[15]]    <- bsll(legend=paste0("n (fail | cens.) = ",fit$n,
				" (",fit$fail," | ",fit$cens,")"))
		}else{
			li[[15]]    <- bsll(legend=paste0("n (f | s | d | i) = ",fit$n,
				" (",fit$fail," | ",fit$cens," | ",fit$discovery," | ",fit$interval,")"))		
		}
        }
        if(opafit$in.legend){
            li[[20]]    <- bsll(legend = paste0(fit$options$dist," (",
                paste0(fit$options$method.fit,collapse=", "),")"),
                col=opafit$col,lwd=opafit$lwd,lty=opafit$lty)
            li[[30]]    <- bsll(legend=ifelse(is.null(fit$rate),NA,
                    paste0("rate = ",si(fit$rate))))
            li[[40]]    <- bsll(legend=ifelse(is.null(fit$meanlog),NA,
                    paste0("mean(log) = ",si(exp(fit$meanlog))," (",
                    si(fit$meanlog),")")))
            li[[50]]    <- bsll(legend=ifelse(is.null(fit$sdlog),NA,
                    paste0("sd(log) = ",si(exp(fit$sdlog))," (",
                    si(fit$sdlog),")")))
            li[[60]]    <- bsll(legend=ifelse(is.null(fit$beta),NA,
                    paste0("beta = ",si(fit$beta))))
            li[[70]]    <- bsll(legend=ifelse(is.null(fit$eta),NA,
                    paste0("eta = ",si(fit$eta))))
            li[[80]]    <- bsll(legend=ifelse(is.null(fit$t0),NA,
                    paste0("t0 = ",si(fit$t0))))
            if(!is.null(fit$gof) && opafit$in.legend.gof){
                if(!is.null(fit$gof$r2)){
                    if(!is.null(fit$gof$ccc2)){
                        li[[100]]    <- bsll(legend=paste0("r^2 | CCC^2 = ",
                            si(fit$gof$r2)," | ",si(fit$gof$ccc2),
                            ifelse(fit$gof$r2>=fit$gof$ccc2," (good)"," (BAD)")))
                    }else{
                        li[[100]]    <- bsll(legend=paste0("r^2 = ",si(fit$gof$r2)))
                    }
                }
                if(!is.null(fit$gof$loglik)){
                    li[[110]]    <- bsll(legend=paste0("loglik = ",si(fit$gof$loglik)))
                }
                li[[120]]    <- bsll(
                    legend=ifelse(is.null(fit$gof$prr),NA,
                        #paste0("prr = ",si(fit$gof$prr)," (S=",
                        #ifelse(is.null(fit$gof$S),"NA",fit$gof$S),")")))
						paste0("prr = ",si(fit$gof$prr)," by corr.")
					)
				)
            }
        }
        #leconfpos <- length(na.omit(unlist(li))) + 1
            # where displaying confidence info begins
        leconf <- legendConf(fit,"blives",opadata=opadata,...)
        if(!is.null(leconf))li[[130]] <- bsll(legend="")
        li <- c(li,leconf)
        removeBadLegendEntries <- function(e){
            if(!is.null(e))!is.na(e$legend) else FALSE
        }
        if(length(li)>0)li <- li[sapply(li,removeBadLegendEntries)]
        else li <- ""
            # remove list items where the legend text = NA
        fu  <- function(x,i){if(i %in% names(x))x[[i]]}
        fu2 <- function(i,x){lapply(x,fu,i=i)}
        items <- c("legend","lty","lwd","pch","col")
        le  <- lapply(items,fu2,li)
        names(le) <- items
        if(identical(label <- opafit$label,""))label <- NULL
        le$rect <- legend(
            "bottomright",
    #                "topright",
            legend=le$legend,
            title=label,
            cex = opafit$legend.text.size,
    #        inset=0.1,
    #        merge = TRUE,
            plot=FALSE)$rect
        le$label <- opafit$label
        le$legend.text.size <- opafit$legend.text.size
    }
    le
}

bsll <- function(...){
    arg <- list(...)
    leline <- list(
        legend= NA,
        lty= NA,
        lwd= NA,
        pch= NA,
        col= NA)
    modifyList(leline,arg)
#    leline <- list(
#        legend= <- ifelse(is.null(arg$legend),NA,arg$legend)
#        title= <- ifelse(is.null(arg$title),NA,arg$title)
#        cex= <- ifelse(is.null(arg$cex),NA,arg$cex)
#        bg= <- ifelse(is.null(arg$bg),NA,arg$bg)
#        lty= <- ifelse(is.null(arg$lty),NA,arg$lty)
#        lwd= <- ifelse(is.null(arg$lwd),NA,arg$lwd)
#        pch= <- ifelse(is.null(arg$pch),NA,arg$pch)
#        col= <- ifelse(is.null(arg$col),NA,arg$col)
}

Blifestring <- function(B,blicon,signif,...){
    # This functions creates a string for displaying the B-lives in the plot's
    # legend. missing input data result in an "NA". For example, the output
    # string could look like:
    #   "B10 = 9.86 | 50.13 | 103.4"
    # or
    #   "B1 = 9.86 | 50.13 | NA"
    si <- function(number)
        if(!is.null(number))signif(number,signif)
        else NA
      # shorthand writing of the signif() function
    qfun <- function(B,...){
        args <- as.list(unlist(...))

        ret <- NULL
        if(!is.null(args$beta) && !is.null(args$eta)){
            # the fit type was weibull
            ret <- qweibull(B,args$beta,args$eta)
            if(!is.null(args$t0)){
                # the fit type was weibull3p
                ret <- ret+args$t0

            }
        }
        if(!is.null(args$meanlog) && !is.null(args$sdlog)){
            # the fit type was lognormal
            ret <- qlnorm(B,args$meanlog,args$sdlog)
        }
        if(!is.null(args$rate)){
            # the fit type was exponential
            ret <- qexp(B,args$rate)
        }
        ret
    }
    id <- function(x,y)isTRUE(all.equal(x,y))
    c1 <- is.null(blicon$bounds) || is.null(blicon$bounds$Lower)
    if(!c1) lo <- si(subset(blicon$bounds,
        sapply(blicon$bounds$unrel,id,B),Lower))
    c2 <- is.null(blicon$bounds) || is.null(blicon$bounds$Upper)
    if(!c2) up <- si(subset(blicon$bounds,
        sapply(blicon$bounds$unrel,id,B),Upper))
    ret <- paste(sep = "","    B",signif(100*B)," = ",
        ifelse(c1,
           "NA",lo),
        " | ",si(qfun(B,...)),
        " | ",ifelse(c2,
           "NA",up))
    ret
}
