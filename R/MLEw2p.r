## MLEw2p.r
## A quick fit method using defaults, with options to provide confidence interval bounds
## and/or a graphical display.

MLEw2p<-function(x,s=NULL, bounds=FALSE, show=FALSE) {

# permit one to forget that second argument was reserved for a suspensions vector
if(is.logical(s)) {
	show<-bounds
	bounds<-s
	s<-NULL
}

obj<-wblr.fit(wblr(x,s), method.fit="mle", col="darkgreen")
fit<-obj$fit[[1]]$fit_vec
attributes(fit)$data_types<-NULL
if(bounds==TRUE) {
obj<-wblr.conf(obj, method.conf="lrb", dq="minitab", col="green3")
bnds<-obj$fit[[1]]$conf[[1]]$bounds
ret<-list(fit,bnds)
stitle<-"MLE fit with 90% double-sided likelihood ratio bounds"
}else{
ret<-fit
stitle<-"Weibull MLE fit"
}

if(show==TRUE) {
plot(obj, in.legend.blives=FALSE, sub=stitle)
}

ret
}