## MRRw2p.r
## A quick fit method using defaults, with options to provide confidence interval bounds
## and/or a graphical display.

MRRw2p<-function(x,s=NULL, bounds=FALSE, show=FALSE) {

# permit one to forget that second argument was reserved for a suspensions vector
if(is.logical(s)) {
	show<-bounds
	bounds<-s
	s<-NULL
}

obj<-wblr.fit(wblr(x,s), col="blue2")
fit<-obj$fit[[1]]$fit_vec
attributes(fit)$data_types<-NULL
## bounds are not prepared for 3p fits, so the bounds argument is simply ignored.
if(bounds==TRUE) {
obj<-wblr.conf(obj, dq="minitab", col="deepskyblue3")
bnds<-obj$fit[[1]]$conf[[1]]$bounds
ret<-list(fit,bnds)
stitle<-"MRR fit with 90% double-sided pivotal bounds"
}else{
ret<-fit
stitle<-"Lognormal RR fit"
}

if(show==TRUE) {
plot(obj, in.legend.blives=FALSE, sub=stitle)
}

ret
}