## MRRw3p.r
## A quick fit method using defaults, with options to provide confidence interval bounds
## and/or a graphical display.

MRRw3p<-function(x,s=NULL, bounds=FALSE, show=FALSE) {

# permit one to forget that second argument was reserved for a suspensions vector
if(is.logical(s)) {
	bounds<-s
	show<-bounds
	s<-NULL
}

obj<-wblr.fit(wblr(x,s), dist="weibull3p", modify.by.t0=T, col="blue2")
fit<-obj$fit[[1]]$fit_vec
attributes(fit)$data_types<-NULL
## bounds are not prepared for 3p fits, so the bounds argument is simply ignored.
#if(bounds==TRUE) {
#obj<-wblr.conf(obj, method.conf="lrb", dq="minitab", col="deepskyblue3")
#bnds<-obj$fit[[1]]$conf[[1]]$bounds
#ret<-list(fit,bnds)
#stitle<-"MLE fit with 90% double-sided likelihood ratio bounds"
#}else{
ret<-fit
stitle<-"Weibull MRR 3-parameter fit"
#}

if(show==TRUE) {
plot(obj, in.legend.blives=FALSE, sub=stitle)
}

ret
}