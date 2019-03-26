## p2y.r
##
## This function originally appeared as F0inv in the abrem package by Jurgen Symynck
## That function used the graphics log option setting to distinguish between weibull 
## cand lognormal anvas selections. That limited its usefulness outside of the S3 plot 
## function since the log argument had an unclear meaning.
## 
## Now with a canvas argument, this makes sense as an exported function. It has proven
## to be very handy when adding points and lines to plots generated from plot.wblr

## This function is equivalent to SPLEDA::qsev when applying the default weibull canvas 

p2y <- function(p,canvas="weibull"){
    # This is the inverse Cumulative Distribution function
	# used to transform a probability value to the
    # y-axis of the plot canvas. 
	# Use of this transformation permits distributions
	# to appear as curves on unrelated canvas
    if(canvas =="weibull")ret <- log(qweibull(p,1,1))
	if(canvas =="lognormal") ret <- qlnorm(p,0,1)
    ret
}