## contour.wblr.R
## inspired by code originally authored by Jurgen Symynck, April 2014
## Copyright 2014-2018 OpenReliability.org

# contour.wblr provides S3 object functionality for plotting a likelihood ratio contour map for
# any single wblr object with just the contour function.
# In order to plot contour maps from multiple objects in a single canvas
# it is necessary to call contour.wblr specifically with a list of wblr objects as primary argument.
#
# contour parameters are drawn from a contour existing in the object(s) passed in, or 
# from the base object (sometimes just defaults) if no contour exists in the object.

contour.wblr<-function(x,...)  {

plot_contour(x,CL=seq(.1,.9,by=.1))


}