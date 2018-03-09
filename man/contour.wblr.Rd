\name{contour.wblr}
\alias{contour.wblr}

\title{
    S3 Likelihood Ratio Contour Maps From \code{wblr} Objects
}
\description{
    This function adds the \code{.wblr} method to \code{\link[graphics]{contour}}
    from the \pkg{graphics} package.
}
\usage{\method{contour}{wblr}(x, \dots)}
\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of \code{wblr} objects.
    }
    \item{\dots}{
        A placeholder required for S3 method conistency intended for setting graphic parameters
        (like col, lty, lwd) but NOT IMPLEMENTED. Such graphic parameters are to be added in the
        \code{.wblr} objects themselves.
    }
}
\details{
    This function provides S3 object functionality for plotting a likelihood ratio contour map for
    any single wblr object with just the contour function.
 
    Contour parameters are drawn from a contour existing in the object(s) passed in, or 
    from the base object (sometimes just defaults) if no contour exists in the object.
 
    In order to plot contour maps from multiple objects onto a single canvas
    it is necessary to call contour.wblr specifically with a list of wblr objects as primary argument.

    Unlike the \code{\link[graphics]{contour}} this map is not prepared from a matrix of z-values.
    Rather the specific x,y points for each CL level (Z-value) are provided for plotting lines
    connecting the points.
	
}
\value{
    The contour.wblr function itself returns no value, however, the full output of contour points and
    parameters from the internally called WeibullR::plot_contour function is displayed.
}

\examples{
set.seed(1234)
da2 <- wblr.conf(wblr.fit(wblr(runif(5,100,1000),col="red")))
da3 <- wblr.conf(wblr.fit(wblr(rweibull(5,3,1000),col="green4")))

% not run due to device dependence
\dontrun{
contour.wblr(list(da2,da3))
}
}
