\name{plot_contour}
\alias{plot_contour}

\title{
    Plotting of Likelihood Ratio Contours from \code{wblr} Object Data
}

\description{
    Currently, this function plots likelihood ratio contours for weibull distribution fitted data, only, as found in an \code{wblr}
    object or a list of \code{wblr} objects. In the future this will be expanded to cover lognormal fitted data also with options to calculate modified
	versions of the contour based on alteration of the fit distribution, confidence level, and degrees of freedom based solely on the data in the listed \code{wblr} object.
}

\usage{plot_contour(x, CL=NULL)}

\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of \code{"wblr"} objects.
    }
    \item{CL}{An optional argument for establishing the confidence limit(s) for calculated contours. A vector of values may be provided to generate concentric contours on a single object.
    }
}

\value{
    This function returns the list of contours that were prepared for plotting.
}

\details{
   If no CL argument is provided, the contours generated in the provided \code{wblr} objects for likelihood confidence bounds will be plotted.
   If a vector is provided for CL, multiple contours for each provided object will be calculated to display concentric contours according to values in the CL vector.
   In all cases options for distribution and degrees of freedom, as well as graphic options will be drawn from the underlying objects.
}

\examples{
set.seed(1234)
fail<-rweibull(5,2,10)
ctest<-wblr(fail)
ctest<-wblr.fit(ctest, method.fit="mle")
ctest<-wblr.conf(ctest, method.conf="lrb",col="red")
fail2<-c(10,40,40,50)
susp<-c(20,60)
ctest2<-wblr(fail2, susp)
ctest2<-wblr.fit(ctest2, method.fit="mle")
ctest2<-wblr.conf(ctest2, method.conf="lrb",col="blue")
% not run due to device dependence
% plot_contour(list(ctest, ctest2))
}