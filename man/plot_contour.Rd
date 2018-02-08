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

\usage{plot_contour(x,dist=NULL, CL=NULL, dof=NULL, \dots)}

\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of \code{"wblr"} objects.
    }
    \item{dist}{An optional argument for establishing the distribution (only "weibull" and "lognormal" are accepted) for a calculated contour. 
    }
    \item{CL}{An optional argument for establishing the confidence limit(s) for calculated contours. A vector of values may be provided to generate concentric contours on a single object.
    }
    \item{dof}{An optional argument for establishing the degrees of freedom for calculated contours. It is most accurate to use dof=2 for dataset comparisons.
    }
    \item{\dots}{
        Graphical options for the contour plot.
    }
}

\value{
    This function returns no value.
}

\details{
	Currently, the calculation of contours based on one or more optional arguments (dist, CL, or dof) has not been implemented.
	If optional arguments are provided the function will return an error indicating that the C2P object cannot be found, because it was never created.
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