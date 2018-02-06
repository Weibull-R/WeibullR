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

\usage{plot_contour(x,\dots)}

\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of same objects.
    }
    \item{\dots}{
        Graphical options for the contour plot.
    }
}

\value{
    This function returns no value.
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