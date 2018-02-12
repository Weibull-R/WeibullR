\name{plot_contour}
\alias{plot_contour}

\title{
    Plotting of Likelihood Ratio Contours from \code{wblr} Objects
}

\description{
    Plots likelihood ratio contours from a \code{wblr} object or a list of \code{wblr} objects. 
}

\usage{plot_contour(x, CL=NULL), AL=TRUE, main="", sub="")}

\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of \code{"wblr"} objects.
    }
    \item{CL}{An optional argument for establishing the confidence limit(s) for calculated contours. A vector of values may be provided to generate concentric contours on a single object.
    }
    \item{AL}{A Logical defining whether parameters of lognormal contours should be presented as AntiLog values.
    }
    \item{main}{A string override for the the plot title that would otherwise be drawn from the contour.main option in the first \code{"wblr"} object provided.
	}
    \item{main}{A string override for the the plot subtitle that would otherwise be drawn from the contour.sub option in the first \code{"wblr"} object provided.
	}
}

\value{
    This function returns the list of contours that were prepared for plotting.
}

\details{
   If no CL argument is provided, the contours generated in the provided \code{wblr} objects for likelihood confidence bounds will be plotted.
   If a vector is provided for CL, multiple contours for each provided object will be calculated to display concentric contours according to values in the CL vector.
   In all cases options for distribution, degrees of freedom and graphic options will be extracted from the underlying objects.
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