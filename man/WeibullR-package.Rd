\name{WeibullR-package}
\alias{WeibullR-package}

\docType{package}
\title{
    Weibull-R : Weibull Analysis on R
}
\description{
    Life data analysis in the graphical tradition of Waloddi Weibull
}

\details{
The WeibullR package provides a flexible data entry capability with three levels of usage.
	\describe{
		\item{\emph{Quick Fit} Functions}{Functions with intuitive names \code{\link{MLEw2p}} through \code{\link{MRRln3p}} for preparing simple fits, bounds, and displays using default options. Only data sets with exact failure times and/or suspensions are processed.
		The quick fit functions return a simple named vector of the fitted parameters with appropriate goodness of fit measure(s).
		 Optional preparation of appropriate interval bounds (at 90\% confidence), or a display of fit and bounds are controlled by two final arguments taking logical entry,
		 Such that  a function call like \code{MLEw2p(input_data, T, T)} will generate a plot with the fitted data and confidence interval bounds.
		 When the first logical for bounds is set to TRUE, the returned object will be a list with the fitted parameter vector first and dataframe of bounds values second.
		}
		\item{wblr Object Model}{Construction of a wblr object is initiated by providing a data set through function \code{\link{wblr}}.
		 Modification of the object with the progressive addition of fits and confidence interval bounds is made via functions \code{\link{wblr.fit}} and \code{\link{wblr.conf}}. 
		 Fine control over many aspects of fit, confidence, and display are made possible using a flexible options mechanism.
		 Display for single object models is via S3 methods \code{plot} or \code{contour}, while multiple objects \emph{(provided as a list)} can be displayed on a single plot using \code{\link{plot.wblr}}, \code{\link{plot_contour}}, or \code{\link{contour.wblr}}.
		}
		\item{Backend Functions}{Access to backend functions providing all the functionality of the upper levels of usage are provided as exported functions.
		These functions may provide advanced users with resources to expand analysis further than has been implemented in the WeibullR package.
		}
	}
Data entry is made through the Quick Fit functions, \code{\link{wblr}}, or on the backend through \code{\link{getPPP}} for rank regression, \code{\link{mleframe}} for mle processing.
 In all cases the primary argument \code{x} can be a vector of exact time failures or a dataframe with \code{time}, and \code{event} columns as a minimum. An additional column \code{qty} may optionally be used to record duplicated data.
 If the dataframe entry is not used (in favor of an exact time failure vector), a second argument, \code{s}, can be used to enter a vector of last observed success times for right censored data (suspensions).\cr
 Beyond the entry of the first two data types, interval data (including discoveries with last known success time=0) are entered via argument \code{interval} as a dataframe with columns \code{left}, and \code{right} as a minimum.
 As with the primary argument dataframe entry, an additional column \code{qty} may optionally be used to record duplicated interval data.
 Such interval data entry is not supported with the Quick Fit functions.
}

\author{
David J. Silkworth <djsilk@openreliability.org>
Jurgen Symynck <email withheld>
}