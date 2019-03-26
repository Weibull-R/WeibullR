\name{BBB}
\alias{BBB}

\title{Beta Binomial Bounds}

\description{
 Generates confidence inference bounds by means of the beta binomial distribution applied to data ranks.}

\usage{
BBB(xdata, xfit,n=NULL, CI=.90, unrel=NULL, type="horizontal", nknots=NULL)
}

\arguments{
  \item{xdata}{A list containing either dpoints and/or dlines such as available from a wblr object.
  }
  \item{xfit}{A list containing dist and fit elements. A suitable object can be extracted from a wblr object containing the resuult of at least one modifying wblr.fit call.
	Function \code{xfit} has been provided to return such a list.
  }
  \item{CI}{A scalar for the double-sided confidence interval of interest. Default = 0.9, for 90% double-sided bounds.
  }
  \item{unrel}{An optional vector of unreliability values to be used as the descriptive percentiles at which the bound quantiles will be calculated.
  }
  \item{type}{A string indicating the adjustment direction(s) for application of the beta binomial. Since both "horizontal" and "vertical" directions
   on the rank positions of failure data are considered valid, a smoothed spline can be created by appling both directions using "valid".
   An extrapolation of the smoothed spline is provided by using "extrapolated".
  }
  \item{nknots}{An optional number of knots to be used to form the smooth spline in the event that disjoints are formed
   at the junction of horizontal and vertical bounds. Usually, applying 75 to 80% of the number of failure points (including intervals) will
   result in a pleasing smoothness.
  } 
}

\value{
  Returns a dataframe holding values of Prob, Lower, Fit, and Upper values for the bounds.
}

\details{
 This non-parametric approach calculates confidence bounds for quantiles and/or
 failure probabilities applying the CFD of the beta distribution at ranks and reverse ranks of the failure data including intervals.
 Although a result is provided for 3-parameter models, it does not reflect uncertainty in the third, translation parameter.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
  
  Tim-Gunnar Hensel, (2017) "weibulltools", A package on CRAN.
}

\examples{
set.seed(4321)
obj<-wblr.fit(wblr(rweibull(10,1,1)))
beta_binomial_bounds<-BBB(obj$data,xfit(obj))
}

\keyword{ regression reliability }