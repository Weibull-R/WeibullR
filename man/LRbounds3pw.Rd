\name{LRbounds3pw}
\alias{LRbounds3pw}

\title{Likelihood Ratio Bounds for the 3-Parameter Weibull}

\description{
 Generates bounds across a double-sided confidence interval based on a series of
 likelihood ratio contours at varied potential values for the third parameter.}

\usage{
LRbounds3pw(x, s=NULL, CL=0.9, DF=1 ,ptDensity=100, tzpoints=10,
 RadLimit=1e-5, listout=FALSE, show=c(FALSE,FALSE))
}

\arguments{
  \item{x}{
	Either a dataframe containing at least \code{$time} and \code{$event}
	columns, or a vector of class \code{"numeric"} or \code{"integer"} with
	(life-)time observations.
  }
  \item{s}{An optional vector of suspension data.}
  \item{CL}{A scalar for the double-sided confidence interval of interest. Default = 0.9, for 90% double-sided bounds.
  }
  \item{DF}{An integer value indicating degrees of freedom to apply to the Chi square test, which defaults to dof=1 for
  confidence interval bound use. Should be set to 2 for comparison of two models each with 2 parameters.
  }
  \item{ptDensity}{an integer value for the number of points to be plotted around the circumference of the contour.}
  \item{tzpoints}{an integer value for the number of potential values for the third parameter to be modeled.
  }
  \item{RadLimit}{A convergence limit for radius determination as contours points are found.
  }
   \item{listout}{A logical determining whether a list of development objects will be displayed.
  }
   \item{show}{A logigal vector determining whether to display [1] a modified data probability plot [2] a plot of all contours for the potential third parameter cases sampled.
  } 
}

\value{
  A dataframe holding values of percentiles, lower, datum, and upper values of the bound.
  If listout is set to TRUE the following development items are returned in a list:
  [[1]]dataframe holding values of percentiles, lower, datum, and upper values of the bound.
  [[2]] list of 2 dataframes holding the lower and upper bound values for all potential values of the third parmameter,
  [[3]] a list of all the contour points determined for each potential value of the third parameter, and 
  [[4]] x and y ranges for placing all of the generated contours in a single plot.
}

\details{
 This function is temporarily placed for evaluation and testing. The functionality will be incorporated into
 LRbounds in future CRAN submission.
 
 As beta and eta parameters form a log-likelihood mound peaking at the maximum estimate fit. Often a profile is
 generated as a vertical cross section trough the mound. A contour is a horizontal slice through the mound at a
 log likelihood value matching the likelihood ratio defined by max(LL)- qchisq(CL,DF)/2. The beta and eta parameter pairs
 forming this contour can generate a series of fitted lines on the weibull canvas. Likelihood ratio bounds are determined
 by capturing the furthest extent of these fitted lines.
 
 The confidence interval bounds reflecting the uncertainty of the third parameter are determined by generating
 a series of bounds across a range of potential third parameter values. Each of the contours underlying these bounds
 match a ratio using the highest log-likelihood from third parameter optimization. The final bound selection is formed 
 by the furthest extent of the compiled bounds from varied third parameter values.
 
 The example below runs, but exceeds execution time permitted by check --as-cran.
}
\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York
  
  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition"
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
\dontrun{
set.seed<-1234
test50<-rweibull(50, shape=2, scale=100)+25
add25<-rweibull(25, shape=.9, scale=50)+9
test75<-c(test50,add25)
lrb75<-LRbounds3pw(test75,show=TRUE)
}
}

\keyword{ regression reliability }