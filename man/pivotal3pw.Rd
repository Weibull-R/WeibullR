\name{pivotal3pw}
\alias{pivotal3pw}

\title{Pivotal (Parametric Bootstrap) Rank Regression Bounds for the 3-Parameter Weibull}

\description{
 Generates bounds across a double-sided confidence interval based on random
 samples generated on all three parameters of an optimized 3p fit.}

\usage{
pivotal3pw(x, s=NULL, CI=0.9, unrel=NULL, S=1000, listout=FALSE, show=FALSE)
}

\arguments{
  \item{x}{A vector of class \code{"numeric"} or \code{"integer"} with (life-)time observations.
  }
  \item{s}{An optional vector of suspension data.
  }
  \item{CI}{A scalar for the double-sided confidence interval of interest. Default = 0.9, for 90% double-sided bounds.
  }
  \item{unrel}{An optional vector of unreliability values to be used as the descriptive quantiles at which the bounds will be calculated.
  }
  \item{S}{The number of random samples to be drawn for Monte Carlo simulation.
  }
  \item{listout}{A logical determining whether a list of development objects will be displayed.
    *** currently not implemented - only the bounds dataframe is returned.
  }
  \item{show}{A logigal vector determining whether to display plots of the bounds.
  } 
} 

\details{
 This function is temporarily placed for evaluation and testing. The functionality will be incorporated into
 pivotal.rr in future CRAN submission.
 
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
piv75<-pivotal3pw(test75,show=TRUE)
}
}

\keyword{ regression reliability }