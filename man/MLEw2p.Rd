\name{MLEw2p}
\alias{MLEw2p}

\title{Maximum Likelihood Estimate for 2-parameter weibull distributions}

\description{
 Determination of Weibull fitting parameters, goodness of fit measures and confidence interval bounds with optional graphical display.
}

\usage{
MLEw2p(x, s=NULL, bounds=FALSE, show=FALSE) 
}

\arguments{
  \item{x}{A vector of failure data, or a dataframe wit time, event and optionally qty columns}
  \item{s}{An optional vector of suspension data.}
  \item{bounds}{A logical argument defining whether confidence interval bounds should be calculated by pivotal analysis.} 
  \item{show}{A logical argument defining whether a simple graphical output is desired.}
}

\details{
 This function is intended to provide a simple casual method of standard weibull analysis based on default methods, without options. 
}

\value{
 When the bounds argument is set to FALSE this function returns a vector with named elements for Eta, Beta, Rsqr,
 AbPval (Abernethy's P-value), and LL (log-likelihood).  When the bounds argument is set to TRUE a list is returned 
 containing the vector as described and a dataframe of confidence interval bound values at a fixed set of descriptive quantiles,
 suitable for comparison with other software.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
fit<-MLEw2p(failures, suspensions)
}

\keyword{ regression reliability }