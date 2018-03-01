\name{MLEln3p}
\alias{MLEln3p}

\title{Maximum Likelihood Estimate for 3-parameter lognormal distributions}

\description{
 Determination of lognormal fitting parameters with third, translation parameter optimization.}

\usage{
MLEln3p(x, s=NULL, bounds=FALSE, show=FALSE) 
}

\arguments{
  \item{x}{A vector of failure data, or a dataframe wit time, event and optionally qty columns}
  \item{s}{An optional vector of suspension data.}
  \item{bounds}{A logical argument defining whether confidence interval bounds should be calculated by pivotal analysis.} 
  \item{show}{A logical argument defining whether a simple graphical output is desired.}
}

\details{
 This function is intended to provide a simple casual method of standard lognormal fitting based on default methods, without options. 
}

\value{
 When the bounds argument is set to FALSE this function returns a vector with named elements for Mulog, Sigmalog,
 and LL (log-likelihood).  Bounds are not determined for 3p models, hence any bounds argument is simply ignored.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
fit<-MLEln2p(failures, suspensions)
}

\keyword{ regression reliability }