\name{MRRln3p}
\alias{MRRln3p}

\title{Median rank regression for log-normal distribution with third parameter optimization}

\description{
 Determination of log-normal fitting parameterswith third parameter optimization.
 Goodness of fit measures and confidence interval bounds can be returned with optional graphical display.
}

\usage{
MRRln3p(x, s=NULL, bounds=FALSE, show=FALSE) 
}

\arguments{
  \item{x}{A vector of failure data.}
  \item{s}{An optional vector of suspension data.}
  \item{bounds}{A logical argument defining whether confidence interval bounds should be calculated by pivotal analysis.}
  \item{show}{A logical argument defining whether a simple graphical output is desired.}
}

\details{
 This function is intended to provide a simple casual method of standard lognoirmal fitting based on default methods,
 without options.  It also provides an example for handling the pivotal values returned from pivotalMC and the likelihood ratio test. 
}

\value{
 When bounds is set to FALSE this function returns a vector with named elements for Eta, Beta, Rsqr,
 AbPval (Abernethy's P-value), and LL (log-likelihood).  When the bounds argument is set to TRUE a list is returned with
 the vector as described and a dataframe of confidence interval bound values at a fixed set of descriptive quantiles,
 dq<-c(.01, .02, .05, .10, .15, .20, .30, .40, .50,  .60, .70, .80, .90, .95, .99), suitable for comparison with other software.
}

\references{
  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition"
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
fit<-MRRln3p(failures, suspensions)
}

\keyword{ regression reliability }