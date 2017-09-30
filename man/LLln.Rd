\name{LLln}
\alias{LLln}

\title{Log Likelihood for log-normal fitted data}

\description{
 Determination of log-likelihood values for the basic distributions covered by R should be ubiquitous, but an implementation is provided here to enable likelihood ratio testing particularly for comparison of 3-parameter optimized models with their 2-parameter counterpart.
}

\usage{
LLln(x, s=NULL, Mulog, Sigmalog)
}

\arguments{
  \item{x}{A vector of failure data.}
  \item{s}{An optional vector of suspension data.}
  \item{Mulog}{The mean parameter from a log-normal fit.}
  \item{Sigmalog}{The standard deviation parameter from a log-normal fit.}
}

\details{
 Function LLln can only handle datasets with failure and suspension data. For data including intervals \code{wblrLikelihood} must be used.
 This function is somewhat unique among likelihood functions in that it will scan the suspension data for
 removal of negative values. Such negative values are likely to be entered as data is provided
 for a 3-parameter fit by explicitly providing the original data with vector subtraction by the optimized t0.
 The primary intention for likelihood determination is to permit likelihood ratio testing for comparison
 of 3-parameter optimized models with their 2-parameter counterpart.}

\value{
Returns a log-likelihood value.
}

\references{
  Dr. Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition"  
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
fit<-lslr(getPPP(failures, suspensions),dist="lnorm")
LL<-LLln(failures, suspensions, fit[1], fit[2])
}

\keyword{ likelihood reliability }