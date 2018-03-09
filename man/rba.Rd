\name{rba}

\alias{rba}

\title{Reduced Bias Adjustment for Weibull and Lognormal MLE} 

\description{
\code{rba} This function generates reduction factors based on the median bias of the Weibull Beta parameter MLE (C4^3.5),
 alternatively based on the mean bias(C4^6), and a single factor for both median and mean of the symmetrical lognormal distribution.  
}

\usage{
 rba(Qx, dist="weibull", basis="median")
}

\arguments{
  \item{Qx}{The quantity of failures in data (right-censored data, or suspensions excluded)}
  \item{dist}{A string defining a distribution to be fit. Implemented distributions are "weibull" (default) and "lognormal".
  }
  \item{basis}{a string value indicating the basis for bias reduction, "median" is the default, but "mean" may be chosen.
}
}

\value{
   A factor to be multiplied to the MLE Beta or log-Sigma for bias reduction.
}

\details{
  Many references discuss the bias reduction in terms of mean reduction.  Abernethy has promoted consideration of the
  RBA in terms of median bias reduction.
}

\references{
  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition"
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
MLEfit<-mlefit(mleframe(failures,suspensions))
MLE_RBA<-c(MLEfit[1],MLEfit[2]*rba(length(failures)))
}



\keyword{ likelihood }