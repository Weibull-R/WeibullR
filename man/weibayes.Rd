\name{weibayes}
\alias{weibayes}

\title{
    Fitting for Minimal Failure Datasets
}

\description{
    A simplistic fitting method also called "1-parameter Weibull" based on a "prior" known beta value. 
}

\usage{
weibayes(x, s=NULL, beta)
}

\arguments{
  \item{x}{
	Either a vector of class \code{"numeric"} or \code{"integer"} with
	(life-)time observations, or a dataframe containing at least \code{$time} and \code{$event}
	columns and optionally a \code{$qty} column.
  }
  \item{s}{An optional vector of suspension data.}
    \item{beta}{
        A weibull shape parameter based on prior like kind study.
    }
}

\value{
    Returns the Weibull scale, Eta, paramater value.
}

\details{
   This function implements a very simplistic relationship defined by two lines of R code:\cr
   	 \code{t_eta<-(times^beta)/nfail}\cr
     \code{out_val<-sum(t_eta)^(1/beta)}\cr
   Note: for a single failure dataset, with beta = 1, the return is simply the sum of the failure time value plus sum of all suspension times.
   Should no failures exist, then the nfail will be set to 1 for a pessimistic solution.\cr
   }

\examples{
fail<-5
susp<-rweibull(10, 1, 10)
eta<-weibayes(fail, susp, beta=1)
}