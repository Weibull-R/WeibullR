\name{wblr.conf}
\alias{wblr.conf}


\title{Add Confidence Interval Bounds to \code{wblr} Objects}
\description{
    This function adds confidence calculations to the last fit entity in a \code{wblr} object.
}
\usage{wblr.conf(x,\dots)}
\arguments{
    \item{x}{Object of class \code{"wblr"}.}
    \item{\dots}{Options for calculating confidence, and for plotting the results.}
}
\details{
    This function adds confidence calculations to the last fit entity in a
    \code{wblr} object and adds the result to the object alongside any pre-existing
    confidence calculations.\cr

	Certain confidence calculations pertain to specific fit methods. The pivotal.rr confidence method
	can only apply to a rank regression fit, while a likelihood ratio calculation can only be
	performed on an mle fit.\cr

    Option arguments suitable for passing into wblr.conf:
	\describe{
        \item{\code{method.conf}}{
            A character string describing the techniques used for calculating confidence interval bounds.
            Implemented methods are "pivotal-rr" (default), "bbb", "fm" "fmbounds", and "lrb".
			Methods must conform to the method.fit in the wblr.fit call immediately preceding the wblr.conf call.
			Method "pivotal-rr" requires a rank regression fit method. Methods "fm", "fmbounds, and "lrb" require a mle based fit.
        }
        \item{\code{dq}}{
			A named series of quantiles at which confidence interval bounds will be calculated.
			\itemize{
				\item{\code{"abrem"}} {
					Default. This is the original default by Jurgen Symynck for predecessor package abrem
					it produces evenly spaced points across the y limits of a weibull canvas
					attempting to hold a constant number of points (see \code{num_dq} below).
				}
				\item{\code{"minitab"}} {
					Quantiles matching Minitab(TM) unchangeable defaults (27 values).
				}
				\item{\code{"supersmith"}} {
					Quantiles for comparison with SuperSMITH(TM) (limit of 15 values)
				}
				\item{\code{"user"}} {
					Provides for a user defined series of quantiles. (see \code{user_dq} below).
				}
			}
		}
        \item{\code{num_dq}}{
            The number of points used for the \code{"abrem"} dq determination.
        }
        \item{\code{user_dq}}{
            A vector of quantiles set by user. Default \code{c(seq(.01,.09,by=.01),seq(.10,.90,by=.10),seq(.91,.99, by=.01))}.
        }
        \item{\code{ci}}{
            The double-sided confidence interval, also chi sq confidence level for likelihood ratio.
			Must be in a range <1 && >0, default is \code{0.9}.
        }
        \item{\code{blife.pts}}{
            The probability points at which to report Blife on legend.
        }
		\item{Specific controls for pivotal analysis only}{
			\itemize{
				\item{\code{seed}}{
					The RNG seed integer such that results are duplicated between runs, default is \code{1234}.
				}
				\item{\code{S}}{
					The number of samples to be run during pivotal analysis, default is \code{1e4}.
				}
			}
		}
		\item{Specific controls for likelihood ratio contour and bounds only}{
			\itemize{
				\item{\code{dof}}{
					 The degrees of freedom, \code{dof=1} (default) for confidence interval, \code{dof=2} for comparison.
				}
				\item{\code{ptDensity}}{
					 The number of points calculated to form the contour outline, default is \code{120}.
				}
			}
		}
	}
Additionally, graphical control options for the confidence interval bounds can be passed in with the following options:
	\itemize{ 
		\item{\code{col}}{ Color defaults to \code{"black"}.
		}
		\item{\code{lty}}{ Line type, defaults to \code{"solid"}.
		}
		\item{\code{lwd}}{ Line width defaults to \code{2}.
		}
	}
    Subsequent calls to \code{wblr.fit} and \code{wblr.conf} will inherit these options.
}
\value{
    The function returns its argument \code{x}, extended with the confidence
    calculations and any optional graphical and calculation arguments
    as passed to the function.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
## full dataset ##
da1 <- wblr(runif(10,100,1e4),label="Complete data")
da1 <- wblr.fit(da1)
da1 <- wblr.conf(da1,method.conf="pivotal-rr",col="red")
da1 <- wblr.conf(da1,method.conf="bbb",col="orange")
da1 <- wblr.fit(da1, method.fit="mle")
da1 <- wblr.conf(da1,method.conf="lrb",col="yellow3")

## plot datasets ##
\dontrun{
plot(da1,main="Comparison between different bound types.")
}

## with right-censored data
da2 <- runif(8,100,1e4)
da2 <- wblr(fail=da2,susp=rep(max(da2),2),label="Type II right-censored data")
da2 <- wblr.fit(da2)
da2 <- wblr.conf(da2,method.conf="pivotal-rr",col="blue1")
da2 <- wblr.conf(da2,method.conf="bbb",col="steelblue")
da2 <- wblr.fit(da2, method.fit="mle")
da2 <- wblr.conf(da2,method.conf="lrb",col="cyan3")
\dontrun{
plot(da2,main="Comparison between different bound types with right-censored data.")
}
}


\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", 
  Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}