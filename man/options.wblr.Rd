\name{options.wblr}
\alias{options.wblr}

\title{
    Options list for \code{wblr} Objects
}
\description{
    This function handles the various calculation, printing and plotting options
    for wblr objects.
}
\usage{options.wblr(\dots)}
\arguments{
    \item{\dots}{
        Options for calculating, printing and plotting \code{wblr} objects using the \pkg{WeibullR}
        package.
    }
}


\value{
   Executing \code{options.wblr} without arguments returns a named list
   containing the currently active global options of the \code{wblr} object.\cr
   When arguments are supplied, these are returned in a named list.
}

\details{
    Typical usage of \code{wblr} objects involves calling a sequence of functions each forming
	the object, then modifying it in this typical manner:

    \code{da <- wblr(c(10,11,27), col="red")}\cr
    \code{da <- wblr.fit(da, col="darkgreen")}\cr
    \code{da <- wblr.conf(da), col="blue"}\cr
    \code{plot(da)}

    The correct time to specify an option is when it is needed for the first
    time. For example, when the color setting option \code{col = "red"} is
    passed as an argument of function \code{wblr}, it will be used for
    data points, fitted lines and confidence interval bounds. If supplied to \code{\link{wblr.conf}},
    only the confidence bounds will have the re-specified color, hereby overriding
    any previously inherited color settings from \code{\link{wblr.fit}}
    or \code{\link{wblr}}.

    Do not call \code{options.wblr} in between these functions because some
    options are locked and cannot be altered further in this chain. This is an
    implication of the way the \code{wblr} object is structured.

    As a function, options.wblr, borrows its internal structure from
    the \code{\link[graphics]{par}} function of package \pkg{graphics}.
    It can be used independently of the typical wblr, wblr.fit, wblr.conf sequence,
    but this is discouraged. However, when used independently the following guidance is
    provided:

        \code{options.wblr()}

        Returns the currently used options and their values.


        \code{options.wblr()$dist }
        \code{options.wblr("dist")}

        Returns the current value of an option.

        \code{options.wblr(ci=0.95) }
        \code{options.wblr(list(ci=0.95,S=5e4))}

        Sets the specified options.

    Currently, there is no way to reset the options to the default values using
    this function.  One might, before changing any options, store the option
    list in a temporary variable like

    \code{wblr.defaults <- options.wblr()}

    for restoring it later by running \code{options.wblr(wblr.defaults)}.

    The function creates a globally accessible list named  \code{options_wblr},
    holding the options.  One should always use the  \code{options.wblr}
    function to access the option list, do not access this list directly.

}

\section{wblr options}{
Options specific to initial \code{\link{wblr}} object creation with data prepared for graphical display
    \describe{
        \item{\code{dist}}{
            A character string defining the distribution target. When used to establish the basis for
			contour mapping (without using wblr.conf with method.conf="lrb") only "weibull" (default)
			and "lognormal" are recognized.\cr
			Also used with \code{\link{wblr.fit}} for specific fitting control.
        }
        \item{\code{pp}}{
            Plotting position method, it is a character string describing the method of determining
            vertical plot positions. Implemented methods are "median" (default), "benard","hazen","mean", "kaplan-meier", and "blom".
        }
        \item{\code{rank.adj}}{
		    The method employed for determining rank of failures when suspensions (right censored data)
			are present in the data set. Implemented methods are "johnson" (default) and "KMestimator".
		}
        \item{\code{ties.handler}}{
		    The method employed for handling duplicate values in the data set.\cr
			Implemented methods are "none" (default) "highest", "lowest", "mean", and "sequential".\cr
			It is expected that ties handling will be applied to large data sets that will be fitted
			using the maximum likelihood estimation method, where the effect is only on the graphical
			presentation. Employing a ties handler on a rank regression model will effectively remove data
			from the data set, which is likely not intended.\cr
			Use of simply \code{ties} as an argument to function \code{wblr} will silently be accepted as \code{ties.handler}.
		}
		\item{Options for graphical control over data points see \code{\link[graphics]{par}}}{
			\itemize{
				\item{\code{pch}}{ Point choice defaults to \code{1}. For more info, see
					\code{\link[graphics]{points}}.
				}
				\item{\code{cex.points}}{ Point size defaults to \code{1}.
				}
				\item{\code{lwd.points}}{ Line width defaults to \code{2}.
				}
			}
		}
		\item{Independent graphical control over interval lines}{
			\itemize{
				\item{\code{interval.col}}{ Color defaults to \code{"black"}.
				}
				\item{\code{interval.lty}}{ Line type, defaults to \code{"dashed"}.
				}
				\item{\code{interval.lwd}}{ Line width defaults to \code{1}.
				}
			}
		}
	}
Options specific to \code{\link{wblr.fit}}
    \describe{
        \item{\code{dist}}{
            A character string defining the distribution target. Recognized values are
			"weibull" (default), "lognormal","lnorm","lognormal2p", "weibull2p","lognormal3p", and "weibull3p".
			}
        \item{\code{method.fit}}{
            A vector of class \code{"character"} with fitting options. Recognized values are "rr-xony" (default),
			"rr","rr-yonx", "mle","mle-rba", and "mle-unbias".
        }
	}

Options specific to \code{\link{wblr.conf}}
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
General graphical options
	\describe{
		\item{Graphical control for fitted lines, confidence bounds, or contours can be set at
		 wblr, wblr.fit, or wblr.conf functions for individual object control. See \code{\link[graphics]{par}}}{
			\itemize{
				\item{\code{col}}{ Color defaults to \code{"black"}.
				}
				\item{\code{lty}}{ Line type, defaults to \code{"solid"}.
				}
				\item{\code{lwd}}{ Line width defaults to \code{2}.
				}
			}
		}
		\item{Graphical controls typically used with \code{plot} or \code{\link{plot.wblr}}}{
			\itemize{
				\item{\code{canvas}}{
					 The plotting canvas to be used. This does not necessarily have to match the fit distribution.
					Only \code{"weibull"} (default) or \code{"lognormal"} are recognized.
				}
				\item{\code{mar}}{ Margins defaults to \code{c(5.1,4.1,5.1,2.1),}.
				}
				\item{\code{main}}{ Title, defaults to \code{"Probability Plot"}.
				}
				\item{\code{main.contour}}{ Contour plot title defaults to \code{"Contour Plot"}.
				}
				\item{\code{sub}}{ Subtitle defaults to \code{NULL}.
				}
				\item{\code{sub.contour}}{ Contour subtitle defaults to \code{NULL}.
				}
				\item{\code{xlim}}{ Plot x limits override to be presented as a vector c(lo,hi), default \code{NULL}.
				}
				\item{\code{ylim}}{ Plot y limits override to be presented as a vector c(lo,hi), default \code{NULL}.
				}
				\item{\code{xlab}}{ X axis label defaults to \code{"Time To Failure"}.
				}
				\item{\code{ylab}}{ Y axis label defaults to \code{"Unreliability [\%]"}.
				}
				\item{\code{coordinate.text.size}}{ default \code{0.7}.
				}
				\item{\code{signif}}{ Used to control display of numbers in Legend, default \code{4}.
				}
				\item{\code{col.grid}}{ Color for chart gridlines defaults \code{"grey"}.
				}
				\item{\code{is.plot.grid}}{ default \code{TRUE}.
				}
				\item{\code{is.plot.fit}}{ default \code{TRUE}.
				}
				\item{\code{is.plot.pp}}{ default \code{TRUE}.
				}
				\item{\code{is.plot.ppcoordinates}}{ default \code{FALSE}.
				}
				\item{\code{is.plot.legend}}{ default \code{TRUE}.
				}
				\item{\code{legend.text.size}}{ default \code{0.7}.
				}
				\item{\code{label}}{ defaults to \code{""}.
				}
				\item{\code{in.legend}}{ default \code{TRUE}.
				}
				\item{\code{in.legend.blives}}{ default \code{TRUE}.
				}
				\item{\code{in.legend.gof}}{ default \code{TRUE}.
				}
				\item{\code{is.plot.cb}}{ default \code{TRUE}.
				}
				\item{\code{persistent}}{ default \code{TRUE}.
				}
			}
		}
	}
}

\examples{
## backup options ##
wblr.defaults <- options.wblr()

## setting new options ##
options.wblr(S=5e5,ci=0.99)
%%options.wblr(sub="Testing options.wblr()")

## listing options ##
options.wblr()
options.wblr()$main

## restore options ##
options.wblr(wblr.defaults)
}
