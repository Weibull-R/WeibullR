\name{wblr}
\alias{wblr}

\title{
    Create a \code{wblr} Object for Life Data Analysis
}
\description{
    This function creates an object of class \code{"wblr"} for further processing
    by the other functions of \pkg{wblr}.
}
\usage{wblr(x, s=NULL, interval=NULL,\dots)}

\arguments{
    \item{x}{
        Either a dataframe containing at least \code{$time} and \code{$event}
        columns and optionally a \code{$qty} column, or a vector of class \code{"numeric"} or \code{"integer"} with
        (life-)time observations.
	}
	\item{s}{
		An optional vector of right-censored data, or suspensions. 
	}
	\item{interval}{
		An optional dataframe of interval data having columns specifically named "left" and "right". 
		Left values are the last time at which no failure was evident and may be zero for discovery.
		Right values are the earliest time at which failure was observed.
	}
    \item{\dots}{
        Graphical options for plotting the \code{wblr} object.
        See section "Details" for other data passing arguments.
    }
}

\details{
    There are several methods to passing arguments for building an \code{wblr}
    object.
    \itemize{
        \item If argument \code{x} is of class \code{"data.frame"},
            then it must contain \code{$time} and \code{$event}
            columns. Additional columns in the dataframe will be ignored.
        \item When a single unnamed vector of class \code{"numeric"}
            or \code{"integer"} is supplied, it is treated as a vector
            of (life-)time observations.
        \item If argument \code{time} or \code{fail} is provided, it is treated as
            a vector of (life-)time observations. Take care NOT to supply both
            \code{time} and \code{fail} in the same function call.
        \item If argument \code{event} is provided, it is treated as
            a vector of event indicators with possible values of
            \code{0} and \code{1}. See section "Value" for more details on
            event vectors.
        \item If the \code{x} argument is not provided as a dataframe and \code{susp} is provided, 
			it is treated as a vector of right-censored (life-)time observations 
            (also called suspended observations or suspensions).
    }
\code{wblr} always generates (probability) plot positions for graphically
    displaying the (life-)time observations and for (possible) later usage
    by \code{\link{wblr.fit}}. The following optiona arguments are most appropriate for
	passing in with \code{wblr}:
    \describe{
        \item{\code{dist}}{
            A character string defining the distribution target. When used to establish the basis for
			contour mapping (without using \code{\link{wblr.conf}} with method.conf="lrb") only "weibull" (default)
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

    Subsequent calls to \code{wblr.fit} and \code{wblr.conf} will inherit these options.
}

\value{
    A named list of class \code{"wblr"}. The first list
    item (\code{$data}) is a list with up to least three items:
    \describe{
        \item{\code{$lrq_frame}}{
            A dataframe containing the provided data formatted with \code{"left"}, \code{"right"}, and \code{"qty"} columns.
			This is the output of WeibullR function \code{"mleframe"}.
        }
        \item{\code{$data$dpoints}}{
            A dataframe contianing graphical data for exact failure point with their probability plotting positions and adjusted ranks.
        }
        \item{\code{$data$dlines}}{
            If interval data has been provided this dataframe will contain the graphical data for display 
			similar to {\code{$data$dpoints}}, but with endpoints t1 and t2 for the interval.
        }
    }
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
## These code lines all generate the same object ##
wblr(c(500,1200,900,1300,510))
wblr(time=c(500,1200,900,1300,510))
## this input format works, but not recommended.
wblr(time=c(500,1200,900,1300,510),event=c(1,1,1,1,1))
wblr(fail=c(500,1200,900,1300,510))
wblr(fail=c(500,1200,900,1300,510),susp=c())
da1 <- data.frame(
    serial=c("S12","S16","S17","S3","S5"),
    time=c(500,1200,900,1300,510),
    event=c(1,1,1,1,1))
## it is best practice set named objects
obj1 <- wblr(da1,label="complete dataset",pch=3,col="orange3")
obj2 <- wblr(da1,label="complete dataset",pch=4,pp="benard",col="red")


## Generate a similar dataset, but with suspensions ##
wblr(time=c(500,1200,900,1300,510),event=c(1,1,1,0,0))
wblr(data.frame(time=c(500,1200,900,1300,510),event=c(1,1,1,0,0)))
wblr(fail=c(500,1200,900),susp=c(1300,510))
wblr(time=c(500,1200,900),susp=c(1300,510))
da3 <- wblr(fail=c(500,1200,900,1300,510),
    event=c(1,1,1,0,0),label="censored dataset",pch=1,col="blue")

## plot datasets ##
\dontrun{
plot.wblr(list(da1,da3))
}

}
