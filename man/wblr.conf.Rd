\name{wblr.conf}
\alias{wblr.conf}


\title{Add Confidence to \code{wblr} Objects}
\description{
    This function adds confidence calculations to
    various entities in \code{wblr} objects.
}
\usage{wblr.conf(x,which="all",\dots)}
\arguments{
    \item{x}{Object of class \code{"wblr"}.}
    \item{which}{Calculate which fit in the \code{wblr} object will be processed.}
    \item{\dots}{Options for calculating confidence, and for plotting the results.}
}
\details{
    This function adds confidence calculations to various entities in
    \code{wblr} objects and adds them to the object alongside any pre-existing
    confidence calculations.
    
    Additional options for calculating B-life confidence are passed with:

    \describe{
        \item{\code{cl}}{
            Confidence level: A single number from the interval \code{[0,[1}
            specifying the confidence level for various confidence calculations.

            Defaults to \code{0.9}.
        }
        \item{\code{conf.blives.sides}}{
            Either \code{"lower"}, \code{"upper"} or \code{"double"},
            specifying the type of bound(s) to be calculated.

            Defaults to \code{c("double")}, the other options are currently
            not implemented.
        }
        \item{\code{unrel.n}}{
            An integer controlling the amount of unreliability levels for
            which B-life confidence bounds are calculated and ultimately plotted.
            
            Higher numbers will result in smoother confidence bounds. In any
            case, confidence intervals will be calculated for:
            \itemize{
                \item the B-lives at unreliability levels specified with option \code{\link{unrel}}
                \item the B-life at \code{50 [\%]} unreliability
                \item the B-life at the calculcate characteristic life
                    or logmean (depending on the fitted distribution)
            }
            
            Note: When plotting fits and confidence bounds that are adjusted with
            a threshold (see option \code{"threshold"}), it is often the case that
            the bounds appear to be cut of on the left. This can be countered by 
            dramatically increasing \code{unrel.n}, resulting in confidence
            bounds that extend to the edge of the plotting area.

            Defaults to \code{25}.
        }
        \item{\code{conf.what}}{
            A vector of class \code{"character"} describing for which entities
            that confidence should be calculated.

            Defaults to \code{c("blives")}, the only type currently supported.
        }
        \item{\code{unrel}}{
            An unordered numeric vector with unreliability levels for which
            B-life confidence will be calculated.

            Defaults to \code{c(0.1,0.05,0.01)}.
        }

        \item{\code{method.conf.blives}}{
            A vector of class \code{"character"} describing the technique to be
            used for calculating confidence for B-lives. Possible values are
            \code{"bbb"} (Beta Binomial confidence bounds),
            \code{"lrb"} (Likelihood Ratio confidence bounds) and
            \code{"mcpivotals"} or \code{"mcpivotal"} (Monte Carlo Pivotal
            confidence bounds).
            
            Monte Carlo Pivotal confidence bounds use a large number of
            simulations to calculate the confidence bounds. See option
            \code{"S"} for more info.
            
            Defaults to \code{c("mcpivotals")}.
        }
        \item{\code{S}}{
            An integer describing the number of Monte Carlo simulations on
            which the Monte Carlo pivotal confidence bounds and calculation
            of the "prr" goodness-of-fit indicator are based.
        
            High values are needed for good confidence bounds at the lower
            end of the fitted model, especially for data with heavy censoring.

            Note that \code{S >= 100} and that \code{S} must be divisible by 10.

            Defaults to \code{10000}.
        }
        \item{\code{in.legend}}{
            Logical value controlling the inclusion of confidence calculation
            results in the legend. 
            
            If \code{in.legend=FALSE} is passed ,
            the resulting confidence calculations will be omitted from the legend.
            
            Defaults to \code{TRUE}.
        }
    }

    Additionally, one can pass any options available from \code{options.wblr},
    such as \code{col} or \code{is.plot.legend}. The graphical options
    will be used when plotting the (life-)time observations using \code{plot.wblr}.
}
\value{
    The function returns its argument \code{x}, extended with the confidence
    calculations and any optional graphical and calculation arguments
    as passed to the function.
}
\author{Jurgen Symynck \email{jusy@openreliability.org}}
\note{
   \itemize{
        \item Currently, only \code{which = "all"} is supported, meaning that a
            call to \code{wblr.conf} attempts calculation of confidence for all
            fits in the \code{wblr} object.
        \item Currently, only \code{conf.what = "blives"} and
            \code{conf.blives.sides = "double"} are supported.
    }
}


\examples{
## full dataset ##
da1 <- wblr(runif(10,100,1e4),label="Complete data")
da1 <- wblr.fit(da1)
da1 <- wblr.conf(da1,method.conf.blives="mcpivotals",col="red")
da1 <- wblr.conf(da1,method.conf.blives="bbb",col="orange")
da1 <- wblr.conf(da1,method.conf.blives="lrb",col="yellow3")
print(da1$fit[[1]]$conf$blives[[1]])
plot(da1,main="Comparison between MC Pivotal bounds and BB Bounds")

## censored dataset: generates a warning for MC Pivotal confidence bounds ##
da2 <- runif(8,100,1e4)
da2 <- wblr(fail=da2,susp=rep(max(da2),2),label="Type II censored data")
    # generate a 'type 2' censored dataset
da2 <- wblr.fit(da2)
da2 <- wblr.conf(da2,method.conf.blives="mcpivotals",col="blue1")
da2 <- wblr.conf(da2,method.conf.blives="bbb",col="steelblue")
da2 <- wblr.conf(da2,method.conf.blives="lrb",col="cyan3")
plot(da2,main="Comparison between different bound types.")

## show variability in Monte Carlo Pivotal bounds with low S ##
da3 <- wblr(rweibull(5,3,1000))
da3 <- wblr.fit(da3)
for(i in 1:20) da3 <- wblr.conf(da3,S=1000,lwd=1,col="red")
    # just keep adding bounds to the wblr object...
plot(da3,is.plot.legend=FALSE,
    main="Variability in MC Pivotal Conf. Bounds for S=1000")
}
