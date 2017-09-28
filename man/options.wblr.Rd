\name{options.wblr}
\alias{options.wblr}

\title{
    Options for wblr Objects
}
\description{
    This function handles the various calculation, printing and plotting options
    for wblr objects.
}
\usage{options.wblr(\dots)}
\arguments{
    \item{\dots}{
        Options for calculating, printing and plotting with the \pkg{wblr}
        package.
    }
}

\details{
    This function borrows its internal structure from
    the \code{\link[graphics]{par}} function of package \pkg{graphics}.
    It can be used in different ways:

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

    \describe{
        \item{\code{cex.points}}{
            A number describing the relative size of the datapoint glyphs.

            Defaults to \code{1}.
        }
        \item{\code{ci}}{
            Confidence interval: A single number from the interval \code{[0,[1}
            specifying the double-sided confidence for various confidence calculations.

            Defaults to \code{0.9}.
        }
        \item{\code{col}}{
            An integer or character string describing the color of a
            graphical entity.
            See function \code{\link[grDevices]{colors}} for available colors.

            Defaults to \code{"black"}.
        }
        \item{\code{col.grid}}{
            An integer or character string describing the color of the grid

            Defaults to \code{"gray"}.
        }

        \item{\code{coordinate.text.size}}{
            A number determining the relative coordinate text label size.

            Defaults to \code{0.7}.
        }
        \item{\code{dist}}{
            A character string with the distribution target for fitting.
            See \code{\link{wblr.fit}} for in-depth discussion of the settings.

            Defaults to \code{"weibull"}.
        }
        \item{\code{in.legend}}{
            Logical value controlling the inclusion of various elements in
            the legend. See \code{\link{wblr.fit}} and \code{\link{wblr.conf}}
            for more details.

            Defaults to \code{TRUE}.
        }
        \item{\code{in.legend.blives}}{
            Logical value controlling the inclusion of B-life levels
            in the legend.

            Defaults to \code{TRUE}.
        }
        \item{\code{in.legend.gof}}{
            Logical value controlling the inclusion of goodness-of-fit
            indicators in the legend. Note that this does \emph{not} prohibit
            calculation of sometimes time consuming \code{prr} values, it only
            prevents inclusion in the legend.

            Defaults to \code{TRUE}.
        }
        \item{\code{is.plot.cb}}{
            Logical value controlling the plotting of various types of
            confidence bounds (if any present in the \code{wblr} object).

            Defaults to \code{TRUE}.
        }
        \item{\code{is.plot.fit}}{
            Logical value controlling the plotting of the fitted line.

            Defaults to \code{TRUE}.
        }
        \item{\code{is.plot.grid}}{
            Logical value controlling the plotting of the grid.

            Defaults to \code{TRUE}.
        }
        \item{\code{is.plot.legend}}{
            Logical value controlling the plotting of the legend.

            Defaults to \code{TRUE}.
        }
        \item{\code{is.plot.pp}}{
            Logical value controlling the plotting of the
            (life-)time observations' probability plot positions.

            Defaults to \code{TRUE}.
        }

        \item{\code{label}}{
            A character string with the title of the legend box, can be used to
            label a dataset or fit.

            Defaults to \code{NULL}.
        }
        \item{\code{legend.text.size}}{
            A number determining the relative legend text size.

            Defaults to \code{0.7}.
        }
        \item{\code{canvas}}{
            A character string describing the type of plotting canvas to use.
            Only \code{"weibull"} (default), and \code{"lognormal"} are implemented.
        }
        \item{\code{lty}}{
            An integer describing the line type of both the fitted
            line and the confidence bounds, if any.
            See function \code{\link[graphics]{par}} for available line types.

            Defaults to \code{1}.
        }
        \item{\code{lwd}}{
            An integer describing the line width of both the fitted
            line and the confidence bounds, if any.

            Defaults to \code{2}.
        }
        \item{\code{lwd.points}}{
            An integer describing the thickness of the datapoint glyphs.

            Defaults to \code{2}.
        }
        \item{\code{main}, \code{main.contour}}{
            A character string; the main title of the plot or contour plot.

            Defaults to \code{"Probability Plot"} or \code{"Contour Plot"} .
        }
        \item{\code{mar}}{
            A numerical vector of the form \code{c(bottom, left, top, right)}
            which gives the number of lines of margin to be specified on
            the four sides of the plot.
            See function \code{\link[graphics]{par}} for more info.

            Defaults to \code{c(5.1,4.1,5.1,2.1)}.
        }
        \item{\code{method.conf}}{
            A vector of class \code{"character"} describing the techniques used
            for calculating confidence interval bounds.
            See \code{\link{wblr.conf}} for in-depth discussion of the settings.

            Defaults to \code{"mcpivotals"} for rank regression fitted models, 
			but will be silently changed to \code{"fm"} for mle fitted models.
        }
        \item{\code{method.fit}}{
            A vector of class \code{"character"} with fitting options.
            See \code{\link{wblr.fit}} for in-depth discussion of the settings.

            Defaults to \code{c("rr","xony")}, but will be silently changed
			to \code{"mle"} when input data includes intervals.
        }
        \item{\code{pch}}{
            An integer or single character describing the plotting symbol, used for
            plotting the datapoints.  For more info, see
            \code{\link[graphics]{points}}.

            Defaults to \code{1}.
        }

        \item{\code{pp}}{
            \code{pp} is short for "probability plot positions", it is a vector of
            class \code{"character"} describing the method of determining plot positions.
            Possible values are \code{"median"} (exact median plot postions
            calculated using the Incomplete Beta function) and  \code{"benard"}
            (a good approximation to the median plot positions). Both options
            can be supplied at the same time but currently, only the first
            element in the \code{pp} vector will be used for
            further fitting by \code{\link{wblr.fit}}

            Defaults to \code{c("median")}.
        }
        \item{\code{S}}{
            An integer describing the number of Monte Carlo simulations on
            which the Monte Carlo pivotal confidence bounds and calculation
            of the "prr" goodness-of-fit indicator are based (See
            \code{\link[pivotals]{prrVal}} of package \pkg{pivotals} for more
            details).

            Defaults to \code{10000}.
        }
        \item{\code{signif}}{
            An integer describing the significant digits of various numbers
            that are displayed in the legend.

            Defaults to \code{4}.
        }
        \item{\code{sub}, \code{sub.contour}}{
            A character string; the subtitle of the plot or contour plot.

            Defaults to \code{""}.
        }
        \item{\code{blife.pts}}{
            An unordered numeric vector with unreliability levels for which
            B-life values will be presented.

            Defaults to \code{c(0.1,0.05,0.01)}.
        }

        \item{\code{xlab}, \code{ylab}}{
            Character strings with the labels for the X and Y axis.

            Default to \code{"Time To Failure"} and
            \code{"Unreliability [\%]"} respectively.
            For contour plots, Default to \code{"Eta"} and \code{"Beta"} respectively.
        }
        \item{\code{xlim}}{
            A numeric vector with two values determining the plotting range
            of the horizontal axis of the plot.
            \code{plot.wblr} calculates horizontal limits automatically
            from its data argument.

            Defaults to \code{NULL}.
        }
        \item{\code{ylim}}{
            A numeric vector with two values determining the plotting range
            of the vertical axis of the plot.  Allowed values come
            from the interval \code{]0,1[}. \code{plot.wblr} calculates
            vertical limits automatically from its data argument; if vertical
            plotting positions do not exceed the standard range of
            \code{c(0.01,0.99)}, then the latter is used.
            Otherwise, the standard range is expanded to include the
            extreme plot positions.

            The standard vertical range is currently hardcoded to \code{c(0.01,0.99)}.

            Defaults to \code{NULL}.
        }
    }
}

\value{
   Executing \code{options.wblr} without arguments returns a named list
   containing the currently active global options of the \code{wblr} object.

   When arguments were supplied, these are returned in a named list.
}
\author{Jurgen Symynck \email{jusy@openreliability.org}}
\note{
    Typical usage of \pkg{wblr} involves calling a sequence of functions like:

    \code{da <- wblr(c(10,11,27))}\cr
    \code{da <- wblr.fit(da)}\cr
    \code{da <- wblr.conf(da)}\cr
    \code{plot(da)}

    Do not call \code{options.wblr} in between these functions because some
    options are locked and cannot be altered further in this chain. This is an
    implication of the way the \code{wblr} object is structured.
    
    The correct time to specify an option is when it is needed for the first
    time. For example, when the color setting option \code{col = "red"} is
    passed as an argument of function \code{wblr}, it will be used for
    datapoints, fits and confidence bounds. If supplied to \code{\link{wblr.conf}},
    only the confidence bounds will have the specified color, hereby overriding
    any previously inherited color settings from \code{\link{wblr.fit}}
    or \code{\link{wblr}}.
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
