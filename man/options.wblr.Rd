\name{options.wblr}
\alias{options.wblr}
\alias{cex.points}
\alias{cl}
\alias{col}
\alias{col.grid}
\alias{conf.blives.sides}
\alias{conf.what}
\alias{coordinate.text.size}
\alias{dist}
\alias{in.legend}
\alias{in.legend.blives}
\alias{in.legend.gof}
\alias{is.plot.cb}
\alias{is.plot.fit}
\alias{is.plot.grid}
\alias{is.plot.legend}
\alias{is.plot.pp}
\alias{is.plot.ppcoordinates}
\alias{label}
\alias{legend.text.size}
\alias{log}
\alias{lty}
\alias{lwd}
\alias{lwd.points}
\alias{main}
\alias{main.contour}
\alias{mar}
\alias{method.conf.blives}
\alias{method.fit}
\alias{pch}
\alias{persistent}
\alias{pivotals}
\alias{pp}
\alias{S}
\alias{signif}
\alias{sub}
\alias{threshold}
\alias{unrel}
\alias{unrel.n}
\alias{verbosity}
\alias{xlab}
\alias{xlim}
\alias{ylab}
\alias{ylim}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    \pkg{wblr} Options
}
\description{
    This function handles the various calculation, printing and plotting options
    of the \pkg{wblr} package.
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

        \code{options.wblr(cl=0.95) }
        \code{options.wblr(list(cl=0.95,S=5e4))}

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
        \item{\code{cl}}{
            Confidence level: A single number from the interval \code{[0,[1}
            specifying the confidence level for various confidence calculations.

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
        \item{\code{conf.blives.sides}}{
            Either \code{"lower"},  \code{"upper"} or \code{"double"},
            specifying the type of bound(s) to be calculated.

            Defaults to \code{c("double")}.
        }
        \item{\code{conf.what}}{
            A vector of class \code{"character"} describing for which entities
            that confidence should be calculated.

            Defaults to \code{c("blives")}.
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
            Logical value controlling the inclusion of B-life confidence bounds
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
%%        \item{\code{is.plot.ppvalues}}{
%%            Logical value controlling the plotting of the text labels
%%            (the coordinates on the plot) next to the plot positions.
%%
%%            Defaults to \code{FALSE}.
%%        }
        \item{\code{label}}{
            A character string with the title of the legend box, can be used to
            label a dataset or fit.

            Defaults to \code{NULL}.
        }
        \item{\code{legend.text.size}}{
            A number determining the relative legend text size.

            Defaults to \code{0.7}.
        }
        \item{\code{log}}{
            A character string describing the type of plotting canvas to use.
            Possible values are \code{"x"} (Weibull canvas), \code{"xy"}
            (Lognormal canvas), \code{""} (Weibull Y axis, linear X axis)
            and \code{"y"} (Lognormal Y axis, linear X axis).


            Defaults to \code{"x"}, Weibull canvas.
            For contour plots, \code{log} is always \code{""} and cannot be set.
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
        \item{\code{method.conf.blives}}{
            A vector of class \code{"character"} describing the techniques used
            for calculating confidence for B-lives.
            See \code{\link{wblr.conf}} for in-depth discussion of the settings.

            Defaults to \code{c("mcpivotals")}.
        }
        \item{\code{method.fit}}{
            A vector of class \code{"character"} with fitting options.
            See \code{\link{wblr.fit}} for in-depth discussion of the settings.

            Defaults to \code{c("rr","xony")}.
        }
        \item{\code{pch}}{
            An integer or single character describing the plotting symbol, used for
            plotting the datapoints.  For more info, see
            \code{\link[graphics]{points}}.

            Defaults to \code{1}.
        }
        \item{\code{persistent}}{
            Experimental parameter, do not use.

            Defaults to \code{TRUE}.
        }
        \item{\code{pivotals}}{
            Currently unused. Defaults to \code{FALSE}.
        }
        \item{\code{pp}}{
            \code{pp} is short for "probability plot positions", it is a vector of
            class \code{"character"} describing the type of plot positions to be
            calculated.
            Possible values are \code{"median"} (exact median plot postions
            calculated using the Incomplete Beta function) and  \code{"benard"}
            (a good approximation to the median plot positions). Both options
            can be supplied at the same time but currently, only the first
            element in the \code{pp} vector will be used for
            further fitting by \code{\link{wblr.fit}}

            Currently, there are no plans to implement
            \code{"hazen"} and \code{"mean"} plotting positions.

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
        \item{\code{threshold}}{
            A logical value specifying whether a threshold time value should be
            substracted from the graphical elements before plotting. This is
            particulary useful when fitting three parameter models like
            \code{weibull3p}, where the third model parameter is a threshold value.

            If \code{TRUE}, the software will subtract the threshold value
            (that is to be calculated by \code{wblr.fit}) from its
            associated graphical
            elements like plot poistions, fits and confidence bounds.

            If a numeric value is passed, the software will subtract this value
            from the appropriate graphical elements.

            Note: currently, \code{"threshold"} should only be used in combination with
            \code{wblr} and \code{wblr.fit} and not with
            \code{plot.wblr}. See \code{wblr.fit} for more info.


            Defaults to \code{FALSE}
        }
        \item{\code{unrel}}{
            An unordered numeric vector with unreliability levels for which
            B-life confidence will be calculated.

            Defaults to \code{c(0.1,0.05,0.01)}.
        }
        \item{\code{unrel.n}}{
            An integer controlling the amount of unreliability levels for
            which confidence bounds are calculated and ultimately plotted.
            See \code{\link{wblr.conf}} for more details.

            Defaults to \code{25}
        }
        \item{\code{verbosity}}{
            An integer from \code{c(1,2,3)} setting the level of verbosity of
            all \pkg{wblr} functions. Use this for debugging.

            Defaults to \code{1}.
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
   containing the currently active global options of the \pkg{wblr} package.

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
\section{To Do}{
    \itemize{
        \item Add the possibility to reset to the default values.
        \item Implement returning of the pivotals dataframe in the
            \code{wblr} object.
        \item Add detection of unsupported arguments.
    }
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{wblr}},
    \code{\link{wblr.fit}},
    \code{\link{wblr.conf}},
    \code{\link{plot.wblr}}
}
\examples{
## backup options ##
wblr.defaults <- options.wblr()

## setting new options ##
options.wblr(S=5e5,cl=0.99)
%%options.wblr(sub="Testing options.wblr()")

## listing options ##
options.wblr()
options.wblr()$main

## restore options ##
options.wblr(wblr.defaults)
}
