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
        \item If argument \code{susp} is provided, it is treated as
            a vector of right-censored (life-)time observations (also called
            suspended observations or suspensions).
            In that case, argument \code{time} or \code{fail} is mandatory and is
            treated as a vector of failure times.
        \item If argument \code{x} is of class \code{"data.frame"},
            then it should al least contain \code{$time} and \code{$event}
            columns. Additional columns in the dataframe will be reused in the
            \code{wblr} object, allowing for extra information like
            serial numbers to be included (see section "Examples").
    }

    Note that is is currently allowed to have \code{NA} values in the
    \code{time} argument. In that case, the vector is expected to be
    ordered and no ordering will be applied by \code{wblr}. This feature is
    useful in combination with the output of internal function \code{params.to.ob}.

    \code{wblr} always generates (probability) plot positions for graphically
    displaying the (life-)time observations and for (possible) later usage
    by \code{\link{wblr.fit}}. The type of plot positions is
    controlled with the \code{pp} argument.
    See \code{\link{options.wblr}} for more details.

    Additionally, one can supply any options available from \code{options.wblr},
    such as \code{col} or \code{is.plot.legend}. Some of these options
    will be used when plotting the (life-)time observations using \code{plot.wblr}.
    Subsequent calls to \code{wblr.fit} and \code{wblr.conf} will inherit these options.

    % \code{\link[wblrout:plot.wblr]{plot.wblr}}.
}
\value{
    A named list of class \code{"wblr"}. The first list
    item (\code{$data}) is a dataframe with at least three columns:
    \describe{
        \item{\code{$time}}{
            An ordered vector with (life-)time observations.
        }
        \item{\code{$event}}{
            A vector of class \code{"numeric"} with right-censoring indicators.
            Values of \code{1} mean "dead" or "failed" while \code{0}
            mean "alive" or "right-censored"/"suspended" observations.
            This censoring indicator scheme is modeled after the
            \code{Surv} function of the \pkg{survival}
            package.
            %% add support for factors
        }
        \item{\code{$rank. \dots }}{
            Depending on the argument \code{pp} (either passed as
            an argument or taken from \code{options.wblr}), a vector
            of class \code{"numeric"} with exact median ranks
            (\code{$rank.median}) or Benard's approximate median ranks
            (\code{$rank.benard}).
        }
    }
}
\author{Jurgen Symynck \email{jusy@openreliability.org}}

%%\section{To Do}{
%%    \itemize{
%%        \item Add support for specifying right censored (life-)time observations
%%        using the format
%%
%%        \code{wblr(fail=c(500,1200,900),susp=c(1300,510))}.
%%    }
%%}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{options.wblr}}
}
\examples{
## These code lines all generate the same object ##
wblr(c(500,1200,900,1300,510))
wblr(time=c(500,1200,900,1300,510))
wblr(time=c(500,1200,900,1300,510),event=c(1,1,1,1,1))
wblr(fail=c(500,1200,900,1300,510))
wblr(fail=c(500,1200,900,1300,510),susp=c())
da1 <- data.frame(
    serial=c("S12","S16","S17","S3","S5"),
    time=c(500,1200,900,1300,510),
    event=c(1,1,1,1,1))
wblr(da1,label="complete dataset",pch=1)%da2 <- wblr(da1,label="complete dataset",pch=4,pp="benard",col="red")
da1 <- wblr(da1,label="complete dataset",pch=3,col="orange3")

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
