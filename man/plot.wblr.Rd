\name{plot.wblr}
\alias{plot.wblr}

\title{
    S3 \code{wblr} Object Plotting on pretty canvax
}
\description{
    This function adds the \code{.wblr} method to \code{\link[graphics]{plot.default}}
    from the \pkg{graphics} package.

    Currently, the function plots the (life-)time observations, fits (if any)
    and confidence bounds for B-lives (if any) of an
    \code{wblr} object or a list of \code{wblr} objects on Weibull or
    Lognormal probability paper. 

    For each fit in the (list of) \code{wblr} object(s), legends are added to
    the plot, displaying the fit parameters and (if available) goodness-of-fit
    indicators and confidence information.
}
\usage{\method{plot}{wblr}(x, \dots)}
\arguments{
    \item{x}{
        Object of class \code{"wblr"} or a list of \code{wblr} objects.
    }
    \item{\dots}{
        Options for plotting the \code{wblr} object; see section "Details".
    }
}
\details{
    The \code{\dots} argument can be any graphical parameter that can be
    supplied to \code{\link[graphics]{plot.default}}, and any option that can be
    set by the function \code{\link{options.wblr}}. The options set
    in this way are applied to all graphical elements of the plot, overriding
    any previously supplied options.\cr

    One can pass a list of \code{wblr} objects to \code{plot.wblr}; in
    that case it is mandatory to use the full method name: \code{plot.wblr(\dots)}
    and not \code{plot(\dots)}.\cr

    The calculated Weibull or Lognormal distribution fits are plotted
    on Weibull probability paper by default, but by passing the argument
    \code{canvas = "lognormal"} to the function, Lognormal paper is used.\cr

    When a \emph{list} of wblr objects is passed, the plot window is generated
    with the options of the first \code{wblr} object in the list.\cr

    Graphical controls typically used with \code{plot} or \code{plot.wblr}:
    \describe{
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
		\item{\code{legend.position}}{ default \code{"bottomright"}. See \code{\link{legend}} Details.
		}
		\item{\code{legend.inset}}{ default \code{c(0,0)}. legend inset values are fractions of graph width and height.
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

\value{
    Currently, the function returns no value.
}

\references{
Jurgen Symynck, Filip De Bal, Weibull analysis using R, in a nutshell (New Technologies and Products in Machine Manufacturing Technology, Stefan cel Mare University of Suceava, 2010).

Jurgen Symynck, Filip De Bal, Monte Carlo pivotal confidence bounds for Weibull analysis, with implementations in R (New Technologies and Products in Machine Manufacturing Technology, Stefan cel Mare University of Suceava, 2011).
}

\examples{
options.wblr(blives=0.1) # make the legend boxes a bit shorter...
da2 <- wblr.conf(wblr.fit(wblr(runif(5,10,100),col="red")))
da3 <- wblr.conf(wblr.fit(wblr(rweibull(5,2,1000),col="green4",pch=3)))
da4 <- wblr.conf(wblr.fit(wblr(rlnorm(5,log(500),log(2)),col="blue3",pch=8),
    dist="lognormal"))
% not run due to device dependence
\dontrun{
plot.wblr(list(da2,da3,da4),xlim=c(1,1e6),
    main="Uniformly distributed observations")
}
}
