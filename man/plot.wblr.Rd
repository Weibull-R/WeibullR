\name{plot.wblr}
\alias{plot.wblr}

\title{
    \code{wblr} Object Plotting
}
\description{
    This function adds the \code{.wblr} method to \code{\link[graphics]{plot}}
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
    any previously supplied options.
    
    One can pass a list of \code{wblr} objects to \code{\link{plot.wblr}}; in
    that case it is mandatory to use the full method name: \code{plot.wblr(\dots)}
    and not \code{plot(\dots)}.
    
    The calculated Weibull or Lognormal distribution fits are plotted
    on Weibull probability paper by default, but by passing the argument
    \code{\link{log} = "xy"} to the function, Lognormal paper is used
    (see \code{\link{options.wblr}}).
    
    When a \emph{list} of wblr objects is passed, the plot window is generated
    with the options of the first \code{wblr} object in the list.
    % TODO: list these options
        
%    A thin dashed line appears close to the fitted line (if any) when confidence
%    are present and plotted.  
}
\value{
    Currently, the function returns no value. 
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
