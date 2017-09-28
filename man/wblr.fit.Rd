\name{wblr.fit}
\alias{wblr.fit}

\title{
    Fit Distributions to \code{wblr} Objects
}
\description{
   This function fits probability distributions to \code{wblr} objects.
}
\usage{wblr.fit(x, modify.by.t0=FALSE, \dots)}
\arguments{
    \item{x}{
        Object of class \code{"wblr"}.
    }
	\item{modify.by.t0}{
		A logical value to signifying whether to revise object data by subtraction of the "t0" (translation) parameter of a 3-parameter fit.
		A value of TRUE generates a linearized view of the fit on its base distribution canvas. It is recommended that
		the resulting object have an altered name perhaps adding a ".3p" suffix to the original wblr object to preserve original data.
	}
	\item{\dots}{
        Options for fitting the (life-)time observations, and for plotting the results.
    }
}
\details{
    This function calculates fits for the (life-)time observations in the
    \code{wblr} object and adds them to the object alongside any
    pre-existing fits.

    Fitting options are passed with the \code{dist} and \code{method.fit}
    arguments:
    
    \describe{
        \item{\code{dist}}{
            A character string with the target distribution for fitting.
            Possible values are \code{"weibull"}, \code{"weibull2p"},
            \code{"weibull3p"} (three parameter Weibull), \code{"lognormal"}
            , \code{"lognormal"}or \code{"lognormal3p"}.
           
            Defaults to \code{"weibull"}.
        }
        \item{\code{in.legend}}{
            Logical value controlling the inclusion of various elements in 
            the legend. 
            
            If \code{in.legend=FALSE} is passed,
            the resulting fit calculations will be omitted from the legend, 
            leaving only observation summary data.
            
            Defaults to \code{TRUE}.
        }
        \item{\code{method.fit}}{
            A vector of class \code{"character"} with fitting options.
            
            Defaults to \code{c("rr","xony")}.
            
            
            \itemize{
                \item \code{"rr"}: Rank Regression (RR). Depending on the method for
                calculating probability plot positions chosen during the creation of the
                \code{wblr} object (see option \code{pp} and
                function \code{\link{wblr}}), this can
                either be "exact median rank regression" or
                "Benard's approximate median rank regression".
                If this method is used then it is mandatory to additionally specify
                either X-on-Y or Y-on-X regression.

                \item \code{"xony"},\code{"yonx"}: Differentiate between X-on-Y and Y-on-X
                regression, respectively. For rank regression in lifetime analysis,
                it is best practice to use the X values ((life-)time observations)
                as the response variables whose horizontal distance to the fit line
                must be minimized, and the Y values (unreliabilities) as the
                explanatory variable.

                \item \code{"mle"}: Maximum Likelihood Estimation (MLE), using
                many functions of the \pkg{debias} package.

                \item \code{"mle-rba"}: Maximum Likelihood Estimation with Reduced Bias Adjustment
                as popularized by Abernethy based on the median bias of MLE fitted distributions.
				
				\item \code{"mle-unbias"}: Maximum Likelihood Estimation with bias adjustment
                as popularized by Reliasoft software based on the mean bias of MLE fitted distributions.
			}
        }

        Additionally, one can pass any options available from \code{options.wblr},
        such as \code{col} or \code{is.plot.legend}. The graphical options
        will be used when plotting the (life-)time observations using \code{plot.wblr}.
        Subsequent calls to \code{wblr.conf} will inherit these options.

        The "prr" goodness-of-fit indicator is also calculated
        here. See the "Examples" section on how to retrieve it, see the "References"
        section for additional information on the prr and pve value.
        % \code{\link[wblrout:plot.wblr]{plot.wblr}}.
        
        Passing the \code{threshold} parameter here will result in plotting the
        fit (and its associated plot positions) with a threshold value subtracted. 
        
        If three parameter models like \code{weibull3p} are used and 
        \code{threshold = TRUE}, the calculated third parameter (\code{t0}) of
        the \emph{last three-parameter fit} in the wblr object will be used for 
        plotting the graphical entities. If a numeric value is passed, then 
        only the value from the \emph{last} call to \code{wblr.fit} is used.
        % what if different t0 values are given with different fits? this should
        % result in only one fit used.
        % test the statement that the t0 of the last three-parameter fit
        % is used
        
        Currently, there is no graceful error recovery after attempting to fit
        lifetime data including negative time observations, for example 
        \code{wblr.fit(wblr(-5:10)).}
        
    }
}
\value{
    The function returns its argument object \code{x}, extended with the
    calculated fit and the optional graphical and calculation arguments as
    provided to the function.
}
\author{Jurgen Symynck \email{jusy@openreliability.org}}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
}

