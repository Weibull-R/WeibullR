\name{xbounds}
\alias{xbounds}

\title{
    Extract a bounds dataframe from a \code{wblr} Object
}
\description{
    Creates a list of bounds and ,if appropriate, t0 modifier from a \code{wblr} object.
}
\usage{xbounds(obj, fit_num=1, conf_num)}

\arguments{
    \item{obj}{
        An object created by wblr and the result of at least one modifying wblr.fit call.
	}
	\item{fit_num}{
		The numbered fit in the wblr object.
	}
	\item{conf_num}{
		The numbered confidence calculation in the specified fit of the wblr object.
	}
}

\value{
	A list holding the bounds dataframe and, if modified.by.t0, the t0 value as extracted from the \code{wblr} object.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
  
  Jurgen Symynck, Filip De Bal, Weibull analysis using R, in a nutshell (New Technologies and Products in Machine Manufacturing Technology, Stefan cel Mare University of Suceava, 2010).
}

\examples{
	obj<-wblr.fit(wblr(rweibull(10,1,1)))
	obj<-wblr.conf(obj)
	xbounds(obj)
}
