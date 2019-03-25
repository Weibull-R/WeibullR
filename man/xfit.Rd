\name{xfit}
\alias{xfit}

\title{
    Extract a Fit Summary from a \code{wblr} Object
}
\description{
    Creates a list of distribution and fit elements from a \code{wblr} object.
}
\usage{xfit(obj, fit_num=1)}

\arguments{
    \item{obj}{
        An object created by wblr and having at least one modifying call to wblr.fit.
	}
	\item{fit_num}{
		The numbered fit in the wblr object.
	}
}

\value{
	A list holding the distribution name and fit vector as extracted from the \code{wblr} object.
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York

  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
  
  Jurgen Symynck, Filip De Bal, Weibull analysis using R, in a nutshell (New Technologies and Products in Machine Manufacturing Technology, Stefan cel Mare University of Suceava, 2010).
}

\examples{
	obj<-wblr.fit(wblr(rweibull(10,1,1)))
	xfit(obj)
}
