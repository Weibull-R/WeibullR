\name{getPPP}
\alias{getPPP}

\title{Determination of percentile plotting positions for linear regression with many optional methods}

\description{
Determination of plotting positions to be used for linear regression of life data analysis is a highly debated topic.  This function implements a wide assortment of options for the Abernethy Reliability Methods project.
}

\usage{
getPPP(x, s=NULL, interval=NULL, ppos="Benard", aranks="Johnson", ties=NULL) 
}

\arguments{
  \item{x}{A vector of failure data.}
  \item{s}{An optional vector of suspension data.}
  \item{interval}{reserved argument for interval censored data - NOT YET IMPLEMENTED.}
  \item{ppos}{A string defining a plotting position method.  Implemented options include "Benard" (as default), "beta" for the incomplete beta function, "mean" also known as Herd-Johnson, "Hazen" or modified Kaplan-Meier, "Kaplan-Meier" with modification for final complete failure, and "Blom.}
  \item{aranks}{A string defining the method for establishing adjusted ranks when suspension data (right censored) are present. Implemented options include "Johnson" (as default) and "KMestimator"}
  \item{ties}{A string defining a method of eliminating ties, or duplicate time valued data, from plotting. Implemented options include "highest" (used as 'Inspection Option #1' by Abernethy), "mean", "lowest", and "sequential".}
}


\value{
Returns a dataframe with the failure data (as potentially reduced by ties argument), the probability plotting positions, and adjusted ranks.
}

\references{
  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition"
  Leonard C. Johnson (1964) "The Statistical Treatment of Fatigue Experiments"  
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
median_percentile_ranks<-getPPP(failures, suspensions)[,2]
}

\keyword{ regression reliability }