\name{p2y}

\alias{p2y}

\title{Transform Probability Value to the Y-Axis of a \code{"plot.wblr"} Canvas} 

\description{
  This function used to transform a probability value to the y-axis of a plot canvas. 
}

\usage{
 p2y(p,canvas="weibull")
}

\arguments{
  \item{p}{The percentile or probability value}
  \item{canvas}{The name of the wblr canvas to plot to. (Does not have to match the fit distribution. 
  	Use of this transformation permits distributions to appear as curves on unrelated canvas.)
  }
}

\value{
   The transformed y-axis value for a \code{"plot.wblr"} Canvas.
}

\details{
  This applies the inverse Cumulative Distribution function. When applying the default \code{"weibull"} canvas
   this function is equivalent to \code{SPLEDA::qsev}. It is particularly handy when adding points and lines
   to plots generated on wblr objects.
}

