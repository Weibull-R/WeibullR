\name{hrbu}

\alias{hrbu}

\title{Hirose and Ross beta unbias factors for Weibull MLE} 

\description{
\code{hrbu} generates the reduction factor based on the mean bias of the weibull MLE beta parameter(roughly C4^6) for complete failure
 samples, modestly increased correction as number of suspensions increases.
}

\usage{
 hrbu(Qx, Qs=NULL)
}

\arguments{
  \item{Qx}{The quantity of actual failures}
  \item{Qs}{An optional quantity of suspensions}
}

\value{
   A factor to be multiplied to the MLE Beta account for known bias.
}

\details{
  This, as many references, discuss the bias reduction in terms of mean.
}

\references{
  Hirose, H. (1999) "Bias Correction for the Maximum Likelihood Estimation in Two-parameter Weibull Distribution" IEEE Transactions on Dielectrics and Electrical Insulation, Vol. 6, No.1
  Ross, R. (1996) "Bias and Standard Deviation Due to Weibull Parameter Estimation for Small Data Sets" IEEE Transactions on Dielectrics and Electrical Insulation, Vol. 3, No.1
}

\examples{
failures<-c(90,96,30,49,82)
suspensions<-c(100,45,10)
MLEfit<-mlefit(mleframe(failures,suspensions))
MLE_Unbiased<-c(MLEfit[1],MLEfit[2]*hrbu(length(failures),length(suspensions)))
}

\keyword{ likelihood bias}