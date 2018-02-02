\name{mleframe}
\alias{mleframe}

\title{Set life(time) data into the format required by mlefit}

\description{
Validate input data by types and build a dataframe to be used as the primary argument object to mlefit
}

\usage{
mleframe(x, s=NULL, interval=NULL) 
}

\arguments{
  \item{x}{An expected vector of failure data. Alternative support is provided for a dataframe holding time and event columns,
  where event markers for failure occurrences are 1 and suspensions (right censored data) are 0.}
  \item{s}{An optional vector of suspension data.}
  \item{interval}{A dataframe holding interval bounds for failure data in columns named 'left' and 'right'.
  The left(early) interval bound for left-censored data must be entered as zero (NA is not accepted).
  An optional column named 'qty' may contain the integer quantity of data entries having same interval values.}  
}


\value{
Returns a dataframe of the life[time] data in columns named 'left', 'right', and 'qty' with an attribute of fsiq set to TRUE.
  Exact failure data (occurrences) have same time entered in both 'left' and 'right' columns.  
  Suspension data has last known life[time] entered in 'left' column and -1 entered in 'right' column.  
  The interval dataframe argument is appended.  The 'qty' field is populated with a value of 1 for all failure
  and suspension entries and interval rows where qty was not provided. 
}

\references{
  William Q. Meeker and Luis A. Escobar, (1998) "Statistical Methods for Reliability Data", Wiley-Interscience, New York
  Robert B. Abernethy, (2008) "The New Weibull Handbook, Fifth Edition" 
  John I. McCool, (2012) "Using the Weibull Distribution: Reliability, Modeling and Inference"
}

\examples{
fail<-c(10,40,40,50)
susp<-c(20,60)
left<-c(0,0,0,20,10)
right<-c(30,70,100,80,85)
qty<-c(2,1,1,2,1)
interval_ex<-data.frame(left,right,qty)
input_frame<-mleframe(fail,susp,interval_ex)
## time_event dataframe argument
failDF<-data.frame(time=fail, event=1)
suspDF<-data.frame(time=susp, event=0)
time_event_frame<-rbind(failDF, suspDF)
input_frame2<-mleframe(time_event_frame)
## now a time_event_qty dataframe argument
fa<-c(10,40,50)
fq<-c(1,2,1)
su<-susp
faDF<-data.frame(time=fa, event=rep(1,length(fa)), qty=fq)
## note: data.frame is forgiving about repeated single column entry
suDF<-data.frame(time=su, event=0, qty=1)
time_event_qty<-rbind(faDF, suDF)
input_frame3<-mleframe(time_event_qty)
}

\keyword{ interval }