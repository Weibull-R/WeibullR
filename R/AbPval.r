AbPval<-function(F,R2,model="weibull")  {
m<-1
if(tolower(model) %in% c("lnorm","lognormal","lognormal2p")) m<-0
##	.Call( "CallgetPvalue", F, R2, m, PACKAGE = "WeibullR" )
	.Call( CallgetPvalue, F, R2, m)
}