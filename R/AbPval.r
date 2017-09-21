AbPval<-function(
F,R2,model="weibull")  {
m<-1
if(model=="lnorm") {m<-0}

	.Call( "CallgetPvalue", F, R2, m, PACKAGE = "WeibullR" )
}