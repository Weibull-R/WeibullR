getCCC2<-function(F, model="weibull")  {
m<-1
if(model=="lnorm") {m<-0}

	.Call( "CallgetCCC2", F, m, PACKAGE = "WeibullR" )
}