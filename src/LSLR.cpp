#include "WeibullR.h"
#include "LSLRmodel.h"
#include <math.h>

SEXP LSLR(SEXP arg1) {
// Construct an LSLRmodel with this input argument
	std::unique_ptr<LSLRmodel> LM(new LSLRmodel(arg1));
	return Rcpp::wrap(LM->LSLRfit());
}