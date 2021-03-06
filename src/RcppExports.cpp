// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// forwardsBackwards
List forwardsBackwards(arma::vec& prior, arma::mat& transmat, arma::mat& f_tk);
RcppExport SEXP _HMMR_forwardsBackwards(SEXP priorSEXP, SEXP transmatSEXP, SEXP f_tkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type transmat(transmatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type f_tk(f_tkSEXP);
    rcpp_result_gen = Rcpp::wrap(forwardsBackwards(prior, transmat, f_tk));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HMMR_forwardsBackwards", (DL_FUNC) &_HMMR_forwardsBackwards, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_HMMR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
