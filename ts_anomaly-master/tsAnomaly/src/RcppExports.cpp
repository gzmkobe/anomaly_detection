// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// segment
NumericVector segment(NumericVector Y, int K, double alpha, bool exact);
RcppExport SEXP tsAnomaly_segment(SEXP YSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP exactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    rcpp_result_gen = Rcpp::wrap(segment(Y, K, alpha, exact));
    return rcpp_result_gen;
END_RCPP
}
