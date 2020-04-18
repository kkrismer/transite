// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// score_sequences
SEXP score_sequences(Rcpp::List sequences, Rcpp::NumericMatrix pwm);
RcppExport SEXP _transite_score_sequences(SEXP sequencesSEXP, SEXP pwmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pwm(pwmSEXP);
    rcpp_result_gen = Rcpp::wrap(score_sequences(sequences, pwm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_local_consistency
Rcpp::List calculate_local_consistency(Rcpp::NumericVector x, int numPermutations, int minPermutations, int e);
RcppExport SEXP _transite_calculate_local_consistency(SEXP xSEXP, SEXP numPermutationsSEXP, SEXP minPermutationsSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type numPermutations(numPermutationsSEXP);
    Rcpp::traits::input_parameter< int >::type minPermutations(minPermutationsSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_local_consistency(x, numPermutations, minPermutations, e));
    return rcpp_result_gen;
END_RCPP
}
// calculate_transcript_mc
Rcpp::List calculate_transcript_mc(Rcpp::NumericVector absoluteHits, Rcpp::NumericVector totalSites, double relHitsForeground, int n, int maxPermutations, int minPermutations, int e);
RcppExport SEXP _transite_calculate_transcript_mc(SEXP absoluteHitsSEXP, SEXP totalSitesSEXP, SEXP relHitsForegroundSEXP, SEXP nSEXP, SEXP maxPermutationsSEXP, SEXP minPermutationsSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absoluteHits(absoluteHitsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type totalSites(totalSitesSEXP);
    Rcpp::traits::input_parameter< double >::type relHitsForeground(relHitsForegroundSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type maxPermutations(maxPermutationsSEXP);
    Rcpp::traits::input_parameter< int >::type minPermutations(minPermutationsSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_transcript_mc(absoluteHits, totalSites, relHitsForeground, n, maxPermutations, minPermutations, e));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_transite_score_sequences", (DL_FUNC) &_transite_score_sequences, 2},
    {"_transite_calculate_local_consistency", (DL_FUNC) &_transite_calculate_local_consistency, 4},
    {"_transite_calculate_transcript_mc", (DL_FUNC) &_transite_calculate_transcript_mc, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_transite(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
