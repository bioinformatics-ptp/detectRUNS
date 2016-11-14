// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// genoConvertCpp
IntegerVector genoConvertCpp(IntegerVector genotype);
RcppExport SEXP detectRUNS_genoConvertCpp(SEXP genotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type genotype(genotypeSEXP);
    __result = Rcpp::wrap(genoConvertCpp(genotype));
    return __result;
END_RCPP
}
// snpInRunCpp
LogicalVector snpInRunCpp(LogicalVector RunVector, const int window, const float threshold);
RcppExport SEXP detectRUNS_snpInRunCpp(SEXP RunVectorSEXP, SEXP windowSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< LogicalVector >::type RunVector(RunVectorSEXP);
    Rcpp::traits::input_parameter< const int >::type window(windowSEXP);
    Rcpp::traits::input_parameter< const float >::type threshold(thresholdSEXP);
    __result = Rcpp::wrap(snpInRunCpp(RunVector, window, threshold));
    return __result;
END_RCPP
}
