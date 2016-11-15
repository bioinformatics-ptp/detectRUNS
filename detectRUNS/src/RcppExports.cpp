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
// homoZygotTestCpp
bool homoZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHet, int maxMiss, int maxGap);
RcppExport SEXP detectRUNS_homoZygotTestCpp(SEXP xSEXP, SEXP gapsSEXP, SEXP maxHetSEXP, SEXP maxMissSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type maxHet(maxHetSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    __result = Rcpp::wrap(homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap));
    return __result;
END_RCPP
}
// heteroZygotTestCpp
bool heteroZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHom, int maxMiss, int maxGap);
RcppExport SEXP detectRUNS_heteroZygotTestCpp(SEXP xSEXP, SEXP gapsSEXP, SEXP maxHomSEXP, SEXP maxMissSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type maxHom(maxHomSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    __result = Rcpp::wrap(heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap));
    return __result;
END_RCPP
}
// slidingWindowCpp
LogicalVector slidingWindowCpp(IntegerVector data, IntegerVector gaps, int windowSize, int step, int maxGap, bool ROHet, int maxOppositeGenotype, int maxMiss);
RcppExport SEXP detectRUNS_slidingWindowCpp(SEXP dataSEXP, SEXP gapsSEXP, SEXP windowSizeSEXP, SEXP stepSEXP, SEXP maxGapSEXP, SEXP ROHetSEXP, SEXP maxOppositeGenotypeSEXP, SEXP maxMissSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    Rcpp::traits::input_parameter< bool >::type ROHet(ROHetSEXP);
    Rcpp::traits::input_parameter< int >::type maxOppositeGenotype(maxOppositeGenotypeSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    __result = Rcpp::wrap(slidingWindowCpp(data, gaps, windowSize, step, maxGap, ROHet, maxOppositeGenotype, maxMiss));
    return __result;
END_RCPP
}
// snpInRunCpp
LogicalVector snpInRunCpp(LogicalVector RunVector, const int windowSize, const float threshold);
RcppExport SEXP detectRUNS_snpInRunCpp(SEXP RunVectorSEXP, SEXP windowSizeSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< LogicalVector >::type RunVector(RunVectorSEXP);
    Rcpp::traits::input_parameter< const int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< const float >::type threshold(thresholdSEXP);
    __result = Rcpp::wrap(snpInRunCpp(RunVector, windowSize, threshold));
    return __result;
END_RCPP
}
