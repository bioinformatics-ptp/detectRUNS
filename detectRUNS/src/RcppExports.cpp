// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_factor
SEXP fast_factor(SEXP x);
RcppExport SEXP _detectRUNS_fast_factor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_factor(x));
    return rcpp_result_gen;
END_RCPP
}
// genoConvertCpp
IntegerVector genoConvertCpp(IntegerVector genotype);
RcppExport SEXP _detectRUNS_genoConvertCpp(SEXP genotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type genotype(genotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(genoConvertCpp(genotype));
    return rcpp_result_gen;
END_RCPP
}
// pedConvertCpp
IntegerVector pedConvertCpp(CharacterVector genotype);
RcppExport SEXP _detectRUNS_pedConvertCpp(SEXP genotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type genotype(genotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(pedConvertCpp(genotype));
    return rcpp_result_gen;
END_RCPP
}
// homoZygotTestCpp
bool homoZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHet, int maxMiss, int maxGap);
RcppExport SEXP _detectRUNS_homoZygotTestCpp(SEXP xSEXP, SEXP gapsSEXP, SEXP maxHetSEXP, SEXP maxMissSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type maxHet(maxHetSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    rcpp_result_gen = Rcpp::wrap(homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap));
    return rcpp_result_gen;
END_RCPP
}
// heteroZygotTestCpp
bool heteroZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHom, int maxMiss, int maxGap);
RcppExport SEXP _detectRUNS_heteroZygotTestCpp(SEXP xSEXP, SEXP gapsSEXP, SEXP maxHomSEXP, SEXP maxMissSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type maxHom(maxHomSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    rcpp_result_gen = Rcpp::wrap(heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap));
    return rcpp_result_gen;
END_RCPP
}
// findOppositeAndMissing
StringVector findOppositeAndMissing(IntegerVector data, bool ROHet);
RcppExport SEXP _detectRUNS_findOppositeAndMissing(SEXP dataSEXP, SEXP ROHetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type ROHet(ROHetSEXP);
    rcpp_result_gen = Rcpp::wrap(findOppositeAndMissing(data, ROHet));
    return rcpp_result_gen;
END_RCPP
}
// slidingWindowCpp
List slidingWindowCpp(IntegerVector data, IntegerVector gaps, int windowSize, int step, int maxGap, bool ROHet, int maxOppositeGenotype, int maxMiss);
RcppExport SEXP _detectRUNS_slidingWindowCpp(SEXP dataSEXP, SEXP gapsSEXP, SEXP windowSizeSEXP, SEXP stepSEXP, SEXP maxGapSEXP, SEXP ROHetSEXP, SEXP maxOppositeGenotypeSEXP, SEXP maxMissSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type gaps(gapsSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    Rcpp::traits::input_parameter< bool >::type ROHet(ROHetSEXP);
    Rcpp::traits::input_parameter< int >::type maxOppositeGenotype(maxOppositeGenotypeSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    rcpp_result_gen = Rcpp::wrap(slidingWindowCpp(data, gaps, windowSize, step, maxGap, ROHet, maxOppositeGenotype, maxMiss));
    return rcpp_result_gen;
END_RCPP
}
// snpInRunCpp
LogicalVector snpInRunCpp(LogicalVector RunVector, const int windowSize, const float threshold);
RcppExport SEXP _detectRUNS_snpInRunCpp(SEXP RunVectorSEXP, SEXP windowSizeSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type RunVector(RunVectorSEXP);
    Rcpp::traits::input_parameter< const int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< const float >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(snpInRunCpp(RunVector, windowSize, threshold));
    return rcpp_result_gen;
END_RCPP
}
// readPOPCpp
DataFrame readPOPCpp(std::string genotypeFile);
RcppExport SEXP _detectRUNS_readPOPCpp(SEXP genotypeFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type genotypeFile(genotypeFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readPOPCpp(genotypeFile));
    return rcpp_result_gen;
END_RCPP
}
// consecutiveRunsCpp
DataFrame consecutiveRunsCpp(IntegerVector indGeno, List individual, DataFrame mapFile, bool ROHet, int minSNP, int maxOppositeGenotype, int maxMiss, int minLengthBps, int maxGap);
RcppExport SEXP _detectRUNS_consecutiveRunsCpp(SEXP indGenoSEXP, SEXP individualSEXP, SEXP mapFileSEXP, SEXP ROHetSEXP, SEXP minSNPSEXP, SEXP maxOppositeGenotypeSEXP, SEXP maxMissSEXP, SEXP minLengthBpsSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type indGeno(indGenoSEXP);
    Rcpp::traits::input_parameter< List >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type mapFile(mapFileSEXP);
    Rcpp::traits::input_parameter< bool >::type ROHet(ROHetSEXP);
    Rcpp::traits::input_parameter< int >::type minSNP(minSNPSEXP);
    Rcpp::traits::input_parameter< int >::type maxOppositeGenotype(maxOppositeGenotypeSEXP);
    Rcpp::traits::input_parameter< int >::type maxMiss(maxMissSEXP);
    Rcpp::traits::input_parameter< int >::type minLengthBps(minLengthBpsSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    rcpp_result_gen = Rcpp::wrap(consecutiveRunsCpp(indGeno, individual, mapFile, ROHet, minSNP, maxOppositeGenotype, maxMiss, minLengthBps, maxGap));
    return rcpp_result_gen;
END_RCPP
}
// snpInsideRunsCpp
DataFrame snpInsideRunsCpp(DataFrame runsChrom, DataFrame mapChrom, std::string genotypeFile);
RcppExport SEXP _detectRUNS_snpInsideRunsCpp(SEXP runsChromSEXP, SEXP mapChromSEXP, SEXP genotypeFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type runsChrom(runsChromSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type mapChrom(mapChromSEXP);
    Rcpp::traits::input_parameter< std::string >::type genotypeFile(genotypeFileSEXP);
    rcpp_result_gen = Rcpp::wrap(snpInsideRunsCpp(runsChrom, mapChrom, genotypeFile));
    return rcpp_result_gen;
END_RCPP
}
// tableRuns
DataFrame tableRuns(DataFrame runs, std::string genotypeFile, std::string mapFile, const float threshold);
RcppExport SEXP _detectRUNS_tableRuns(SEXP runsSEXP, SEXP genotypeFileSEXP, SEXP mapFileSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type runs(runsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genotypeFile(genotypeFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type mapFile(mapFileSEXP);
    Rcpp::traits::input_parameter< const float >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(tableRuns(runs, genotypeFile, mapFile, threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_detectRUNS_fast_factor", (DL_FUNC) &_detectRUNS_fast_factor, 1},
    {"_detectRUNS_genoConvertCpp", (DL_FUNC) &_detectRUNS_genoConvertCpp, 1},
    {"_detectRUNS_pedConvertCpp", (DL_FUNC) &_detectRUNS_pedConvertCpp, 1},
    {"_detectRUNS_homoZygotTestCpp", (DL_FUNC) &_detectRUNS_homoZygotTestCpp, 5},
    {"_detectRUNS_heteroZygotTestCpp", (DL_FUNC) &_detectRUNS_heteroZygotTestCpp, 5},
    {"_detectRUNS_findOppositeAndMissing", (DL_FUNC) &_detectRUNS_findOppositeAndMissing, 2},
    {"_detectRUNS_slidingWindowCpp", (DL_FUNC) &_detectRUNS_slidingWindowCpp, 8},
    {"_detectRUNS_snpInRunCpp", (DL_FUNC) &_detectRUNS_snpInRunCpp, 3},
    {"_detectRUNS_readPOPCpp", (DL_FUNC) &_detectRUNS_readPOPCpp, 1},
    {"_detectRUNS_consecutiveRunsCpp", (DL_FUNC) &_detectRUNS_consecutiveRunsCpp, 9},
    {"_detectRUNS_snpInsideRunsCpp", (DL_FUNC) &_detectRUNS_snpInsideRunsCpp, 3},
    {"_detectRUNS_tableRuns", (DL_FUNC) &_detectRUNS_tableRuns, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_detectRUNS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
