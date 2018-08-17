// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// aggregation_resamplingSum
void aggregation_resamplingSum(String originalMapPath, String newMapPath);
RcppExport SEXP _rasteraggregation_aggregation_resamplingSum(SEXP originalMapPathSEXP, SEXP newMapPathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type originalMapPath(originalMapPathSEXP);
    Rcpp::traits::input_parameter< String >::type newMapPath(newMapPathSEXP);
    aggregation_resamplingSum(originalMapPath, newMapPath);
    return R_NilValue;
END_RCPP
}
// aggregation_resamplingAverage
void aggregation_resamplingAverage(String originalMapPath, String newMapPath);
RcppExport SEXP _rasteraggregation_aggregation_resamplingAverage(SEXP originalMapPathSEXP, SEXP newMapPathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type originalMapPath(originalMapPathSEXP);
    Rcpp::traits::input_parameter< String >::type newMapPath(newMapPathSEXP);
    aggregation_resamplingAverage(originalMapPath, newMapPath);
    return R_NilValue;
END_RCPP
}
// aggregation_gdalVersion
String aggregation_gdalVersion();
RcppExport SEXP _rasteraggregation_aggregation_gdalVersion() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(aggregation_gdalVersion());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rasteraggregation_aggregation_resamplingSum", (DL_FUNC) &_rasteraggregation_aggregation_resamplingSum, 2},
    {"_rasteraggregation_aggregation_resamplingAverage", (DL_FUNC) &_rasteraggregation_aggregation_resamplingAverage, 2},
    {"_rasteraggregation_aggregation_gdalVersion", (DL_FUNC) &_rasteraggregation_aggregation_gdalVersion, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_rasteraggregation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
