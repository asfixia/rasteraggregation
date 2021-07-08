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
// aggregation_fieldUniqueValuesString
CharacterVector aggregation_fieldUniqueValuesString(String shapePath, String fieldName);
RcppExport SEXP _rasteraggregation_aggregation_fieldUniqueValuesString(SEXP shapePathSEXP, SEXP fieldNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type shapePath(shapePathSEXP);
    Rcpp::traits::input_parameter< String >::type fieldName(fieldNameSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_fieldUniqueValuesString(shapePath, fieldName));
    return rcpp_result_gen;
END_RCPP
}
// aggregation_sumAreaWhere
double aggregation_sumAreaWhere(String shapePath, String where);
RcppExport SEXP _rasteraggregation_aggregation_sumAreaWhere(SEXP shapePathSEXP, SEXP whereSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type shapePath(shapePathSEXP);
    Rcpp::traits::input_parameter< String >::type where(whereSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_sumAreaWhere(shapePath, where));
    return rcpp_result_gen;
END_RCPP
}
// aggregation_sumAreaByValue
double aggregation_sumAreaByValue(String shapePath, String fieldName, String fieldValue);
RcppExport SEXP _rasteraggregation_aggregation_sumAreaByValue(SEXP shapePathSEXP, SEXP fieldNameSEXP, SEXP fieldValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type shapePath(shapePathSEXP);
    Rcpp::traits::input_parameter< String >::type fieldName(fieldNameSEXP);
    Rcpp::traits::input_parameter< String >::type fieldValue(fieldValueSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_sumAreaByValue(shapePath, fieldName, fieldValue));
    return rcpp_result_gen;
END_RCPP
}
// aggregation_sumAreaGroupedByColumn
DataFrame aggregation_sumAreaGroupedByColumn(String shapePath, String fieldName);
RcppExport SEXP _rasteraggregation_aggregation_sumAreaGroupedByColumn(SEXP shapePathSEXP, SEXP fieldNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type shapePath(shapePathSEXP);
    Rcpp::traits::input_parameter< String >::type fieldName(fieldNameSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_sumAreaGroupedByColumn(shapePath, fieldName));
    return rcpp_result_gen;
END_RCPP
}
// aggregation_init
String aggregation_init();
RcppExport SEXP _rasteraggregation_aggregation_init() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(aggregation_init());
    return rcpp_result_gen;
END_RCPP
}
// aggregation_interpretExpressionDouble
bool aggregation_interpretExpressionDouble(String curFilter, double ruleValue);
RcppExport SEXP _rasteraggregation_aggregation_interpretExpressionDouble(SEXP curFilterSEXP, SEXP ruleValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type curFilter(curFilterSEXP);
    Rcpp::traits::input_parameter< double >::type ruleValue(ruleValueSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_interpretExpressionDouble(curFilter, ruleValue));
    return rcpp_result_gen;
END_RCPP
}
// aggregation_interpretExpressionString
bool aggregation_interpretExpressionString(String curFilter, String ruleValue);
RcppExport SEXP _rasteraggregation_aggregation_interpretExpressionString(SEXP curFilterSEXP, SEXP ruleValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type curFilter(curFilterSEXP);
    Rcpp::traits::input_parameter< String >::type ruleValue(ruleValueSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregation_interpretExpressionString(curFilter, ruleValue));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rasteraggregation_aggregation_resamplingSum", (DL_FUNC) &_rasteraggregation_aggregation_resamplingSum, 2},
    {"_rasteraggregation_aggregation_resamplingAverage", (DL_FUNC) &_rasteraggregation_aggregation_resamplingAverage, 2},
    {"_rasteraggregation_aggregation_fieldUniqueValuesString", (DL_FUNC) &_rasteraggregation_aggregation_fieldUniqueValuesString, 2},
    {"_rasteraggregation_aggregation_sumAreaWhere", (DL_FUNC) &_rasteraggregation_aggregation_sumAreaWhere, 2},
    {"_rasteraggregation_aggregation_sumAreaByValue", (DL_FUNC) &_rasteraggregation_aggregation_sumAreaByValue, 3},
    {"_rasteraggregation_aggregation_sumAreaGroupedByColumn", (DL_FUNC) &_rasteraggregation_aggregation_sumAreaGroupedByColumn, 2},
    {"_rasteraggregation_aggregation_init", (DL_FUNC) &_rasteraggregation_aggregation_init, 0},
    {"_rasteraggregation_aggregation_interpretExpressionDouble", (DL_FUNC) &_rasteraggregation_aggregation_interpretExpressionDouble, 2},
    {"_rasteraggregation_aggregation_interpretExpressionString", (DL_FUNC) &_rasteraggregation_aggregation_interpretExpressionString, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rasteraggregation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
