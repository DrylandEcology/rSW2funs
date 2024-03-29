// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GISSM_germination_wait_times
IntegerVector GISSM_germination_wait_times(const IntegerVector& time_to_germinate, const IntegerVector& duration_fave_cond);
RcppExport SEXP _rSW2funs_GISSM_germination_wait_times(SEXP time_to_germinateSEXP, SEXP duration_fave_condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type time_to_germinate(time_to_germinateSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type duration_fave_cond(duration_fave_condSEXP);
    rcpp_result_gen = Rcpp::wrap(GISSM_germination_wait_times(time_to_germinate, duration_fave_cond));
    return rcpp_result_gen;
END_RCPP
}
// GISSM_get_KilledBySoilLayers
LogicalVector GISSM_get_KilledBySoilLayers(const IntegerVector& relevantLayers, const LogicalMatrix& kill_conditions);
RcppExport SEXP _rSW2funs_GISSM_get_KilledBySoilLayers(SEXP relevantLayersSEXP, SEXP kill_conditionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type relevantLayers(relevantLayersSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type kill_conditions(kill_conditionsSEXP);
    rcpp_result_gen = Rcpp::wrap(GISSM_get_KilledBySoilLayers(relevantLayers, kill_conditions));
    return rcpp_result_gen;
END_RCPP
}
// GISSM_kill_seedling
LogicalVector GISSM_kill_seedling(LogicalVector& ss1s, const IntegerVector& ry_year_day, const IntegerVector& ry_useyrs, int y, int doy);
RcppExport SEXP _rSW2funs_GISSM_kill_seedling(SEXP ss1sSEXP, SEXP ry_year_daySEXP, SEXP ry_useyrsSEXP, SEXP ySEXP, SEXP doySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector& >::type ss1s(ss1sSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ry_year_day(ry_year_daySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ry_useyrs(ry_useyrsSEXP);
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type doy(doySEXP);
    rcpp_result_gen = Rcpp::wrap(GISSM_kill_seedling(ss1s, ry_year_day, ry_useyrs, y, doy));
    return rcpp_result_gen;
END_RCPP
}
