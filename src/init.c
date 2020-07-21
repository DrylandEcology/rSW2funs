#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls: Rcpp registers these correctly */
extern SEXP _rSW2funs_GISSM_germination_wait_times(SEXP, SEXP);
extern SEXP _rSW2funs_GISSM_get_KilledBySoilLayers(SEXP, SEXP);
extern SEXP _rSW2funs_GISSM_kill_seedling(SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"_rSW2funs_GISSM_germination_wait_times", (DL_FUNC) &_rSW2funs_GISSM_germination_wait_times, 2},
    {"_rSW2funs_GISSM_get_KilledBySoilLayers", (DL_FUNC) &_rSW2funs_GISSM_get_KilledBySoilLayers, 2},
    {"_rSW2funs_GISSM_kill_seedling", (DL_FUNC) &_rSW2funs_GISSM_kill_seedling, 5},
    {NULL, NULL, 0}
};

void R_init_rSW2funs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
