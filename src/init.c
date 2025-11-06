// src/init.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// forward declaration of your .Call function
SEXP pow_matrix_c(SEXP x, SEXP nSEXP);
extern SEXP ARMA_forecast_c(SEXP A_, SEXP X0_, SEXP b_, SEXP h_, SEXP intercept_);
extern SEXP ARMA_filter_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"pow_matrix_c", (DL_FUNC) &pow_matrix_c, 2},
  {"ARMA_forecast_c", (DL_FUNC) &ARMA_forecast_c,   5},
  {"ARMA_filter_c",  (DL_FUNC) &ARMA_filter_c,   6},
  {NULL, NULL, 0}
};

void R_init_yourpkg(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);  // no dynamic symbol lookup
}
