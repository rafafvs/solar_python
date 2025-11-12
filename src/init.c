// src/init.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// forward declaration of your .Call function
extern SEXP pow_matrix_c(SEXP x, SEXP nSEXP);
extern SEXP ARMA_forecast_c(SEXP A_, SEXP X0_, SEXP b_, SEXP h_, SEXP intercept_);
extern SEXP ARMA_filter_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sGARCH_filter_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sGARCH_forecast_sigma4_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_, SEXP S0_, SEXP E_S0_, SEXP H_, SEXP m2_, SEXP m4_);
extern SEXP sGARCH_forecast_sigma2_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_, SEXP S0_, SEXP H_, SEXP m2_);
extern SEXP sGARCH_forecast_covariance_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_, SEXP E_S0_, SEXP E2_S0_, SEXP H_, SEXP m2_);
extern SEXP solarModel_quasi_loglik_c(SEXP YtS, SEXP wS, SEXP tS, SEXP aS, SEXP phiS, SEXP thetaS, SEXP bS, SEXP omegaS, SEXP alphaS, SEXP betaS, SEXP negS, SEXP perObsS);


static const R_CallMethodDef CallEntries[] = {
  {"pow_matrix_c", (DL_FUNC) &pow_matrix_c, 2},
  {"ARMA_forecast_c", (DL_FUNC) &ARMA_forecast_c,   5},
  {"ARMA_filter_c",  (DL_FUNC) &ARMA_filter_c,   6},
  {"sGARCH_filter_c", (DL_FUNC) &sGARCH_filter_c, 6},
  {"sGARCH_forecast_sigma2_c", (DL_FUNC) &sGARCH_forecast_sigma2_c, 6},
  {"sGARCH_forecast_sigma4_c", (DL_FUNC) &sGARCH_forecast_sigma4_c, 6},
  {"sGARCH_forecast_covariance_c", (DL_FUNC) &sGARCH_forecast_covariance_c, 6},
  {"solarModel_quasi_loglik_c", (DL_FUNC) &solarModel_quasi_loglik_c, 6},
  {NULL, NULL, 0}
};

void R_init_yourpkg(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);  // no dynamic symbol lookup
}
