// src/sgarch_filter.c
#include <R.h>
#include <Rinternals.h>
#include <math.h>

static inline double max2(double a, double b) { return (a > b) ? a : b; }

SEXP sGARCH_filter_c(SEXP x_, SEXP omega_, SEXP alpha_, SEXP beta_,
                     SEXP eps0_, SEXP sigma20_) {

  // ---- Extract and check inputs ------------------------------------------------
  if (!isReal(x_) || LENGTH(x_) < 1) error("x must be a non-empty numeric vector.");
  if (!isReal(omega_) || LENGTH(omega_) != 1) error("omega must be a length-1 numeric.");
  if (!isReal(alpha_)) error("alpha must be numeric.");
  if (!isReal(beta_))  error("beta must be numeric.");

  const double *x      = REAL(x_);
  const double  omega  = REAL(omega_)[0];
  const double *alpha  = REAL(alpha_);
  const double *beta   = REAL(beta_);

  const R_xlen_t n_x = XLENGTH(x_);
  const R_xlen_t p   = XLENGTH(alpha_);
  const R_xlen_t q   = XLENGTH(beta_);

  if (p < 0 || q < 0) error("Internal error: negative p or q.");

  // Optional initials: allow NULL / R_NilValue => no prepend
  int have_init = (eps0_ != R_NilValue) && (sigma20_ != R_NilValue);
  const double *eps0 = NULL;
  const double *s20  = NULL;
  R_xlen_t n_eps0 = 0, n_s20 = 0;

  if (have_init) {
    if (!isReal(eps0_) || !isReal(sigma20_))
      error("eps0 and sigma20 must be numeric or NULL.");
    eps0   = REAL(eps0_);
    s20    = REAL(sigma20_);
    n_eps0 = XLENGTH(eps0_);
    n_s20  = XLENGTH(sigma20_);
    if (n_eps0 != n_s20)
      error("Length mismatch: length(eps0) must equal length(sigma20).");
  }

  // ---- Setup lengths and lags ---------------------------------------------------
  const R_xlen_t k = (R_xlen_t)max2((double)p, (double)q);

  // Final working length (possibly prepended with initials)
  const R_xlen_t n = n_x + (have_init ? n_eps0 : 0);

  // ---- Build eps vector: [eps0, x] or just x -----------------------------------
  SEXP epsS = PROTECT(allocVector(REALSXP, n));
  double *eps = REAL(epsS);

  if (have_init) {
    for (R_xlen_t i = 0; i < n_eps0; ++i) eps[i] = eps0[i];
    for (R_xlen_t i = 0; i < n_x;    ++i) eps[n_eps0 + i] = x[i];
  } else {
    for (R_xlen_t i = 0; i < n_x; ++i) eps[i] = x[i];
  }

  // ---- Long-run variance and sigma2 init ---------------------------------------
  // sigma2_inf = omega / (1 - sum(alpha) - sum(beta))
  double sum_a = 0.0, sum_b = 0.0;
  for (R_xlen_t i = 0; i < p; ++i) sum_a += alpha[i];
  for (R_xlen_t j = 0; j < q; ++j) sum_b += beta[j];
  const double denom = 1.0 - (sum_a + sum_b);
  if (denom <= 0.0) {
    UNPROTECT(1);
    error("Non-stationary parameters: 1 - sum(alpha) - sum(beta) must be > 0.");
  }
  const double sigma2_inf = omega / denom;

  // sigma2: if initials provided, prepend sigma20; otherwise fill with sigma2_inf
  SEXP sigma2S = PROTECT(allocVector(REALSXP, n));
  double *sigma2 = REAL(sigma2S);

  if (have_init) {
    for (R_xlen_t i = 0; i < n_s20; ++i) sigma2[i] = s20[i];
    for (R_xlen_t i = n_s20; i < n;   ++i) sigma2[i] = sigma2_inf;
  } else {
    for (R_xlen_t i = 0; i < n; ++i) sigma2[i] = sigma2_inf;
  }

  // ---- Recursion: for t = k .. n-1 (0-based) -----------------------------------
  // sigma2[t] = omega + sum_{i=1..p} alpha[i-1]*eps[t-i]^2 + sum_{j=1..q} beta[j-1]*sigma2[t-j]
  if (n > 0 && (p > 0 || q > 0)) {
    for (R_xlen_t t = k; t < n; ++t) {
      double arch = 0.0, garch = 0.0;
      for (R_xlen_t i = 1; i <= p; ++i) {
        const double e = eps[t - i];
        arch += alpha[i - 1] * e * e;
      }
      for (R_xlen_t j = 1; j <= q; ++j) {
        garch += beta[j - 1] * sigma2[t - j];
      }
      sigma2[t] = omega + arch + garch;
    }
  }

  // Return sigma2 (same length as working eps)
  UNPROTECT(2);
  return sigma2S;
}
