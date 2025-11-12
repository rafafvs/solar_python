#include <R.h>
#include <Rinternals.h>
#include <math.h>

static inline int imax4(int a, int b, int c, int d) {
  int m = (a > b) ? a : b;
  m = (m > c) ? m : c;
  m = (m > d) ? m : d;
  return m;
}

SEXP solarModel_quasi_loglik_c(SEXP YtS, SEXP wS, SEXP tS,
                               SEXP aS, SEXP phiS, SEXP thetaS,
                               SEXP bS, SEXP omegaS, SEXP alphaS, SEXP betaS,
                               SEXP negS, SEXP perObsS)

{
  const double *Yt = REAL(YtS);
  const double *w  = REAL(wS);
  const double *tt = REAL(tS);

  const double *a = REAL(aS);        // length 3: a0, a1, a2
  const double *b = REAL(bS);        // length 3: b0, b1, b2
  const double omega0 = REAL(omegaS)[0];

  const double *phi   = (LENGTH(phiS)   > 0) ? REAL(phiS)   : NULL;
  const double *theta = (LENGTH(thetaS) > 0) ? REAL(thetaS) : NULL;
  const double *alpha = (LENGTH(alphaS) > 0) ? REAL(alphaS) : NULL;
  const double *beta  = (LENGTH(betaS)  > 0) ? REAL(betaS)  : NULL;

  const int n          = LENGTH(YtS);
  const int ar_order   = LENGTH(phiS);
  const int ma_order   = LENGTH(thetaS);
  const int arch_order = LENGTH(alphaS);
  const int garch_order= LENGTH(betaS);

  if (n != LENGTH(wS) || n != LENGTH(tS))
    error("Lengths of Yt, w, and t must match.");

  const int i_star = imax4(ar_order, ma_order, arch_order, garch_order) + 1; // 1-based like R

  const int neg_loglik = asLogical(negS);
  const int per_obs    = asLogical(perObsS);

  // workspace
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, n));
  double *ansp = REAL(ans);
  for (int i = 0; i < n; ++i) ansp[i] = 0.0; // default zeros

  const double LOG2PI = 1.83787706640934533908;
  const double TWO_PI = 6.28318530717958647692;
  const double omega_t = TWO_PI / 365.0;

  // ---- Precompute seasonal quantities ----
  double *sin_t = (double*) R_alloc(n, sizeof(double));
  double *cos_t = (double*) R_alloc(n, sizeof(double));
  double *Ybar  = (double*) R_alloc(n, sizeof(double));
  double *sbar  = (double*) R_alloc(n, sizeof(double));
  double *Ytil  = (double*) R_alloc(n, sizeof(double));
  double *eps   = (double*) R_alloc(n, sizeof(double));
  double *h     = (double*) R_alloc(n, sizeof(double));

  for (int i = 0; i < n; ++i) {
    const double x = omega_t * tt[i];
    sin_t[i] = sin(x);
    cos_t[i] = cos(x);
    Ybar[i]  = a[0] + a[1]*sin_t[i] + a[2]*cos_t[i];
    const double vb = b[0] + b[1]*sin_t[i] + b[2]*cos_t[i];
    if (vb <= 0.0) error("Non-positive seasonal variance at index %d.", i+1);
    sbar[i] = sqrt(vb);
    Ytil[i] = Yt[i] - Ybar[i];
    eps[i]  = 0.0;
    h[i]    = 1.0;
  }

  // ---- Main recursion ----
  for (int j = i_star-1; j < n; ++j) {
    double mu = Ybar[j];

    if (ar_order > 0) {
      for (int k = 1; k <= ar_order; ++k)
        mu += phi[k-1] * Ytil[j-k];
    }
    if (ma_order > 0) {
      for (int k = 1; k <= ma_order; ++k)
        mu += theta[k-1] * eps[j-k];
    }

    eps[j] = Yt[j] - mu;

    double hv = omega0;
    if (arch_order > 0) {
      for (int k = 1; k <= arch_order; ++k) {
        const double z = eps[j-k] / sbar[j-k];
        hv += alpha[k-1] * (z * z);
      }
    }
    if (garch_order > 0) {
      for (int k = 1; k <= garch_order; ++k)
        hv += beta[k-1] * h[j-k];
    }
    if (hv <= 0.0) hv = 1e-12;
    h[j] = hv;

    // ---- log-likelihood ----
    const double sigma_t = sqrt(hv) * sbar[j];
    const double z = eps[j] / sigma_t;
    const double wi = (w[j] > 0.0) ? 1.0 : 0.0;
    double ll = -0.5 * LOG2PI - 0.5 * (z * z) - log(sigma_t);
    ll *= wi;

    ansp[j] = neg_loglik ? -ll : ll;
  }

  // ---- Sum if requested ----
  if (!per_obs) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
      sum += ansp[i];
    SEXP out;
    PROTECT(out = allocVector(REALSXP, 1));
    REAL(out)[0] = sum;
    UNPROTECT(2);
    return out;
  }

  UNPROTECT(1);
  return ans;
}
