// src/ARMA_filter.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

/* y <- A %*% x  (A is m x m, R column-major) */
static inline void dgemv_Ax(int m, const double *A, const double *x, double *y) {
  const double alpha = 1.0, beta = 0.0;
  BLAS_INT M = (BLAS_INT)m, inc = 1;
  char trans = 'N';
  F77_CALL(dgemv)(&trans, &M, &M,
           &alpha, A, &M,
           x, &inc,
           &beta, y, &inc
                    FCONE);
}

/* .Call entry: ARMA_filter_c(A, b, x, p, q, intercept) */
SEXP ARMA_filter_c(SEXP A_, SEXP b_, SEXP x_, SEXP p_, SEXP q_, SEXP intercept_) {
  if (!isReal(A_) || !isReal(b_) || !isReal(x_))
    error("A, b, x must be type 'double'.");

  // dims
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (isNull(dimA) || LENGTH(dimA) != 2) error("A must be a square matrix.");
  int m = INTEGER(dimA)[0];
  int nA1 = INTEGER(dimA)[1];
  if (m != nA1) error("A must be square.");

  int p = asInteger(p_);
  int q = asInteger(q_);
  if (p < 0 || q < 0) error("p and q must be >= 0.");
  if (p + q != m) error("nrow(A) must equal p + q.");
  if (LENGTH(b_) != m) error("length(b) must equal p + q.");

  R_xlen_t n = XLENGTH(x_);
  if (n == 0) error("x must have positive length.");

  double intercept = asReal(intercept_);

  const double *A = REAL(A_);
  const double *b = REAL(b_);
  const double *x = REAL(x_);

  // Output x_hat (copy x and overwrite from i = k..n-1)
  SEXP xhat_ = PROTECT(allocVector(REALSXP, n));        // (1)
  double *xhat = REAL(xhat_);
  for (R_xlen_t i = 0; i < n; ++i) xhat[i] = x[i];

  int k = (p > q) ? p : q;
  if (m == 0) { UNPROTECT(1); return xhat_; } // nothing to filter (p=q=0)

  if (n <= k) { UNPROTECT(1); return xhat_; } // no loop iterations

  // State vector x_t of length m = p+q
  double *xt   = (double*) R_alloc(m, sizeof(double));
  double *next = (double*) R_alloc(m, sizeof(double)); // holds A %*% xt

  // Build initial state:
  // xt[0:(p-1)] = x[k-1], x[k-2], ..., x[k-p]
  // xt[p:(p+q-1)] = 0
  for (int i = 0; i < m; ++i) xt[i] = 0.0;
  for (int j = 0; j < p; ++j) {
    R_xlen_t idx = (R_xlen_t)(k - 1 - j);
    if (idx < 0 || idx >= n) error("Initial state indexing out of bounds.");
    xt[j] = x[idx];
  }

  // Intercept vector c_ = [intercept, 0, ..., 0]^T
  double *cvec = (double*) R_alloc(m, sizeof(double));
  for (int i = 0; i < m; ++i) cvec[i] = 0.0;
  cvec[0] = intercept;

  // Main loop: i in (k+1):n  (R 1-based)  => C: i = k .. n-1
  for (R_xlen_t i = k; i < n; ++i) {
    // xt <- A %*% xt + c_
    dgemv_Ax(m, A, xt, next);
    for (int j = 0; j < m; ++j) next[j] += cvec[j];

    // fitted value
    double xhat_i = next[0];
    xhat[i] = xhat_i;

    // residual
    double e_i = x[i] - xhat_i;

    // update state: xt <- next + b * e_i
    for (int j = 0; j < m; ++j) xt[j] = next[j] + b[j] * e_i;
  }

  UNPROTECT(1);
  return xhat_;
}
