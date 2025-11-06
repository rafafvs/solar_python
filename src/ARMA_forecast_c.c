// src/arma_forecast.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

// y <- A %*% x  (A is m x m, column-major from R)
static inline void dgemv_Ax(const int m, const double *A, const double *x, double *y) {
  const double alpha = 1.0, beta = 0.0;
  BLAS_INT inc = 1, M = (BLAS_INT)m;
  char trans = 'N';
  F77_CALL(dgemv)(&trans, &M, &M,
           &alpha,
           A, &M,
           x, &inc,
           &beta,
           y, &inc
                    /* hidden length for 'trans' */ FCONE);
}

SEXP ARMA_forecast_c(SEXP A_, SEXP X0_, SEXP b_, SEXP h_, SEXP intercept_) {
  if (!isReal(A_) || !isReal(X0_) || !isReal(b_))
    error("A, X0, b must be numeric (double).");

  // Dimensions (read BEFORE any allocations)
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (isNull(dimA) || LENGTH(dimA) != 2) error("A must be a square matrix.");
  int m = INTEGER(dimA)[0];
  int n = INTEGER(dimA)[1];
  if (m != n) error("A must be square.");
  if (LENGTH(X0_) != m) error("X0 length must equal nrow(A).");
  if (LENGTH(b_)  != m) error("b length must equal nrow(A).");

  int h = asInteger(h_);
  if (h < 1) error("h must be >= 1.");
  double intercept = asReal(intercept_);

  const double *A = REAL(A_);
  const double *X0 = REAL(X0_);
  const double *b  = REAL(b_);

  // Work buffers (no PROTECT needed)
  double *v  = (double*) R_alloc(m, sizeof(double));
  double *nv = (double*) R_alloc(m, sizeof(double));
  double *u  = (double*) R_alloc(m, sizeof(double));
  double *nu = (double*) R_alloc(m, sizeof(double));
  double *w  = (double*) R_alloc(m, sizeof(double));
  double *nw = (double*) R_alloc(m, sizeof(double));

  // Outputs
  SEXP Yt_tilde_hat_   = PROTECT(allocVector(REALSXP, h)); // 1
  SEXP psi_            = PROTECT(allocVector(REALSXP, h)); // 2
  SEXP psi2_           = PROTECT(allocVector(REALSXP, h)); // 3
  SEXP intercept_steps_= PROTECT(allocVector(REALSXP, h)); // 4
  double *Yt_tilde_hat = REAL(Yt_tilde_hat_);
  double *psi          = REAL(psi_);
  double *psi2         = REAL(psi2_);
  double *int_steps    = REAL(intercept_steps_);

  // Initialize iterates
  dgemv_Ax(m, A, X0, v);                 // v = A X0
  for (int i=0;i<m;++i) u[i]=b[i];       // u = b
  for (int i=0;i<m;++i) w[i]=0.0; w[0]=1.0; // w = e1

  // Precompute s_ab[k] = (A^k b)_1, k=0..h-1
  SEXP s_ab_ = PROTECT(allocVector(REALSXP, h));          // 5
  double *s_ab = REAL(s_ab_);
  s_ab[0] = u[0];
  for (int k=1; k<h; ++k) {
    dgemv_Ax(m, A, u, nu);
    double *tmp = u; u = nu; nu = tmp;
    s_ab[k] = u[0];
  }

  // Fill outputs for j=1..h
  for (int j=1; j<=h; ++j) {
    Yt_tilde_hat[j-1] = v[0];

    dgemv_Ax(m, A, w, nw);                    // w <- A w
    double *tmpw = w; w = nw; nw = tmpw;
    int_steps[j-1] = intercept * w[0];

    if (j < h) {                              // v <- A v
      dgemv_Ax(m, A, v, nv);
      double *tmpv = v; v = nv; nv = tmpv;
    }

    psi[j-1]  = s_ab[h - j];
    psi2[j-1] = psi[j-1] * psi[j-1];
  }

  // Return list
  SEXP res   = PROTECT(allocVector(VECSXP, 4));           // 6
  SEXP names = PROTECT(allocVector(STRSXP, 4));           // 7
  SET_STRING_ELT(names, 0, mkChar("Yt_tilde_hat"));
  SET_STRING_ELT(names, 1, mkChar("psi_j"));
  SET_STRING_ELT(names, 2, mkChar("psi2_j"));
  SET_STRING_ELT(names, 3, mkChar("intercept"));
  SET_VECTOR_ELT(res, 0, Yt_tilde_hat_);
  SET_VECTOR_ELT(res, 1, psi_);
  SET_VECTOR_ELT(res, 2, psi2_);
  SET_VECTOR_ELT(res, 3, intercept_steps_);
  setAttrib(res, R_NamesSymbol, names);

  UNPROTECT(7);
  return res;
}
