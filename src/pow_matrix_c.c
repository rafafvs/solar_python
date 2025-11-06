// src/pow_matrix_c.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>   // BLAS_INT, F77_CALL, FCONE
#include <math.h>   // add this

static void dgemm_nn(const int m, const int k, const int n,
                     const double *A, const double *B, double *C)
{
  const double one = 1.0, zero = 0.0;
  char trans = 'N';
  BLAS_INT M = (BLAS_INT)m, N = (BLAS_INT)n, K = (BLAS_INT)k;
  BLAS_INT lda = M, ldb = K, ldc = M;

  F77_CALL(dgemm)(&trans, &trans,
           &M, &N, &K,
           &one,
           A, &lda,
           B, &ldb,
           &zero,
           C, &ldc
                    /* hidden lengths of the two char* args: */
                    FCONE FCONE);
}

SEXP pow_matrix_c(SEXP x, SEXP nSEXP) {
  if (!isReal(x)) x = PROTECT(coerceVector(x, REALSXP)); else PROTECT(x);
  if (!isInteger(nSEXP)) nSEXP = PROTECT(coerceVector(nSEXP, INTSXP)); else PROTECT(nSEXP);

  SEXP dim = getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || LENGTH(dim) != 2)
    error("x must be a numeric matrix.");

  const int p = INTEGER(dim)[0];
  const int q = INTEGER(dim)[1];
  if (p != q) error("x must be square, got %d x %d.", p, q);

  int n = INTEGER(nSEXP)[0];
  if (n < 0) error("n must be >= 0 (nonnegative powers only).");

  // 1x1 fast path
  if (p == 1) {
    SEXP ans = PROTECT(allocMatrix(REALSXP, 1, 1));
    double base = REAL(x)[0];
    REAL(ans)[0] = (n == 0) ? 1.0 : pow(base, (double)n);
    UNPROTECT(3);
    return ans;
  }

  // n == 0 -> Identity
  if (n == 0) {
    SEXP I = PROTECT(allocMatrix(REALSXP, p, p));
    double *ri = REAL(I);
    for (size_t j = 0; j < (size_t)p*(size_t)p; ++j) ri[j] = 0.0;
    for (int i = 0; i < p; ++i) ri[i + (size_t)i*(size_t)p] = 1.0;
    UNPROTECT(3);
    return I;
  }

  // result := I
  SEXP result = PROTECT(allocMatrix(REALSXP, p, p));
  double *Rres = REAL(result);
  for (size_t j = 0; j < (size_t)p*(size_t)p; ++j) Rres[j] = 0.0;
  for (int i = 0; i < p; ++i) Rres[i + (size_t)i*(size_t)p] = 1.0;

  // base := x, TMP workspace
  double *base = (double*) R_alloc((size_t)p*(size_t)p, sizeof(double));
  memcpy(base, REAL(x), (size_t)p*(size_t)p*sizeof(double));
  double *TMP  = (double*) R_alloc((size_t)p*(size_t)p, sizeof(double));

  // exponentiation by squaring
  while (n > 0) {
    if (n & 1) {                         // result <- result %*% base
      dgemm_nn(p, p, p, Rres, base, TMP);
      memcpy(Rres, TMP, (size_t)p*(size_t)p*sizeof(double));
    }
    n >>= 1;
    if (n > 0) {                         // base <- base %*% base
      dgemm_nn(p, p, p, base, base, TMP);
      memcpy(base, TMP, (size_t)p*(size_t)p*sizeof(double));
    }
  }

  UNPROTECT(3); // x, nSEXP, result
  return result;
}
