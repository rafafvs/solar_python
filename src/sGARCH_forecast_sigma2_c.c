// sGARCH_forecast_sigma2.c
#include <R.h>
#include <Rinternals.h>

#define IDX(i,j,nrow) ((i) + (j)*(nrow))

// .Call sGARCH_forecast_sigma2_c(A, b, d, e1, S0, H, m2)
SEXP sGARCH_forecast_sigma2_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_, SEXP S0_, SEXP H_, SEXP m2_) {
  const int d = INTEGER(getAttrib(A_, R_DimSymbol))[0];
  const int dC = INTEGER(getAttrib(A_, R_DimSymbol))[1];
  if (d != dC) error("A must be square (d x d).");
  if (LENGTH(b_) != d || LENGTH(d_) != d || LENGTH(e1_) != d || LENGTH(S0_) != d)
    error("b, d, e1, S0 must have length d.");

  int H = asInteger(H_);
  if (H < 1) error("H must be >= 1.");

  const double *A = REAL(A_);
  const double *b = REAL(b_);
  const double *dv = REAL(d_);
  const double *e1 = REAL(e1_);
  const double *S0 = REAL(S0_);

  const double *m2_raw = REAL(m2_);
  const int m2_len = LENGTH(m2_);

  // outputs
  SEXP e_sigma2_ = PROTECT(allocVector(REALSXP, H));               // [1]
  double *e_sigma2 = REAL(e_sigma2_);
  SEXP E_S0_ = PROTECT(allocVector(VECSXP, H));                    // [2]

  // work
  SEXP E_prev_ = PROTECT(allocVector(REALSXP, d));                 // [3]
  double *E_prev = REAL(E_prev_);
  for (int i=0;i<d;++i) E_prev[i] = S0[i];

  SEXP Fh_ = PROTECT(allocVector(REALSXP, d));                     // [4]
  double *Fh = REAL(Fh_);

  for (int h=0; h<H; ++h) {
    // Fh = A %*% E_prev + d
    for (int i=0;i<d;++i){
      double acc = 0.0;
      for (int j=0;j<d;++j) acc += A[IDX(i,j,d)] * E_prev[j];
      Fh[i] = acc + dv[i];
    }
    // e_sigma2[h] = t(e1) %*% Fh
    double Es2 = 0.0;
    for (int i=0;i<d;++i) Es2 += e1[i] * Fh[i];
    e_sigma2[h] = Es2;

    // Eh = Fh + b * (m2[h] * e_sigma2[h])
    const double m2h = (m2_len==1) ? m2_raw[0] : m2_raw[h];
    const double gamma = m2h * Es2;

    SEXP Eh_ = PROTECT(allocVector(REALSXP, d));                    // [loop+1]
    double *Eh = REAL(Eh_);
    for (int i=0;i<d;++i) Eh[i] = Fh[i] + b[i] * gamma;

    SET_VECTOR_ELT(E_S0_, h, Eh_);
    UNPROTECT(1); // Eh_ now owned by list
    for (int i=0;i<d;++i) E_prev[i] = Eh[i];
  }

  SEXP out = PROTECT(allocVector(VECSXP, 2));                       // [5]
  SET_VECTOR_ELT(out, 0, e_sigma2_);
  SET_VECTOR_ELT(out, 1, E_S0_);
  SEXP nms = PROTECT(allocVector(STRSXP, 2));                       // [6]
  SET_STRING_ELT(nms, 0, mkChar("e_sigma2"));
  SET_STRING_ELT(nms, 1, mkChar("E_S0"));
  setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(6); // e_sigma2_, E_S0_, E_prev_, Fh_, out, nms
  return out;
}
