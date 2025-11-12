// sGARCH_forecast_sigma4.c
#include <R.h>
#include <Rinternals.h>

#define IDX(i,j,nrow) ((i) + (j)*(nrow))

// .Call sGARCH_forecast_sigma4_c(A, b, d, e1, S0, E_S0, H, m2, m4)
SEXP sGARCH_forecast_sigma4_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_, SEXP S0_, SEXP E_S0_, SEXP H_, SEXP m2_, SEXP m4_) {
  const int d  = INTEGER(getAttrib(A_, R_DimSymbol))[0];
  const int dC = INTEGER(getAttrib(A_, R_DimSymbol))[1];
  if (d != dC) error("A must be square (d x d).");
  if (LENGTH(b_) != d || LENGTH(d_) != d || LENGTH(e1_) != d || LENGTH(S0_) != d) error("b,d,e1,S0 must have length d.");

  if (!isNewList(E_S0_)) error("E_S0 must be a list.");
  int H = asInteger(H_);
  if (H < 1 || LENGTH(E_S0_) != H) error("H must match length(E_S0).");

  const double *A = REAL(A_);
  const double *b = REAL(b_);
  const double *dv = REAL(d_);
  const double *e1 = REAL(e1_);
  const double *S0 = REAL(S0_);
  const double *m2 = REAL(m2_);
  const double *m4 = REAL(m4_);
  const int m2_len = LENGTH(m2_);
  const int m4_len = LENGTH(m4_);

  // dd = d d^T, bb = b b^T
  SEXP dd_ = PROTECT(allocMatrix(REALSXP, d, d));                   // [1]
  double *dd = REAL(dd_);
  SEXP bb_ = PROTECT(allocMatrix(REALSXP, d, d));                   // [2]
  double *bb = REAL(bb_);
  for (int j=0;j<d;++j)
    for (int i=0;i<d;++i) {
      dd[IDX(i,j,d)] = dv[i]*dv[j];
      bb[IDX(i,j,d)] = b[i]*b[j];
    }

    // Sigma_prev = S0 S0^T
    SEXP Sigma_prev_ = PROTECT(allocMatrix(REALSXP, d, d));           // [3]
  double *Sigma_prev = REAL(Sigma_prev_);
  for (int j=0;j<d;++j)
    for (int i=0;i<d;++i)
      Sigma_prev[IDX(i,j,d)] = S0[i]*S0[j];

  // outputs
  SEXP theta_    = PROTECT(allocVector(REALSXP, H));                // [4]
  SEXP e_sigma4_ = PROTECT(allocVector(REALSXP, H));                // [5]
  double *theta    = REAL(theta_);
  double *e_sigma4 = REAL(e_sigma4_);
  SEXP E2_S0_ = PROTECT(allocVector(VECSXP, H));                    // [6]

  // work buffers
  SEXP T_  = PROTECT(allocMatrix(REALSXP, d, d));                   // [7]
  double *T = REAL(T_);
  SEXP Av_ = PROTECT(allocVector(REALSXP, d));                      // [8]
  double *Av = REAL(Av_);
  SEXP tmp_ = PROTECT(allocMatrix(REALSXP, d, d));                  // [9]
  double *tmp = REAL(tmp_);

  for (int h=0; h<H; ++h) {
    const double *E_prev = (h==0) ? S0 : REAL(VECTOR_ELT(E_S0_, h-1));

    // tmp = A %*% Sigma_prev
    for (int j=0;j<d;++j)
      for (int i=0;i<d;++i) {
        double acc=0.0;
        for (int k=0;k<d;++k) acc += A[IDX(i,k,d)] * Sigma_prev[IDX(k,j,d)];
        tmp[IDX(i,j,d)] = acc;
      }
      // T = tmp %*% A^T
      for (int i=0;i<d;++i)
        for (int j=0;j<d;++j) {
          double acc=0.0;
          for (int k=0;k<d;++k) acc += tmp[IDX(i,k,d)] * A[IDX(j,k,d)];
          T[IDX(i,j,d)] = acc;
        }

        // Av = A %*% E_prev
        for (int i=0;i<d;++i){
          double acc=0.0;
          for (int k=0;k<d;++k) acc += A[IDX(i,k,d)] * E_prev[k];
          Av[i] = acc;
        }
        // T += A E_prev d^T + d E_prev^T A^T + d d^T
        for (int j=0;j<d;++j)
          for (int i=0;i<d;++i) {
            T[IDX(i,j,d)] += Av[i]*dv[j] + dv[i]*Av[j] + dd[IDX(i,j,d)];
          }

          // theta[h] = e1^T T e1
          double th = 0.0;
    for (int j=0;j<d;++j){
      double col=0.0;
      for (int i=0;i<d;++i) col += e1[i]*T[IDX(i,j,d)];
      th += col * e1[j];
    }
    theta[h] = th;

    const double m2h = (m2_len==1)? m2[0] : m2[h];
    const double m4h = (m4_len==1)? m4[0] : m4[h];

    // v = T e1, w = e1^T T
    SEXP v_ = PROTECT(allocVector(REALSXP, d));                     // [loop+1]
    double *v = REAL(v_);
    SEXP w_ = PROTECT(allocVector(REALSXP, d));                     // [loop+2]
    double *w = REAL(w_);
    for (int i=0;i<d;++i){
      double acc=0.0; for (int k=0;k<d;++k) acc += T[IDX(i,k,d)] * e1[k];
      v[i] = acc;
    }
    for (int j=0;j<d;++j){
      double acc=0.0; for (int k=0;k<d;++k) acc += e1[k] * T[IDX(k,j,d)];
      w[j] = acc;
    }

    // Sigma_h = T + m2h ( v b^T + b w ) + m4h * th * (b b^T)
    SEXP Sigma_h_ = PROTECT(allocMatrix(REALSXP, d, d));            // [loop+3]
    double *Sigma_h = REAL(Sigma_h_);
    for (int j=0;j<d;++j)
      for (int i=0;i<d;++i) {
        Sigma_h[IDX(i,j,d)] = T[IDX(i,j,d)]
        + m2h * ( v[i]*b[j] + b[i]*w[j] )
        + m4h * th * (b[i]*b[j]);
      }

      // e_sigma4[h] = e1^T Sigma_h e1
      double Es4 = 0.0;
    for (int j=0;j<d;++j){
      double col=0.0;
      for (int i=0;i<d;++i) col += e1[i]*Sigma_h[IDX(i,j,d)];
      Es4 += col * e1[j];
    }
    e_sigma4[h] = Es4;

    SET_VECTOR_ELT(E2_S0_, h, Sigma_h_);
    UNPROTECT(3); // v_, w_, Sigma_h_
    for (int i=0;i<d*d;++i) Sigma_prev[i] = Sigma_h[i];
  }

  SEXP out = PROTECT(allocVector(VECSXP, 3));                       // [10]
  SET_VECTOR_ELT(out, 0, theta_);
  SET_VECTOR_ELT(out, 1, e_sigma4_);
  SET_VECTOR_ELT(out, 2, E2_S0_);
  SEXP nms = PROTECT(allocVector(STRSXP, 3));                       // [11]
  SET_STRING_ELT(nms, 0, mkChar("theta"));
  SET_STRING_ELT(nms, 1, mkChar("e_sigma4"));
  SET_STRING_ELT(nms, 2, mkChar("E2_S0"));
  setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(11); // dd_, bb_, Sigma_prev_, theta_, e_sigma4_, E2_S0_, T_, Av_, tmp_, out, nms
  return out;
}
