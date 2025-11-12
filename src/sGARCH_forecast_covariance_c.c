// sGARCH_forecast_covariance_fast.c
#include <R.h>
#include <Rinternals.h>
#include <math.h>

#define IDX(i,j,nrow) ((i) + (j)*(nrow))

// .Call sGARCH_forecast_covariance_c(A, b, d, e1, E_S0, E2_S0, H, m2)
SEXP sGARCH_forecast_covariance_c(SEXP A_, SEXP b_, SEXP d_, SEXP e1_,
                                  SEXP E_S0_, SEXP E2_S0_, SEXP H_, SEXP m2_) {
  const int d = INTEGER(getAttrib(A_, R_DimSymbol))[0];
  const int dC= INTEGER(getAttrib(A_, R_DimSymbol))[1];
  if (d != dC) error("A must be square (d x d).");
  if (LENGTH(b_) != d || LENGTH(d_) != d || LENGTH(e1_) != d) error("b,d,e1 must have length d.");
  if (!isNewList(E_S0_) || !isNewList(E2_S0_)) error("E_S0 and E2_S0 must be lists.");
  const int H = asInteger(H_);
  if (H < 1 || LENGTH(E_S0_) != H || LENGTH(E2_S0_) != H) error("H must match list lengths.");

  const double *A  = REAL(A_);
  const double *b  = REAL(b_);
  const double *dv = REAL(d_);
  const double *e1 = REAL(e1_);
  const double *m2 = REAL(m2_);
  const int m2_len = LENGTH(m2_);

  // Precompute scalar projections
  double e1_dot_d = 0.0, e1_dot_b = 0.0;
  for (int i=0;i<d;++i) { e1_dot_d += e1[i]*dv[i]; e1_dot_b += e1[i]*b[i]; }

  // Precompute vectors: a = A^T e1  (length d)
  SEXP a_ = PROTECT(allocVector(REALSXP, d));
  double *a = REAL(a_);
  for (int i=0;i<d;++i){ // a[i] = sum_j A[j,i]*e1[j]
    double acc=0.0; for (int j=0;j<d;++j) acc += A[IDX(j,i,d)] * e1[j];
    a[i] = acc;
  }

  // Precompute e1^T E_h for all h, and diag terms t_{h,h} = e1^T Sigma_h e1
  SEXP e1Eh_ = PROTECT(allocVector(REALSXP, H));
  SEXP diag_t_ = PROTECT(allocVector(REALSXP, H));
  double *e1Eh = REAL(e1Eh_), *diag_t = REAL(diag_t_);
  for (int h=0; h<H; ++h){
    const double *Eh = REAL(VECTOR_ELT(E_S0_, h));
    double e = 0.0; for (int i=0;i<d;++i) e += e1[i]*Eh[i];
    e1Eh[h] = e;
    const double *Sh = REAL(VECTOR_ELT(E2_S0_, h)); // d x d
    double t = 0.0;
    for (int j=0;j<d;++j){
      double col = 0.0; for (int i=0;i<d;++i) col += e1[i]*Sh[IDX(i,j,d)];
      t += col * e1[j];
    }
    diag_t[h] = t;
  }

  // Output covariance
  SEXP C_ = PROTECT(allocMatrix(REALSXP, H, H));
  double *C = REAL(C_);

  // Lower triangle (including diagonal), then mirror
  // Temporary vector x := M_{h,l} e1
  SEXP x_ = PROTECT(allocVector(REALSXP, d));
  double *x = REAL(x_);

  // For each column l, seed from diagonal (h=l) using Sigma_l
  for (int l=0; l<H; ++l){
    // Diagonal
    C[IDX(l,l,H)] = diag_t[l] - e1Eh[l]*e1Eh[l];

    if (l == H-1) continue;

    // Initialize x_{l,l} = Sigma_l * e1  (d×d times e1)
    const double *Sl = REAL(VECTOR_ELT(E2_S0_, l));
    for (int i=0;i<d;++i){
      double acc=0.0; for (int k=0;k<d;++k) acc += Sl[IDX(i,k,d)] * e1[k];
      x[i] = acc;
    }

    // Recurse forward h = l+1..H-1
    const double e1_dot_El = e1Eh[l];
    for (int h=l+1; h<H; ++h){
      // rho = (e1^T A) x_{h-1,l} + (e1^T d) (e1^T E_l)
      double e1TA_x = 0.0; // (e1^T A) x
      for (int k=0;k<d;++k) e1TA_x += a[k] * x[k]; // since a = A^T e1, (e1^T A) x = a^T x
      double rho = e1TA_x + e1_dot_d * e1_dot_El;

      // t_{h,l} = a^T x_{h-1,l} + e1^T d * e1^T E_l + e1^T b * m2[h] * rho
      double t_hl = e1TA_x + e1_dot_d * e1_dot_El + e1_dot_b * ((m2_len==1)? m2[0] : m2[h]) * rho;

      // Cov entry
      C[IDX(h,l,H)] = t_hl - e1Eh[h]*e1Eh[l];

      // Update x: x_{h,l} = A x_{h-1,l} + d (e1^T E_l) + b m2[h] * rho
      // tmp = A x
      double tmp_i;
      for (int i=0;i<d;++i){
        double acc=0.0; for (int k=0;k<d;++k) acc += A[IDX(i,k,d)] * x[k];
        tmp_i = acc + dv[i]*e1_dot_El + b[i]*(((m2_len==1)? m2[0] : m2[h]) * rho);
        x[i] = tmp_i;
      }
    }
  }

  // Symmetrize
  for (int l=0; l<H; ++l)
    for (int h=l+1; h<H; ++h)
      C[IDX(l,h,H)] = C[IDX(h,l,H)];

  UNPROTECT(5); // a_, e1Eh_, diag_t_, C_, x_
  return C_;
}
