// src/init.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// forward declaration of your .Call function
SEXP pow_matrix_c(SEXP x, SEXP nSEXP);

static const R_CallMethodDef CallEntries[] = {
  {"pow_matrix_c", (DL_FUNC) &pow_matrix_c, 2},
  {NULL, NULL, 0}
};

void R_init_yourpkg(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);  // no dynamic symbol lookup
}
