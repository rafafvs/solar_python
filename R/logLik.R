#' Compute the standard errors and robust standard errors
logLik_std.errors <- function(params, logLik, robust = TRUE){
  # Dimension of the vector of parameters
  k <- length(params)
  # Parameters names
  params_names <- names(params)
  # ***************************************************************
  # Numerical Hessian matrix at QMLE
  H <- numDeriv::hessian(func = logLik, x = params, per_obs = FALSE)
  # Numerical Score at QMLE
  S <- numDeriv::jacobian(func = logLik, x = params, per_obs = TRUE)
  # Cross products of the score
  B <- matrix(0, k, k)
  for(i in 1:nrow(S)) {
    B <- B + S[i,] %*% t(S[i,])
  }
  # ***************************************************************
  # Var-cov matrix
  V <- solve(-H)
  rownames(V) <- colnames(V) <-params_names
  # Standard errors
  if (robust) {
    # Sandwitch var-cov matrix (unbounded)
    V <- V %*% B %*% V
    rownames(V) <- colnames(V) <- params_names
  }
  std.errors <- sqrt(diag(V))
  names(std.errors) <- params_names

  structure(
    list(
      H = H,
      S = S,
      B = B,
      V = V,
      std.errors = std.errors,
      robust = robust
    )
  )
}
