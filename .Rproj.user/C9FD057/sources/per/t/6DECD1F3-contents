#' QMLE Estimate
#'
#' @rdname solarModel_QMLE
#' @name solarModel_QMLE
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_QMLE <- function(model, maxrestarts = 1, seed = 1, quiet = TRUE){
  # Extract data
  data <- dplyr::filter(model$data, isTrain)
  # Time series
  Yt <- data$Yt
  # Time index
  t = data$n
  # Weights
  w = data$weights
  # Initialize unbounded parameters
  init_params <- solarModel_params_to_star_L1(model)
  init_loglik <- -solarModel_quasiLogLik(init_params, Yt, w, t, neg_loglik = FALSE, per_obs = FALSE, quiet = TRUE)
  # Optimize QMLE
  opt <- optim(init_params, solarModel_quasiLogLik, Yt = Yt, t = t, w = w,
               neg_loglik = TRUE, per_obs = FALSE, quiet = TRUE)
  if (!quiet) print(paste0("Log-lik improved by: ", init_loglik - opt$value))
  # QMLE parameters
  theta_star_qml <- opt$par
  # Random change the initial parameters and choose the best one
  if (maxrestarts > 1) {
    best_loglik <- -solarModel_quasiLogLik(theta_star_qml, Yt, w, t, neg_loglik = FALSE, per_obs = FALSE, quiet = TRUE)
    n.restarts <- 1
    set.seed(seed)
    while(n.restarts < maxrestarts) {
      rand_params <- init_params*runif(length(theta_star_qml))
      if (!quiet) print(paste0("Restarting: ", n.restarts, "/", maxrestarts))
      # Optimize the parameters
      opt <- optim(rand_params, solarModel_quasiLogLik, Yt = Yt, t = t, w = w,
                   neg_loglik = TRUE, per_obs = FALSE, quiet = quiet)
      if (!quiet) message("New-loglik: ", opt$value, " Old: ", best_loglik)
      if (opt$value < best_loglik) {
        if (!quiet) print(paste0("Log-lik improved by: ",  abs(opt$value)-abs(best_loglik)))
        best_loglik <- opt$value
        theta_star_qml <- opt$par
      }
      n.restarts <- n.restarts + 1
    }
  }
  # Jacobian
  J_qml <- solarModel_params_from_star_L1(theta_star_qml)
  # Extract converted parameters
  theta_qml <- J_qml$theta
  J_qml <- J_qml$J
  # Number of observations
  n <- length(Yt)
  # Numerical Hessian matrix at QMLE
  H <- numDeriv::hessian(solarModel_quasiLogLik, theta_star_qml, Yt = Yt, t = data$n,
                         neg_loglik = FALSE, per_obs = FALSE, w = data$weights, quiet = TRUE)
  # Numerical Score at QMLE
  S <- numDeriv::jacobian(solarModel_quasiLogLik, theta_star_qml, Yt = Yt, t = data$n,
                          neg_loglik = FALSE, per_obs = TRUE, w = data$weights, quiet = TRUE)
  # Cross products of the score
  p <- length(theta_star_qml)
  B <- matrix(0, p, p)
  for(i in 1:nrow(S)) {
    B <- B + S[i,] %*% t(S[i,])
  }
  # Information matrix (star)
  I <- -H
  # Var-cov matrix (star)
  V_star <- solve(I)
  # Standard errors (star)
  std.error_star <- sqrt(diag(V_star) / n)
  names(std.error_star) <- names(theta_star_qml)
  # Var-cov matrix (original)
  V_orig <- J_qml %*% V_star %*% t(J_qml)
  # Standard errors (star)
  std.error_orig <- sqrt(diag(V_orig))
  names(std.error_orig) <- names(theta_qml)
  # Sandwitch var-cov matrix (star)
  V_rob_star <- V_star %*% B %*% V_star
  std.error_rob_star <- sqrt(diag(V_rob_star))
  names(std.error_rob_star) <- names(theta_star_qml)
  # Sandwitch var-cov matrix (original)
  V_rob_orig <- J_qml %*% V_rob_star %*% t(J_qml)
  std.error_rob_orig <- sqrt(diag(V_rob_orig))
  names(std.error_rob_orig) <- names(theta_qml)
  # Update model's parameters
  model_upd <- model$clone(TRUE)
  model_upd$update(theta_qml[!(names(theta_qml) == "omega")])
  model_upd$seasonal_model_Yt$update_std.errors(std.error_rob_orig[names(model_upd$seasonal_model_Yt$coefficients)])
  model_upd$seasonal_variance$update_std.errors(std.error_rob_orig[names(model_upd$seasonal_variance$coefficients)])

  std.errors.ARMA <- std.error_rob_orig[names(model_upd$ARMA$coefficients)[-1]]
  if (!purrr::is_empty(std.errors.ARMA)) {
    model_upd$ARMA$update_std.errors(std.errors.ARMA)
  }
  std.errors.GARCH <- std.error_rob_orig[names(model_upd$GARCH$coefficients)]
  if (!purrr::is_empty(std.errors.GARCH)) {
    model_upd$GARCH$update_std.errors(std.errors.GARCH)
  }
  model_upd$filter()
  model_upd$fit_NM_model()
  model_upd$update_moments()
  model_upd$update_logLik()
  return(model_upd)
}
