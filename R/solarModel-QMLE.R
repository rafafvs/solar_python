#' From unconstrained to constraint parameters
#'
#' @keywords solarModel
#' @noRd
#' @export
solarModel_params_to_zeta <- function(model){

  # Original parameter's names
  params_names <- c()
  # Unconstraint parameters
  params <- c()

  # 1) seasonal mean
  params <- c(params, unlist(model$coefficients$seasonal_model_Yt))
  params_names <- names(params)

  # 2) ARMA parameters
  phi <- model$ARMA$phi
  if (is.null(phi[1])){
    phi <- c(phi_1 = 0)
  }
  theta <- model$ARMA$theta
  if (is.null(theta[1])){
    theta <- c(theta_1 = 0)
  }
  params <- c(params, unlist(purrr::flatten(ARMA_params_to_zeta(phi, theta))))
  params_names <- c(params_names, names(phi), names(theta))

  # 3) Seasonal variance
  b <- unlist(model$coefficients$seasonal_variance)
  b_star <- seasonalModel_params_to_zeta(b)
  names(b_star) <- paste0("b_star_", 1:length(b_star))
  params <- c(params, b_star)
  params_names <- c(params_names, names(b))

  # GARCH parameters
  if (model$spec$garch_variance){
    if (model$GARCH$order[1] != 0 & model$GARCH$order[2] != 0){
      coefs <- c(model$GARCH$alpha, model$GARCH$beta)
    } else if (model$GARCH$order[1] != 0 & model$GARCH$order[2] == 0){
      coefs <- c(model$GARCH$alpha)
    }
    coefs <- c(model$GARCH$omega, coefs)
    coefs_star <- sGARCH_params_to_zeta(coefs, model$GARCH$archOrder, model$GARCH$garchOrder)
    coefs_star <- head(coefs_star, length(coefs_star)-1)
    #names(coefs_star) <- paste0(names(coefs)[-1], "_star")
    params <- c(params, coefs_star)
    params_names <- c(params_names, names(coefs))
  }
  list(
    params = params,
    orig_names = params_names,
    GARCH_order = model$GARCH$order
  )
}

#' From constraint to unconstrained parameters
#'
#' @keywords solarModel
#' @noRd
#' @export
solarModel_params_to_phi <- function(params, orig_names, GARCH_order){
  # Names of the bounded parameters
  params_names <- names(params)
  # Number of parameters
  d <- length(params)
  # Initialize a Jacobian matrix
  J_qmle <- matrix(0, nrow = d+1, ncol = d)
  # Diagonal entries all equal to 1
  diag(J_qmle) <- 1
  # Col names
  col_names <- params_names
  colnames(J_qmle) <- col_names
  # Row Names
  row_names <- orig_names
  rownames(J_qmle) <- row_names
  # Initialization
  # Parameters' vector
  coefs <- c()
  # ************************************************************
  # Seasonal mean
  a <- params[startsWith(params_names, "a_")]
  coefs <- c(coefs, a)
  # ************************************************************
  # AR parameters
  idx_zeta_phi <- startsWith(params_names, "zeta_phi")
  zeta_phi <- params[idx_zeta_phi]
  # MA parameters
  idx_zeta_theta <- startsWith(params_names, "zeta_theta")
  zeta_theta <- params[idx_zeta_theta]
  # Jacobian and ARMA parameters
  J_zeta <- ARMA_params_to_phi(zeta_phi, zeta_theta)
  # AR parameters
  phi <- J_zeta$phi
  # MA parameters
  theta <- J_zeta$theta
  coefs <- c(coefs, phi, theta)
  if (purrr::is_empty(phi)) {
    phi <- c(phi_1 = 0)
  }
  if (purrr::is_empty(theta)) {
    theta <- c(theta_1 = 0)
  }
  # Jacobian ARMA
  idx_row <- stringr::str_detect(row_names, "phi|theta")
  idx_col <- c(idx_zeta_phi | idx_zeta_theta)
  J_qmle[idx_row, idx_col] <- J_zeta$J
  # ************************************************************
  # Seasonal variance
  J_b_star <- NA
  idx_b_star <- startsWith(params_names, "b_star_")
  b_star <- params[idx_b_star]
  b <- seasonalModel_params_to_phi(b_star)
  names(b) <- c("c_0", "c_sin_1_365", "c_cos_1_365")
  coefs <- c(coefs, b)
  # Jacobian b_star
  idx_row <- startsWith(params_names, "c_")
  J_qmle[idx_row, idx_b_star] <- seasonalModel_params_to_zeta_jacobian(b_star)
  # ************************************************************
  # GARCH
  coef_star <- params[stringr::str_detect(params_names, "eta0|kappa[0-9]")]
  omega <- c(omega = 1)
  alpha <- c(alpha1 = 0)
  beta  <- c(beta1 = 0)
  if (!purrr::is_empty(coef_star)) {
    coef_star <- c(coef_star, 0)
    # GARCH parameters
    par <- sGARCH_params_to_phi(coef_star, GARCH_order[1], GARCH_order[2])
    coefs <- c(coefs, par)
    par_names <- names(par)
    omega <- par[startsWith(par_names, "omega")]
    alpha <- par[startsWith(par_names, "alpha")]
    beta  <- par[startsWith(par_names, "beta")]
    # Jacobian
    J_garch <- sGARCH_params_to_zeta_jacobian(coef_star, GARCH_order[1], GARCH_order[2])
    # Jacobian GARCH
    idx_row <- stringr::str_detect(row_names, "alpha|beta|omega")
    idx_col <- stringr::str_detect(col_names, "eta0|kappa[0-9]")
    J_qmle[idx_row, idx_col] <- J_garch
  }
  # ************************************************************
  structure(
    list(
      params = list(a = a, phi = phi, theta = theta, b = b, omega = omega, alpha = alpha, beta = beta),
      theta_star = params,
      theta = coefs,
      J = J_qmle
    )
  )
}

#' Quasi-likelihood function
#'
#' @keywords solarModel
#' @noRd
#' @export
solarModel_quasi_loglik <- function(params, Yt, w, t, orig_names, GARCH_order,
                                   neg_loglik = FALSE, per_obs = FALSE) {

  # Convert uncontraint to contraint parameters
  conv <- solarModel_params_to_phi(params, orig_names, GARCH_order)
  # Inputs
  par <- conv$params

  .Call("solarModel_quasi_loglik_c",
        as.numeric(Yt),
        as.numeric(w),
        as.numeric(t),
        as.numeric(par$a),
        as.numeric(par$phi[par$phi != 0]),
        as.numeric(par$theta[par$theta != 0]),
        as.numeric(par$b),
        as.numeric(par$omega),
        as.numeric(par$alpha[par$alpha != 0]),
        as.numeric(par$beta[par$beta != 0]),
        as.logical(neg_loglik),
        as.logical(per_obs))
}

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
  # Number of observations
  n <- length(Yt)
  # Time index
  t <- data$n
  # Weights
  w <- data$weights
  # Initialize unbounded parameters
  init <- solarModel_params_to_zeta(model)
  init_params <- init$params
  orig_names <- init$orig_names
  GARCH_order <- init$GARCH_order
  # Log-likelihood function
  logLik <- function(params, neg_loglik = FALSE, per_obs = FALSE, quiet = TRUE){
    solarModel_quasi_loglik(params, Yt = Yt, w = w, t = t, orig_names = orig_names,
                           GARCH_order = GARCH_order, neg_loglik = neg_loglik, per_obs = per_obs)
  }
  # Initial log-likelihood
  init_loglik <- logLik(init_params, neg_loglik = TRUE, per_obs = FALSE, quiet = TRUE)
  # QMLE optimization
  opt <- optim(init_params, logLik, neg_loglik = TRUE, per_obs = FALSE, quiet = TRUE)

  # Improvement of the log-likelihood
  if (!quiet) print(paste0("Log-lik improved by: ", init_loglik - opt$value))
  # QMLE parameters
  theta_star_qml <- opt$par
  # Random change the initial parameters and choose the best one
  if (maxrestarts > 1) {
    best_loglik <- logLik(theta_star_qml, Yt, w, t, neg_loglik = TRUE, per_obs = FALSE, quiet = TRUE)
    n.restarts <- 1
    set.seed(seed)
    while(n.restarts < maxrestarts) {
      rand_params <- init_params*runif(length(theta_star_qml))
      if (!quiet) print(paste0("Restarting: ", n.restarts, "/", maxrestarts))
      # Optimize the parameters
      opt <- optim(rand_params, logLik, neg_loglik = TRUE, per_obs = FALSE, quiet = TRUE)
      if (!quiet) message("New-loglik: ", opt$value, " Old: ", best_loglik)
      if (opt$value < best_loglik) {
        if (!quiet) print(paste0("Log-lik improved by: ",  abs(opt$value)-abs(best_loglik)))
        best_loglik <- opt$value
        theta_star_qml <- opt$par
      }
      n.restarts <- n.restarts + 1
    }
  }
  # Constraint parameters
  par_qmle <- solarModel_params_to_phi(theta_star_qml, orig_names, GARCH_order)
  # Jacobian
  J_qmle <- par_qmle$J
  # Numerical Hessian matrix at QMLE
  H <- numDeriv::hessian(logLik, theta_star_qml, neg_loglik = FALSE, per_obs = FALSE, quiet = TRUE)
  # Numerical Score at QMLE
  S <- numDeriv::jacobian(logLik, theta_star_qml, neg_loglik = FALSE, per_obs = TRUE, quiet = TRUE)
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
  # Sandwitch var-cov matrix (star)
  V_rob_star <- V_star %*% B %*% V_star
  # Sandwitch var-cov matrix (original)
  V_rob_orig <- J_qmle %*% V_rob_star %*% t(J_qmle)

  # Extract QMLE parameters
  theta_qml <- par_qmle$theta
  # Robust standard errors
  std.errors <- sqrt(diag(V_rob_orig))
  names(std.errors) <- names(theta_qml)
  # *****************************************************************
  # Update model's parameters
  model_upd <- model$clone(TRUE)
  # Update the parameters
  model_upd$update(theta_qml[!(names(theta_qml) == "omega")])
  # Update standard errors
  suppressMessages(model_upd$seasonal_model_Yt$update_std.errors(std.errors))
  suppressMessages(model_upd$ARMA$update_std.errors(std.errors))
  suppressMessages(model_upd$seasonal_variance$update_std.errors(std.errors))
  suppressMessages(model_upd$GARCH$update_std.errors(std.errors))
  # Filter the data with new parameters
  model_upd$filter()
  # Fit again the Gaussian mixture
  model_upd$fit_NM_model()
  # Update conditional moments
  model_upd$update_moments()
  # Update log-likelihood
  model_upd$update_logLik()

  return(model_upd)
}

