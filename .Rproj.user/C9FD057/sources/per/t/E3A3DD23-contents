# ************************************************
#                  ARMA functions
# ************************************************
# From unconstraint to contraint parameters
ARMA_zeta_to_phi_jacobian <- function(zeta){
  p <- length(zeta)
  kappa <- tanh(zeta)
  names(kappa) <- paste0("kappa_", 1:length(kappa))

  # store phi^(k) along the way
  phis <- vector("list", p)
  G <- vector("list", p)  # G[[k]] is k x k matrix of d phi^(k) / d kappa[1:k]

  # k = 1
  phis[[1]] <- c(kappa[1])
  G1 <- matrix(0, 1, 1)
  G1[1,1] <- 1
  G[[1]] <- G1

  if (p >= 2) {
    for (k in 2:p) {
      phi_prev <- phis[[k-1]]
      phi_new  <- numeric(k)

      # build new phi^(k)
      for (j in 1:(k-1)) {
        phi_new[j] <- phi_prev[j] - kappa[k] * phi_prev[k-j]
      }
      phi_new[k] <- kappa[k]
      phis[[k]] <- phi_new

      # build new sensitivity matrix G^{(k)} from G^{(k-1)}
      G_prev <- G[[k-1]]         # (k-1) x (k-1)
      Gk <- matrix(0, k, k)      # k x k

      # columns m=1..k-1 propagate; column m=k is special
      # rows j=1..k-1 use the boxed recurrence; row k is [0 ... 0 | 1]
      for (m in 1:(k-1)) {
        for (j in 1:(k-1)) {
          Gk[j, m] <- G_prev[j, m] - kappa[k] * G_prev[k-j, m]
        }
      }
      # column m=k
      for (j in 1:(k-1)) {
        Gk[j, k] <- - phi_prev[k-j]
      }
      Gk[k, 1:(k-1)] <- 0
      Gk[k, k] <- 1

      G[[k]] <- Gk
    }
  }

  phi <- phis[[p]]
  names(phi) <- names(zeta)
  # From phi to kappa
  J_phi_kappa <- G[[p]]
  # From phi to zeta
  sech2 <- 1 - kappa^2   # elementwise
  J_phi_zeta <- J_phi_kappa %*% diag(sech2, nrow = p, ncol = p)
  list(phi = phi,
       kappa = kappa,
       J_phi_kappa = J_phi_kappa,
       J_phi_zeta  = J_phi_zeta)
}
# From constraint to uncontraint parameters
ARMA_phi_to_zeta <- function(phi){
  p <- length(phi)
  if (p == 0) return(numeric(0))
  # work with a list of coefficient vectors phi^(k)
  phik <- vector("list", p)
  phik[[p]] <- phi
  kappa <- numeric(p)
  for (k in seq.int(p, 1)) {
    kappa[k] <- phik[[k]][k]                 # last reflection
    if (abs(kappa[k]) >= 1) {
      # tiny clipping for numerical safety
      kappa[k] <- sign(kappa[k]) * (1 - 1e-12)
    }
    if (k > 1) {
      num <- phik[[k]][1:(k-1)] + kappa[k] * rev(phik[[k]][1:(k-1)])
      den <- 1 - kappa[k]^2
      phik[[k-1]] <- num / den              # back out phi^(k-1)
    }
  }
  coefs_star <- atanh(kappa)
  names(coefs_star) <- names(phi)
  return(coefs_star)
}
# ************************************************
#               Seasonal functions
# ************************************************
# From unconstraint to contraint parameters
seasonal_model_from_b_star_to_b <- function(b_star){
  # Compute the reparametrized coefficients
  b0 <- exp(b_star[1])
  b1 <- b0 * tanh(b_star[2]) * sin(b_star[3])
  b2 <- b0 * tanh(b_star[2]) * cos(b_star[3])
  # Vector of parameters
  b <- c(b0, b1, b2)
  names(b) <- stringr::str_remove_all(names(b_star), "_star")
  return(b)
}
# From constraint to uncontraint parameters
seasonal_model_from_b_to_b_star <- function(b){
  # Compute the reparametrized coefficients
  b0_star <- log(b[1])
  b1_star <- atanh(sqrt(b[2]^2 + b[3]^2)/b[1])
  b2_star <- atan2(b[2], b[3])
  # Vector of parameters
  b_star <- c(b0_star, b1_star, b2_star)
  names(b_star) <- paste0(names(b), "_star")
  return(b_star)
}
# Jacobian
seasonal_model_Jacobian_b_star <- function(b_star){
  # Extract the coefficients
  b0_star <- b_star[1]
  b1_star <- b_star[2]
  b2_star <- b_star[3]
  # Initialize
  J <- matrix(0, 3, 3)
  # Derivatives of b0 wrt b0, b1, b2
  J[1,1] <- exp(b0_star)                                  # d_b0_d_b0 = b0
  J[1,2] <- 0                                             # d_b0_d_b1
  J[1,3] <- 0                                             # d_b0_d_b2
  # Derivatives of b1 wrt b0, b1, b2
  J[2,1] <- J[1,1] * tanh(b1_star) * sin(b2_star)         # d_b1_d_b0 = b1
  J[2,2] <- J[1,1] * (1 - tanh(b1_star)^2) * sin(b2_star) # d_b1_d_b1
  J[2,3] <- J[1,1] * tanh(b1_star) * cos(b2_star)         # d_b1_d_b2 = b2
  # Derivatives of b2 wrt b0, b1, b2
  J[3,1] <- J[2,3]                                        # d_b2_d_b0 = b2
  J[3,2] <- J[1,1] * (1 - tanh(b1_star)^2) * cos(b2_star) # d_b2_d_b1
  J[3,3] <- -J[2,1]                                       # d_b2_d_b2 = -b1
  colnames(J) <- c("d_b0_star", "d_b1_star", "d_b2_star")
  rownames(J) <- c("b0", "b1", "b2")
  return(J)
}
# Gradient with respect to b
seasonal_model_gradient_b <- function(t, seasonal_model){
  # Extract the parameters
  coefs <- seasonal_model$coefficients
  # Extract the periodicity
  omega <- 2*base::pi / seasonal_model$period

  # Compute the gradient with respect to original parameters
  nabla <- matrix(0, length(t), 3)
  nabla[,1] <- 1
  nabla[,2] <- sin(omega * t)
  nabla[,3] <- cos(omega * t)
  colnames(nabla) <- names(coefs)
  rownames(nabla) <- paste0("t = ", as.character(t))
  nabla
}
# Gradient with respect to b_star
seasonal_model_gradient_b_star <- function(t, seasonal_model){
  # Original parameters
  coefs <- seasonal_model$coefficients
  # Reparametrize in terms of b_star
  b_star <- seasonal_model_from_b_to_b_star(coefs)
  # Gradient of original params
  nabla_b <- seasonal_model_gradient_b(t, seasonal_model)
  # Jacobian of b_star
  J_b_star <- seasonal_model_Jacobian_b_star(b_star)
  # Gradient of b_star
  nabla_b_star <- nabla_b
  for(i in 1:nrow(nabla_b)){
    nabla_b_star[i, ] <- t(t(J_b_star) %*% c(nabla_b[i,]))
  }
  colnames(nabla_b_star) <- paste0(colnames(nabla_b), "_star")
  return(nabla_b_star)
}
# ************************************************
#                  GARCH functions
# ************************************************
# From unconstraint to contraint parameters
GARCH_from_zeta_to_phi <- function(zeta){
  # Sigmoid function
  sigmoid <- function(x) 1/(1+exp(-x))
  # Parameters
  alpha <- sigmoid(zeta["alpha_star"])
  beta <- (1 - alpha) * sigmoid(zeta["beta_star"])
  omega <- 1 - alpha - beta
  c(omega = omega[[1]], alpha = alpha[[1]], beta = beta[[1]])
}
# From constraint to uncontraint parameters
GARCH_from_phi_to_zeta <- function(phi){
  # Logit function
  logit <- function(x) log(x/(1-x))
  # Parameters
  alpha_star <- logit(phi["alpha"])
  beta_star <-  logit(phi["beta"]/(1 - phi["alpha"]))
  c(alpha_star = alpha_star[[1]], beta_star = beta_star[[1]])
}
# Jacobian
GARCH_from_zeta_to_phi_jacobian <- function(zeta){
  alpha_star <- zeta["alpha_star"]
  beta_star <- zeta["beta_star"]
  # Compute sigmoids
  s_alpha <- 1/(1+exp(-alpha_star))
  s_beta <- 1/(1+exp(-beta_star))
  # Initialize jacobian
  J <- matrix(0, 3, 2)
  colnames(J) <- c("alpha_star", "beta_star")
  rownames(J) <- c("alpha", "beta", "omega")
  # Derivative of alpha star
  J[1,1] <- s_alpha * (1 - s_alpha)            # d_alpha_star / d_alpha
  J[1,2] <- 0                                  # d_alpha_star / d_beta
  # Derivative of beta star
  J[2,1] <- - s_beta * J[1,1]                  # d_beta_star / d_alpha
  J[2,2] <- - s_alpha * s_beta * (1 - s_beta)  # d_beta_star / d_beta
  # Derivative of omega star
  J[3,1] <- - -(J[1,1] + J[2,1])               # d_omega_star / d_alpha
  J[3,2] <- - -J[2,2]                          # d_omega_star / d_beta
  return(J)
}
# ************************************************
#                  Solar Model
# ************************************************
solarModel_params_to_star_L1 <- function(model){
  # Initialization: seasonal mean
  params <- unlist(model$coefficients$seasonal_model_Yt)
  # AR parameters
  if (model$ARMA$order[1] != 0){
    phi <-  unlist(model$ARMA$phi)
    zeta <- ARMA_phi_to_zeta(phi)
    names(zeta) <- paste0("zeta_", 1:length(zeta))
    params <- c(params, zeta)
  }
  # MA parameters
  if (model$ARMA$order[2] != 0){
    theta <- unlist(model$ARMA$theta)
    xi <- ARMA_phi_to_zeta(theta)
    names(xi) <- paste0("xi_", 1:length(xi))
    params <- c(params, xi)
  }
  # Seasonal variance
  b <- unlist(model$coefficients$seasonal_variance)
  b_star <- seasonal_model_from_b_to_b_star(b)
  names(b_star) <- paste0("b_star_", 1:length(b_star))
  params <- c(params, b_star)
  # GARCH parameters
  if (model$spec$garch_variance){
    if (model$GARCH$order[1] != 0 & model$GARCH$order[2] != 0){
      coefs <- c(model$GARCH$alpha, model$GARCH$beta)
    } else if (model$GARCH$order[1] != 0 & model$GARCH$order[2] == 0){
      coefs <- c(model$GARCH$alpha)
    }
    coefs_star <- ARMA_phi_to_zeta(coefs)
    names(coefs_star) <- paste0(names(coefs), "_star")
    params <- c(params, coefs_star)
  }
  params
}
solarModel_params_from_star_L1 <- function(params){
  # Parameter's names
  params_names <- names(params)
  # Jacobian of the transform
  p <- length(params)
  # Col names
  col_names <- params_names
  # Row Names
  row_names <- params_names
  J_qml <- matrix(0, nrow = p, ncol = p)
  if (sum(stringr::str_detect(params_names, "alpha|beta")) != 0){
    row_names <- c(row_names, "omega_star")
    J_qml <- matrix(0, nrow = p+1, ncol = p)
  }
  # Diagonal entries all equal to 1
  diag(J_qml) <- 1
  colnames(J_qml) <- params_names
  rownames(J_qml) <- row_names
  # Initialization
  #
  # Parameters' vector
  coefs <- c()
  # Seasonal mean
  a <- params[stringr::str_detect(params_names, "a_0|a_sin|a_cos")]
  coefs <- c(coefs, a)
  # AR parameters
  zeta <- params[stringr::str_detect(params_names, "zeta_")]
  if (!purrr::is_empty(zeta)) {
    J_zeta <- ARMA_zeta_to_phi_jacobian(zeta)
    phi <- J_zeta$phi
    names(phi) <- paste0("phi_", 1:length(phi))
    coefs <- c(coefs, phi)
    # Jacobian AR
    idx_row <- stringr::str_detect(row_names, "zeta")
    idx_col <- stringr::str_detect(col_names, "zeta")
    J_qml[idx_row, idx_col] <- J_zeta$J_phi_zeta
  }
  # MA parameters
  xi <- params[stringr::str_detect(params_names, "xi_")]
  if (!purrr::is_empty(xi)) {
    J_xi <- ARMA_zeta_to_phi_jacobian(xi)
    theta <- J_xi$phi
    names(theta) <- paste0("theta_", 1:length(theta))
    coefs <- c(coefs, theta)
    # Jacobian MA
    idx_row <- stringr::str_detect(row_names, "xi")
    idx_col <- stringr::str_detect(col_names, "xi")
    J_qml[idx_row, idx_col] <- J_xi$J_phi_zeta
  }
  # Seasonal variance
  J_b_star <- NA
  b_star <- params[stringr::str_detect(params_names, "b_star_")]
  b <- seasonal_model_from_b_star_to_b(b_star)
  names(b) <- c("c_0", "c_sin_1_365", "c_cos_1_365")
  coefs <- c(coefs, b)
  # Jacobian b_star
  idx_row <- stringr::str_detect(row_names, "b_star")
  idx_col <- stringr::str_detect(col_names, "b_star")
  J_qml[idx_row, idx_col] <- seasonal_model_Jacobian_b_star(b_star)
  # GARCH
  coefs_star <- params[stringr::str_detect(params_names, "alpha[0-9]_star|beta[0-9]_star")]
  if (!purrr::is_empty(coefs_star)) {
    J_garch <- ARMA_zeta_to_phi_jacobian(coefs_star)
    par <- c(J_garch$phi, omega = 1-sum(J_garch$phi))
    names(par) <- stringr::str_remove_all(names(par), "_star")
    coefs <- c(coefs, par)
    # Jacobian GARCH
    idx_row <- stringr::str_detect(row_names, "alpha|beta|omega")
    idx_col <- stringr::str_detect(col_names, "alpha|beta")
    J_qml[idx_row, idx_col] <- rbind(J_garch$J_phi_zeta, 1-rowSums(J_garch$J_phi_zeta))
  }

  structure(
    list(
      theta_star = params,
      theta = coefs,
      J = J_qml
    )
  )
}
# Quasi-likelihood function
solarModel_quasiLogLik <- function(params, Yt, w, t, neg_loglik = FALSE, per_obs = FALSE, quiet = FALSE){
  # Convert uncontraint to contraint parameters
  coefs <- solarModel_params_from_star_L1(params)$theta
  # Parameter's names
  params_names  <- names(coefs)
  # Seasonal mean parameters
  a <- coefs[stringr::str_detect(params_names, "a_0|a_sin|a_cos")]
  # AR parameters
  phi <- coefs[stringr::str_detect(params_names, "phi_")]
  ar_order <- length(phi)
  if (purrr::is_empty(phi)) {
    phi <- 0
  }
  # MA parameters
  theta <- coefs[stringr::str_detect(params_names, "theta_")]
  ma_order <- length(theta)
  if (purrr::is_empty(theta)) {
    theta <- 0
  }
  # Seasonal variance parameters
  b <- coefs[stringr::str_detect(params_names, "c_0|c_sin|c_cos")]
  # ARCH parameters
  alpha <- coefs[stringr::str_detect(params_names, "alpha")]
  arch_order <- length(alpha)
  if (purrr::is_empty(alpha)) {
    alpha <- 0
  }
  # GARCH parameters
  beta <- coefs[stringr::str_detect(params_names, "beta")]
  garch_order <- length(beta)
  if (purrr::is_empty(beta)) {
    beta <- 0
  }
  omega <- (1 - alpha - beta)
  # Number of observations
  n <- length(Yt)
  # Starting point
  i_star <- max(c(ar_order, ma_order, arch_order, garch_order)) + 1
  # Periodicity seasonal functions
  omega_t <- base::pi*2/365
  # Time index
  t <- t
  # Seasonal functions
  sin_t <- sin(omega_t * t)
  cos_t <- cos(omega_t * t)
  # Seasonal mean
  Yt_bar <- a[1] + a[2] * sin_t + a[3] * cos_t
  # Seasonal std. deviation
  sigma_bar <- sqrt(b[1] + b[2] * sin_t + b[3] * cos_t)
  # Deseasonalized time series
  Yt_tilde <- Yt - Yt_bar
  # Initialization
  eps <- rep(0, n)        # ARMA residuals
  h_t <- rep(1, n)        # GARCH variance
  mu_t <- Yt_bar          # Conditional mean
  # ARMA-seasonal-GARCH routine
  for(j in i_star:n){
    # Conditional mean
    if (ar_order != 0) {
      idx_j <- (j-1):(j-ar_order)
      mu_t[j] <- mu_t[j] + (phi %*% Yt_tilde[idx_j])[[1]]
    }
    if (ma_order != 0) {
      idx_j <- (j-1):(j-ma_order)
      mu_t[j] <- mu_t[j] + (theta %*% eps[idx_j])[[1]]
    }
    # ARMA residuals
    eps[j] <- Yt[j] - mu_t[j]
    # Conditional variance
    h_t[j] <- omega
    # ARCH recursion
    if (arch_order != 0) {
      idx_j <- (j-1):(j-arch_order)
      h_t[j] <- h_t[j] + (alpha %*% (eps[idx_j]/sigma_bar[idx_j])^2)[[1]]
    }
    # GARCH recursion
    if (garch_order != 0) {
      idx_j <- (j-1):(j-garch_order)
      h_t[j] <- h_t[j] + (beta %*% h_t[idx_j])[[1]]
    }
  }
  # Conditional variance
  sigma_t <- sqrt(h_t) * sigma_bar
  # Standardized residuals
  eps_tilde <- eps / sigma_t
  # Pseudo likelihoods
  lik <- dnorm(eps_tilde) / sigma_t

  # Pseudo Log-likelihoods
  w <- ifelse(w == 0, 0, 1)
  log_lik <- (w*log(lik))[-c(1:i_star)]
  # Compute total likelihood
  if (!per_obs) {
    log_lik <- sum(log_lik)
  }
  # Return negative likelihood for optimizers
  if (neg_loglik){
    log_lik <- -log_lik
  }
  if(!quiet) print(log_lik)
  return(log_lik)
}


