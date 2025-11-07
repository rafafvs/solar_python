#' Construct the vector b of an ARMA model
#'
#' @param arOrder Numeric scalar, order of AR model.
#' @param maOrder Numeric scalar, order of MA model.
#' @examples
#' # ARMA(2,2)
#' ARMA_vector_b(2,2)
#' # AR(2)
#' ARMA_vector_b(2,0)
#' # MA(2)
#' ARMA_vector_b(0,2)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_vector_b
#' @name ARMA_vector_b
#' @export
#' @noRd
ARMA_vector_b <- function(arOrder, maOrder){
  # Companion vector AR
  e_p <- c()
  if (arOrder > 1){
    e_p <- c(1, rep(0, arOrder-1))
  } else if (arOrder == 1) {
    e_p <- 1
  }
  # Companion vector MA
  e_q <- c()
  if (maOrder > 1){
    e_q <- c(1, rep(0, maOrder-1))
  } else if (maOrder == 1) {
    e_q <- 1
  }
  # Combine the two parts
  b <- matrix(c(e_p, e_q), ncol = 1)
  return(b)
}

#' Construct the companion matrix of an ARMA model
#'
#' @param phi Numeric vector with length `p`, AR parameters.
#' @param theta Numeric vector with length `q`, MA parameters.
#' @return A square matrix with dimension (p+q)
#' @examples
#' # AR(1) / MA(1) ~ No companion
#' ARMA_companion_matrix(c(0.4))
#' ARMA_companion_matrix(theta = c(0.4))
#' # Only AR
#' ARMA_companion_matrix(c(0.4, 0.3, 0.1))
#' # Only MA
#' ARMA_companion_matrix(theta = c(0.4, 0.3, 0.1))
#' # ARMA
#' ARMA_companion_matrix(c(0.4, 0.2), c(0.3))
#' ARMA_companion_matrix(c(0.1, 0.02, 0.01), c(0.3, 0.1, 0.05))
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_companion_matrix
#' @name ARMA_companion_matrix
#' @export
ARMA_companion_matrix <- function(phi, theta){
  if (missing(theta)){
    theta <- NULL
  }
  if (missing(phi)){
    phi <- NULL
  }
  # AR order
  p <- length(phi)
  # MA order
  q <- length(theta)
  I_p <- L_p <- NULL
  if (p > 0){
    # Identity matrix
    I_p <- diag(1, nrow = p-1)
    # Add colum of zeros
    L_p <- cbind(I_p, rep(0, p - 1))
    # Matrix of zeros
    L_p <- cbind(L_p, matrix(0, nrow = p-1, ncol = q))
  }
  I_q <- L_q <- NULL
  if (q > 0) {
    # Identity matrix
    I_q <- diag(1, nrow = q-1)
    # Add colum of zeros
    L_q <- cbind(I_q, rep(0, q-1))
    # Matrix of zeros
    L_q <- cbind(matrix(0, nrow = q-1, ncol = p), L_q)
  }
  zeros <- NULL
  if (p >= 1 & q >= 1) {
    zeros <- rep(0, p + q)
  }
  # Companion matrix
  A <- rbind(c(phi, theta), L_p, zeros, L_q)
  rownames(A) <- colnames(A) <- NULL
  # Add attributes
  attr(A, "arOrder") <- length(phi)
  attr(A, "maOrder") <- length(theta)
  return(A)
}

#' Compute the conditional mean of an ARMA model
#'
#' @param h Numeric scalar, number of steps ahead.
#' @param X0 Numeric vector with length `p + q`, state vector of past values.
#' @param A Matrix with dimension `(p+q) x (p+q)`, companion matrix for ARMA model. See the function [ARMA_companion_matrix()].
#' @param b Numeric vector with length `p + q`. See the function [ARMA_vector_b()].
#' @param intercept Numeric scalar, intercept parameter.
#' @examples
#' A <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' intercept <- 0.2
#' X0 <- c(0.9, 0.2, -0.1, 0.3)
#' ARMA_expectation(h = 3, X0, A, b, intercept)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_expectation
#' @name ARMA_expectation
#' @export
#' @noRd
ARMA_expectation <- function(h = 10, X0, A, b, intercept = 0){
  # Dimension
  pq <- ncol(A)
  # Build intercept vector
  c_ <- c(intercept, rep(0, pq - 1))
  # Power of A
  A_pow_h <- purrr::map(1:h, ~pow_matrix(A, .x))
  # Identity matrix
  I <- diag(1, nrow = pq)
  # Compute expectation
  I_A_c_inv <- solve(I - A) %*% c_
  pow_sum <- purrr::map_dbl(A_pow_h, ~((I - .x) %*% I_A_c_inv  + .x %*% X0)[[1]])
  names(pow_sum) <- paste0("t+", 1:h)
  pow_sum
}

#' Compute the long-term variance of an ARMA model
#'
#' @inheritParams ARMA_expectation
#' @param sigma2 Numeric scalar, std. deviation of the residuals.
#' @examples
#' A <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' ARMA_variance(h = 3, A, b, sigma2 = 1)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_variance
#' @name ARMA_variance
#' @export
#' @noRd
ARMA_variance <- function(h = 1, A, b, sigma2 = 1){
  # Detect AR-order
  bb <- b %*% t(b)
  pow_sum <- list(a = diag(1, nrow(A), nrow(A)))
  A_pow_j <- diag(1, nrow(A))
  if (h > 1) {
    for(i in 1:(h-1)){
      A_pow_j <- A %*% A_pow_j
      pow_sum[[i+1]] <- pow_sum[[i]] + A_pow_j %*% bb %*% t(A_pow_j)
    }
  }
  # Forecasted variances
  pow_sum <- purrr::map_dbl(pow_sum, ~.x[1,1] * sigma2)
  names(pow_sum) <- paste0("t+", 1:h)
  pow_sum
}

#' Compute the ARMA conditional variance / covariance
#'
#' @inheritParams ARMA_variance
#' @param k Numeric scalar, number of steps ahead for the second lag.
#' @examples
#' A <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' ARMA_covariance(3, 1, A, b, sigma2 = 1)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_covariance
#' @name ARMA_covariance
#' @export
#' @noRd
ARMA_covariance <- function(h, k, A, b, sigma2 = 1){
  cv_t <- c(0)
  cv_x <- 0
  bb <- b %*% t(b)
  for(j in 0:(min(h, k)-1)){
    A_hj <- pow_matrix(A, h-1-j)
    A_kj <- pow_matrix(A, k-1-j)
    cv_x <- cv_x + A_hj %*% bb %*% t(A_kj)
    cv_t[j+1] <- cv_x[1,1]
  }
  sigma2 * cv_x
}

#' Next-step value of an ARMA process
#'
#' @inheritParams ARMA_expectation
#' @param eps Numeric, vector of new residuals.
#' @examples
#' # Companion matrix and vector b
#' A <- ARMA_companion_matrix(c(0.4), c())
#' b <- ARMA_vector_b(1,0)
#' # Initial value
#' X0 <- c(0.9)
#'
#' # Next step forecast
#' ARMA_next_step(1, X0, A, b, intercept = 0.2, eps = 0)
#' # 3-steps ahead forecast
#' ARMA_next_step(3,X0, A, b, intercept = 0.2, eps = 0)
#' # 10-step ahead forecast
#' ARMA_next_step(10, X0, A, b, intercept = 0.2, eps = 0)
#'
#' # Next step simulation
#' eps <- rnorm(1)
#' ARMA_next_step(1, X0, A, b, intercept = 0.2, eps = eps)
#' # 3-steps ahead simulation
#' eps <- rnorm(3)
#' ARMA_next_step(3, X0, A, b, intercept = 0.2, eps = eps)
#' # 10-step ahead simulation
#' eps <- rnorm(10)
#' ARMA_next_step(10, X0, A, b, intercept = 0.2, eps = eps)
#'
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_next_step
#' @name ARMA_next_step
#' @export
#' @noRd
ARMA_next_step <- function(h = 1, X0, A, b, intercept = 0, eps = 0){
  if (h > 1 & (length(eps) == 1 && eps == 0)) {
    eps <- rep(0, h)
  } else if (length(eps) != h){
    stop("The length of `eps` must be equal to `h` when specified!")
  } else if (length(X0) != ncol(A)){
    stop("The length of `X0` must be equal to `ncol(A)` when specified!")
  }
  # Initialize state vector
  x_t <- X0
  # Forecasting loop
  for(step in 1:h){
    # Next step state space vector
    x_t <-  A %*% x_t + b * eps[step]
    x_t[1] <- x_t[1] + intercept
  }
  return(x_t)
}

#' @description
#' Filter the time-series and compute fitted values and residuals.
#' @param x Numeric vector, time series to filter.
#' @inheritParams ARMA_expectation
#' @rdname ARMA_filter
#' @name ARMA_filter
#' @export
#' @noRd
ARMA_filter <- function(x, A, b, intercept = 0) {
  # AR order
  p <- attr(A, "arOrder")
  # MA order
  q <- attr(A, "maOrder")
  # Companion matrix
  A <- as.matrix(A); storage.mode(A) <- "double"
  # Fitted values
  x_hat <- .Call("ARMA_filter_c",
                 A,
                 as.numeric(b),
                 as.numeric(x),
                 as.integer(p),
                 as.integer(q),
                 as.numeric(intercept))
  return(x_hat)
}

#' Fast ARMA state-space h-step forecast and weights (C/BLAS)
#'
#' @param h Integer, steps ahead.
#' @param X0 Numeric vector of length p+q (state).
#' @param A  Numeric (p+q) x (p+q) companion matrix.
#' @param b  Numeric vector length p+q (shocks selector).
#' @param intercept Scalar intercept (0 if none).
#' @examples
#' h <- 1000
#' X0 <- c(0.2, 0.1)
#' A <- ARMA_companion_matrix(0.2, 0.1)
#' b <- ARMA_vector_b(1,1)
#' ARMA_forecast(h, X0, A, b, intercept = 0)
#'
#' @export
ARMA_forecast <- function(h, X0, A, b, intercept = 0) {
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  df_tT <- .Call("ARMA_forecast_c", A, as.numeric(X0), as.numeric(b), as.integer(h), as.numeric(intercept))
  df_tT <- dplyr::bind_rows(df_tT)
  # Add step
  df_tT$step <- 1:h
  df_T <- tail(df_tT, 1)[, c(5, 1)]
  df_T$weights <- list(df_tT)
  df_T
}

#' Compute the long-term variance of an AR model
#'
#' @param phi Numeric vector, AR parameters.
#' @param sigma2 Numeric scalar, std. deviation of the residuals.
#' @examples
#' AR_variance(c(0.5, 0.2))
#' AR_variance(c(0.2, 0.3, 0.2), 2)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @export
#' @noRd
AR_variance <- function(phi, sigma2 = 1){
  # Detect AR-order
  arOrder <- length(phi)
  if (arOrder == 1) {
    var <- sigma2/(1-phi[1]^2)
  } else if (arOrder == 2) {
    var <- sigma2*(1-phi[2])/((1 - phi[2])*(1 - phi[1]^2 - phi[2]^2) - 2*phi[1]^2*phi[2])
  } else if (arOrder == 3) {
    phi_tilde_1 <- (phi[1] + phi[2]*phi[3])/(1 - phi[2] - phi[3]^2 - phi[1]*phi[3])
    phi_tilde_2 <- (phi[1] + phi[3])*phi_tilde_1 + phi[2]
    phi_tilde_3 <- (phi[1]*phi_tilde_2 + phi[2]*phi_tilde_1 + phi[3])
    phi_tilde_0 <- 1/(1 - phi[1]*phi_tilde_1 - phi[2]*phi_tilde_2 - phi[3]*phi_tilde_3)
    var <- phi_tilde_0*sigma2
  } else if (arOrder == 4) {
    phi_1 <- (phi[1] + phi[3])/(1 - phi[4])
    phi_0 <- phi[2]/(1 - phi[4])
    psi_1 <- (phi_1*phi[3] + phi[1]*phi_0*phi[4]+ phi[1] + phi[3]*phi[4])
    psi_1 <- psi_1/(1-phi_1*(phi[3] + phi[1]*phi[4]) - phi[2]*(1 + phi[4]) - phi[4]^2)
    psi_2 <- phi_1*psi_1 + phi_0
    psi_3 <- phi[1]*psi_2 + phi[2]*psi_1 + phi[4]*psi_1 + phi[3]
    psi_4 <- phi[1]*psi_3 + phi[2]*psi_2 + phi[3]*psi_1 + phi[4]
    var <- sigma2/(1 - phi[1]*psi_1 - phi[2]*psi_2 - phi[3]*psi_3 - phi[4]*psi_4)
  } else {
    # Companion matrix
    phi <- matrix(phi, nrow = 1)
    A <- rbind(phi, diag(1, arOrder)[-c(arOrder),])
    b <- c(1, rep(0, arOrder-1))
    var <- ARMA_variance(A, b, sigma2)
  }
  return(var)
}

#' From constraint to uncontraint parameters
#' @note Version 1.0.0.
#' @noRd
AR_params_to_zeta <- function(phi){
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

#' From constraint to uncontraint parameters
#' @examples
#' phi <- c(0.3, 0.4, 0.1)
#' theta <- c(0.1, 0.8)
#' zeta <- ARMA_params_to_zeta(phi, theta)
#'
#' @note Version 1.0.0.
#' @export
#' @noRd
ARMA_params_to_zeta <- function(phi, theta){
  zeta_phi <- c()
  # Ensure not missing
  if (!missing(phi) && phi[1] != 0){
    zeta_phi <- AR_params_to_zeta(phi)
    names(zeta_phi) <- paste0("zeta_phi_", 1:length(phi))
  }

  zeta_theta <- c()
  # Ensure not missing
  if (!missing(theta) && theta[1] != 0){
    zeta_theta <- AR_params_to_zeta(theta)
    names(zeta_theta) <- paste0("zeta_theta_", 1:length(theta))
  }

  structure(
    list(
      zeta_phi = zeta_phi,
      zeta_theta = zeta_theta
    )
  )
}

#' From unconstraint to contraint parameters
#'
#' @note Version 1.0.0.
#' @noRd
AR_params_to_phi <- function(zeta){
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

#' From constraint to uncontraint parameters
#' @examples
#' zeta_phi <- c(0.7043836, 0.4652377, 0.1003353)
#' zeta_theta <- c(0.5493061, 1.0986123)
#' ARMA_params_to_phi(zeta_phi, zeta_theta)
#'
#' @note Version 1.0.0.
#' @export
#' @noRd
ARMA_params_to_phi <- function(zeta_phi, zeta_theta){
  # AR order
  p <- 0
  phi <- c()
  if (!missing(zeta_phi) && !is.null(zeta_phi)){
    p <- length(zeta_phi)
    J_phi <- AR_params_to_phi(zeta_phi)
    phi <- J_phi$phi
    names(phi) <- paste0("phi_", 1:p)
  }
  # MA order
  q <- 0
  theta <- c()
  if (!missing(zeta_theta) && !is.null(zeta_theta)){
    q <- length(zeta_theta)
    J_theta <- AR_params_to_phi(zeta_theta)
    theta <- J_theta$phi
    names(theta) <- paste0("theta_", 1:q)
  }
  # Jacobian matrix
  J <- matrix(0, p+q, p+q)
  if (p > 0){
    J[1:p, 1:p] <- J_phi$J_phi_zeta
  }
  if (q > 0){
    J[(p+1):(p+q), (p+1):(p+q)] <- J_theta$J_phi_zeta
  }

  structure(
    list(
      phi = phi,
      theta = theta,
      J = J
    )
  )
}

#' Log-likelihood with contraint parameters
ARMA_loss_logLik <- function(params, p = 0, q = 0, y, per_obs = FALSE){
  # Intercept
  intercept <- c(0)
  if (!is.na(params["intercept"])){
    intercept <- params["intercept"]
    params[which(names(params) != "intercept")]
  }
  # Extract AR parameters
  zeta_phi <- NULL
  if (p != 0){
    zeta_phi <- params[1:p]
  }
  # Extract MA parameters
  zeta_theta <- NULL
  if (q != 0){
    zeta_theta <- params[(p+1):(p+q)]
  }
  # Converted parameters
  phi_theta <- ARMA_params_to_phi(zeta_phi, zeta_theta)
  # Companion matrix and vector b
  A <- ARMA_companion_matrix(phi_theta$phi, phi_theta$theta)
  b <- ARMA_vector_b(p, q)
  # Fitted values
  y_hat <- ARMA_filter(y, A, b, intercept = intercept)
  # Fitted sigma
  idx_excluded <- 1:max(c(p, q))
  # Errors
  eps_hat <- (y - y_hat)[-idx_excluded]
  # Fitted std. deviation
  sigma_hat <- sqrt(mean(eps_hat^2))
  # Standardized residuals
  z_hat <- eps_hat / sigma_hat
  # Log-likelihoods
  loglik <- log(dnorm(z_hat) / sigma_hat)
  # Loss
  if (per_obs){
    loglik
  } else {
    sum(loglik)
  }
}

#' Sum of the errors with contraint parameters
ARMA_loss_CSS <- function(params, p = 0, q = 0, y){
  # Intercept
  intercept <- c(0)
  if (!is.na(params["intercept"])){
    intercept <- params["intercept"]
    params[which(names(params) != "intercept")]
  }
  # Extract AR parameters
  zeta_phi <- NULL
  if (p != 0){
    zeta_phi <- params[1:p]
  }
  # Extract MA parameters
  zeta_theta <- NULL
  if (q != 0){
    zeta_theta <- params[(p+1):(p+q)]
  }
  # Converted parameters
  phi_theta <- ARMA_params_to_phi(zeta_phi, zeta_theta)
  # Companion matrix and vector b
  A <- ARMA_companion_matrix(phi_theta$phi, phi_theta$theta)
  b <- ARMA_vector_b(p, q)
  # Fitted values
  y_hat <- ARMA_filter(y, A, b, intercept = intercept)
  # Fitted sigma
  idx_excluded <- 1:max(c(p, q))
  # Errors
  eps_hat <- (y - y_hat)[-idx_excluded]
  # SSE
  sum(eps_hat^2)
}

#' Fit with contraint parameters
ARMA_fit <- function(y, arOrder = 1, maOrder = 0, method = c("CSS", "ML")){
  # Fit method
  method <- match.arg(method, choices = c("CSS", "ML"))
  # Number of observations
  n <- length(y)
  # Total order
  pq <- arOrder + maOrder

  # Initialize unconstraint parameters
  zeta_phi <- runif(arOrder, -0.1, 0.1)
  if (arOrder > 0){
    names(zeta_phi) <- paste0("zeta_phi_", 1:arOrder)
  }
  zeta_theta <- runif(maOrder, -0.1, 0.1)
  if (maOrder > 0){
    names(zeta_theta) <- paste0("zeta_theta_", 1:maOrder)
  }
  # Initial parameters
  params <- c(zeta_phi, zeta_theta)
  # Loss function
  if (method == "CSS"){
    loss <- function(params) ARMA_loss_CSS(params, p = arOrder, q = maOrder, y = y)
  } else if (method == "ML"){
    loss <- function(params) -ARMA_loss_logLik(params, p = arOrder, q = maOrder, y = y, per_obs = FALSE)
  }
  # Optimization
  opt <- optim(params, loss)
  # Convert parameters
  opt_par <- opt$par
  # Initialize unconstraint parameters
  zeta_phi <- c()
  if (arOrder > 0){
    zeta_phi <- opt_par[1:arOrder]
    names(zeta_phi) <- paste0("zeta_phi_", 1:arOrder)
  }
  zeta_theta <- c()
  if (maOrder > 0){
    zeta_theta <- opt_par[(arOrder+1):pq]
    names(zeta_theta) <- paste0("zeta_theta_", 1:maOrder)
  }
  # QMLE parameters
  theta_star_qmle <- c(zeta_phi, zeta_theta)
  # Constraint parameters
  res <- ARMA_params_to_phi(zeta_phi, zeta_theta)
  theta_qmle <- c(res$phi, res$theta)
  # Numerical Hessian matrix at QMLE
  H <- numDeriv::hessian(func = ARMA_loss_logLik, x = theta_star_qmle,
                         p = arOrder, q = maOrder, y = y, per_obs = FALSE)
  # Numerical Score at QMLE
  S <- numDeriv::jacobian(func = ARMA_loss_logLik, x = theta_star_qmle,
                         p = arOrder, q = maOrder, y = y, per_obs = TRUE)
  # Cross products of the score
  B <- matrix(0, pq, pq)
  for(i in 1:nrow(S)) {
    B <- B + S[i,] %*% t(S[i,])
  }
  # Var-cov matrix (unbounded)
  V_star <- solve(-H)
  std.errors_zeta <- sqrt(diag(V_star))
  names(std.errors_zeta) <- names(theta_star_qmle)
  # Var-cov matrix (constraint)
  V_orig <- res$J %*% V_star %*% t(res$J)
  std.errors_phi <- sqrt(diag(V_orig))
  names(std.errors_phi) <- names(theta_qmle)

  # Sandwitch var-cov matrix (unbounded)
  V_rob_star <- V_star %*% B %*% V_star
  std.errors_rob_zeta <- sqrt(diag(V_rob_star))
  names(std.errors_rob_zeta) <- names(theta_star_qmle)
  # Sandwitch var-cov matrix (constraint)
  V_rob_orig <- res$J %*% V_rob_star %*% t(res$J)
  std.errors_rob_phi <- sqrt(diag(V_rob_orig))
  names(std.errors_rob_phi) <- names(theta_qmle)
  # ************************************************
  # Store loss
  res$loss <- opt$value
  names(res$loss) <- method
  # Store hessian
  res$H <- H
  # Store Vcov matrices
  res$vcov <- list(V_star = V_star, V_orig = V_orig,V_rob_star = V_rob_star, V_rob_orig = V_rob_orig)
  res$vcov_rob <- list(V_star = V_rob_star, V_orig = V_rob_orig)
  # Std.errors
  res$std.errors <- list(zeta = std.errors_zeta, phi = std.errors_phi)
  res$std.errors_rob <- list(zeta = std.errors_rob_zeta, phi = std.errors_rob_phi)
  res
}

