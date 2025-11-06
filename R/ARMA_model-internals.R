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

#' Construct the vector b of an ARMA model
#'
#' @param arOrder Numeric scalar, order of AR model.
#' @param maOrder Numeric scalar, order of MA model.
#' @examples
#' ARMA_vector_b(2,2)
#' ARMA_vector_b(2,0)
#' ARMA_vector_b(0,2)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_vector_b
#' @name ARMA_vector_b
#' @export
#' @noRd
ARMA_vector_b <- function(arOrder, maOrder){
  # Companion vector AR
  if (arOrder > 1){
    e_p <- c(1, rep(0, arOrder-1))
  } else {
    e_p <- 1
  }
  # Companion vector MA
  e_q <- c()
  if (maOrder == 1) {
    e_q <- 1
  } else if (maOrder > 1) {
    e_q <- c(1, rep(0, maOrder-1))
  }
  # Combine the two parts
  b <- matrix(c(e_p, e_q), ncol = 1)
  return(b)
}

#' Construct the companion matrix of an ARMA model
#'
#' @param phi Numeric vector, vector with AR parameters.
#' @param theta Numeric vector, vector with MA parameters.
#' @examples
#' ARMA_companion_matrix(c(0.4, 0.2), c(0.3, 0.1))
#' ARMA_companion_matrix(c(0.4, 0.2, 0.2), c(0.3, 0.1, 0.05))
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_companion_matrix
#' @name ARMA_companion_matrix
#' @export
#' @noRd
ARMA_companion_matrix <- function(phi, theta){
  if (missing(theta)){
    theta <- NULL
  }
  # AR order
  p <- length(phi)
  # MA order
  q <- length(theta)

  # Identity matrix
  I_p <- diag(1, nrow = p-1)
  # Add colum of zeros
  L_p <- cbind(I_p, rep(0, p - 1))
  # Matrix of zeros
  L_p <- cbind(L_p, matrix(0, nrow = p-1, ncol = q))

  if (q > 0) {
    # Identity matrix
    I_q <- diag(1, nrow = q-1)
    # Add colum of zeros
    L_q <- cbind(I_q, rep(0, q-1))
    # Matrix of zeros
    L_q <- cbind(matrix(0, nrow = q-1, ncol = p), L_q)
    A <- rbind(c(phi, theta), L_p, rep(0, p+q), L_q)
  } else {
    A <- rbind(c(phi), L_p)
  }
  return(A)
}

#' Compute the conditional mean of an ARMA model
#'
#' @param X0 Numeric vector, state vector of past values.
#' @param h Numeric scalar, number of steps ahead.
#' @param Phi Matrix, companion matrix for ARMA parameters. See the function [ARMA_companion_matrix()].
#' @param b Numeric vector, vector of 1 and 0 specific of the model. See the function [ARMA_vector_b()]
#' @param c Numeric vector, vector with the intercept at first position and zero otherwise.
#' @examples
#' Phi <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' c <- c(0.2, rep(0, 3))
#' X0 <- c(0.9, 0.2, -0.1, 0.3)
#' ARMA_expectation(X0, h = 3, Phi, b, c)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_expectation
#' @name ARMA_expectation
#' @export
#' @noRd
ARMA_expectation <- function(X0, h = 10, Phi, b, c){
  A_h <- pow_matrix(Phi, h, eigen_dec = FALSE)
  I_p <- diag(1, nrow = nrow(A_h))
  (I_p - A_h) %*% solve(I_p - Phi) %*% c  + A_h %*% X0
}

#' Compute the long-term variance of an ARMA model
#'
#' @inheritParams ARMA_expectation
#' @param sigma2 Numeric scalar, std. deviation of the residuals.
#' @examples
#' Phi <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' ARMA_variance(Phi, b, sigma2 = 1, h = 3)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_variance
#' @name ARMA_variance
#' @export
#' @noRd
ARMA_variance <- function(Phi, b, sigma2 = 1, h = 100){
  # Detect AR-order
  bb <- b %*% t(b)
  pow_sum <- sigma2
  A <- Phi
  for(i in 1:h){
    A_j <- pow_matrix(A, i, eigen_dec = FALSE)
    pow_sum <- pow_sum + A_j %*% bb %*% t(A_j)
  }
  pow_sum[1,1] * sigma2
}

#' Compute the ARMA conditional variance / covariance
#'
#' @inheritParams ARMA_variance
#' @param k Numeric scalar, number of steps ahead for the second lag.
#' @examples
#' Phi <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' ARMA_covariance(3, 1, Phi, b, sigma2 = 1)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_covariance
#' @name ARMA_covariance
#' @export
#' @noRd
ARMA_covariance <- function(h, k, Phi, b, sigma2 = 1){
  cv_t <- c(0)
  cv_x <- 0
  for(j in 0:(min(h, k)-1)){
    A_hj <- pow_matrix(Phi, h-1-j, eigen_dec = FALSE)
    A_kj <- pow_matrix(Phi, k-1-j, eigen_dec = FALSE)
    cv_x <- cv_x + A_hj %*% b %*% t(b) %*% t(A_kj)
    cv_t[j+1] <- cv_x[1,1]
  }
  sigma2 * cv_x
}

#' Compute the ARMA next-step value
#'
#' @inheritParams ARMA_expectation
#' @examples
#' Phi <- ARMA_companion_matrix(c(0.4, 0.1), c(0.1, 0.05))
#' b <- ARMA_vector_b(2,2)
#' X0 <- c(0.9, 0.2, -0.1, 0.3)
#' ARMA_next_step(X0, Phi, b, h = 3, intercept = 0.2, eps = 0)
#' @keywords ARMA
#' @note Version 1.0.0.
#' @rdname ARMA_next_step
#' @name ARMA_next_step
#' @export
#' @noRd
ARMA_next_step <- function(X0, Phi, b, h = 1, intercept = 0, eps = 0){
  if (h > 1 & (length(eps) == 1 && eps == 0)) {
    eps <- rep(0, h)
  } else if (length(eps) != h){
    stop("The length of `eps` must be equal to `h` when specified!")
  } else if (length(X0) != ncol(Phi)){
    stop("The length of `X0` must be equal to `ncol(Phi)` when specified!")
  }
  # Initialize state vector
  x_t <- X0
  # Forecasting loop
  for(step in 1:h){
    # Next step state space vector
    x_t <-  Phi %*% x_t + b * eps[step]
    x_t[1] <- x_t[1] + intercept
  }
  return(x_t)
}

#' Compute the robust std. errors for an ARMA model
#'
#' @keywords ARMA
#' @note Version 1.0.0.
#' @export
#' @noRd
ARMA_HAC_error <- function(Yt, eps, ar_params = NULL, ma_params = NULL) {
  stopifnot(length(Yt) == length(eps))
  n <- length(Yt)

  p <- length(ar_params)
  q <- length(ma_params)
  max_lag <- max(p, q)

  # Build regression frame
  df <- data.frame(y = Yt[(max_lag + 1):n])

  # Add AR terms
  if (p > 0) {
    for (i in 1:p) {
      df[[paste0("ar", i)]] <- Yt[(max_lag + 1 - i):(n - i)]
    }
  }

  # Add MA terms
  if (q > 0) {
    for (j in 1:q) {
      df[[paste0("ma", j)]] <- eps[(max_lag + 1 - j):(n - j)]
    }
  }

  # Create formula with all regressors
  rhs <- c(
    if (p > 0) paste0("ar", 1:p),
    if (q > 0) paste0("ma", 1:q)
  )
  formula <- as.formula(paste("y ~ -1 + ", paste(rhs, collapse = " + ")))

  # Fit lm with NO optimization (coefficients will be ignored)
  fit <- lm(formula, data = df)
  fit$residuals <- eps[(max_lag + 1):n]
  # Replace coefficients with fixed ones
  fixed_coef <- c(`(Intercept)` = 0)
  if (p > 0) fixed_coef[paste0("ar", 1:p)] <- ar_params
  if (q > 0) fixed_coef[paste0("ma", 1:q)] <- ma_params
  fit$coefficients <- fixed_coef
  # Standard errors
  std.errors <- sqrt(diag(sandwich::NeweyWest(fit)))
  names(std.errors) <- names(fixed_coef)[-1]
  std.errors
}

#' Compute the Jacobian for ARMA coefficients
#'
#' @keywords ARMA
#' @note Version 1.0.0.
#' @export
#' @noRd
ARMA_jacobian_params <- function(Yt, eps, ARMA){

  # Length Yt
  n <- length(Yt)
  # Order AR
  arOrder <- ARMA$order[1]
  # Order MA
  maOrder <- ARMA$order[2]
  # Starting lag
  i_start <- max(c(arOrder, maOrder)) + 1
  # Extract ARMA parameter
  theta <- ARMA$theta
  # Initialization
  d_mu_d_phi <- matrix(0, n, arOrder)
  d_mu_d_theta <- matrix(0, n, maOrder)
  for(t in i_start:n){
    # Derivative AR parameters
    if (arOrder != 0){
      for(i in 1:arOrder){
        d_mu_d_phi[t, i] <- Yt[t-i]
        if (maOrder != 0){
          for(j in 1:maOrder) {
            d_mu_d_phi[t, i] <- d_mu_d_phi[t, i] - theta[j] * d_mu_d_phi[t-j, i]
          }
        }
      }
    }
    # Derivative MA parameters
    if (maOrder != 0) {
      for(i in 1:maOrder){
        d_mu_d_theta[t, i] <- eps[t-i]
        for(j in 1:maOrder) {
          d_mu_d_theta[t, i] <- d_mu_d_theta[t, i] - theta[j] * d_mu_d_theta[t-j, i]
        }
      }
    }
  }

  J_phi_theta <- cbind(d_mu_d_phi, d_mu_d_theta)
  colnames(J_phi_theta) <- names(ARMA$coefficients[-1])
  return(J_phi_theta)
}



