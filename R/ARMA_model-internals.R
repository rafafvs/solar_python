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
  A_pow <- pow_matrix(A, h)
  # Identity matrix
  I <- diag(1, nrow = pq)
  # Compute expectation
  (I - A_pow) %*% solve(I - A) %*% c_  + A_pow %*% X0
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
  pow_sum <- sigma2
  A_pow_j <- diag(1, nrow(A))
  for(i in 1:h){
    A_pow_j <- A %*% A_pow_j
    pow_sum <- pow_sum + A_pow_j %*% bb %*% t(A_pow_j)
  }
  pow_sum[1,1] * sigma2
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
  .Call("ARMA_filter_c",
        A,
        as.numeric(b),
        as.numeric(x),
        as.integer(p),
        as.integer(q),
        as.numeric(intercept))
}

#' Compute the Jacobian for ARMA coefficients
#'
#' @keywords ARMA
#' @note Version 1.0.0.
#' @export
#' @noRd
ARMA_jacobian_params <- function(x, eps, ARMA){
  # Length of the time series
  n <- length(x)
  # Order AR
  arOrder <- ARMA$order[1]
  # Order MA
  maOrder <- ARMA$order[2]
  # Starting lag
  i_start <- max(c(arOrder, maOrder)) + 1
  # Extract MA parameters
  theta <- ARMA$theta
  # Initialization
  d_mu_d_phi <- matrix(0, n, arOrder)
  d_mu_d_theta <- matrix(0, n, maOrder)
  for(t in i_start:n){
    # Derivative AR parameters
    if (arOrder != 0){
      for(i in 1:arOrder){
        d_mu_d_phi[t, i] <- x[t-i]
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
