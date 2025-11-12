#' Standard GARCH expected value formula
#'
#' @param h Numeric scalar, number of steps ahead.
#' @param omega Numeric scalar, intercept of the GARCH model.
#' @param alpha Numeric scalar, ARCH parameter of the GARCH model.
#' @param beta Numeric scalar, GARCH parameter of the GARCH model.
#' @param e_x2 Numeric vector, second moment of the residuals with length equal to `h`.
#' @param sigma2_t Numeric scalar, GARCH variance at time t.
#'
#' @examples
#' # Forecast horizon
#' h <- 10
#' # GARCH parameters
#' alpha <- 0.08
#' beta <- 0.35
#' omega <- 1.2*(1 - alpha - beta)
#' # Moments
#' e_x2 = 1
#' e_x4 = 3
#' # Initial values for variance
#' sigma2_t <- 1.1
#' sigma4_t <- sigma2_t^2
#'
#' e_sigma2_h(10, omega, alpha, beta, e_x2, sigma2_t)
#' e_sigma2_h_mix(10, omega, alpha, beta, e_x2, sigma2_t)
#' e_sigma4_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' v_sigma2_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' v_sigma_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' e_sigma12_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' e_sigma32_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' @name e_sigma2_h
#' @rdname e_sigma2_h
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
e_sigma2_h <- function(h, omega, alpha, beta, e_x2 = 1, sigma2_t){
  # Derived quantities
  lambda <- alpha * e_x2[1] + beta
  # Long term expectation
  sigma2_inf <- omega / (1 - lambda)
  # Forecasted second moment
  sigma2_h <- sigma2_inf + lambda^(1:h) * (sigma2_t - sigma2_inf)
  sigma2_h <- c(sigma2_t, sigma2_h)
  names(sigma2_h) <- paste0("t+", 0:h)
  sigma2_h
}

#' Iterative GARCH expected value formula
#'
#' @inheritParams e_sigma2_h
#' @name e_sigma2_h_mix
#' @rdname e_sigma2_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
e_sigma2_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, sigma2_t){
  # Initialize a vector
  sigma2_h <- c(`t+0` = sigma2_t)
  # Check forecast horizon
  if (h == 0) return(sigma2_h)
  # Check second moment
  m2 <- e_x2
  if (length(e_x2) == 1 & h > 1){
    m2 <- rep(e_x2, h)
  }
  # Derived quantities
  lambda <- alpha * m2 + beta
  # Add next step variance
  sigma2_h <- c(sigma2_h, `t+1` = omega + lambda[1] * sigma2_t)
  # Check forecast horizon
  if (h == 1) return(sigma2_h)
  # Iterate forecast
  for(h_ahead in 2:h){
    sigma2_h <- c(sigma2_h, omega + lambda[h_ahead] * sigma2_h[h_ahead])
  }
  # Assign standard names
  names(sigma2_h) <- paste0("t+", 0:(length(sigma2_h)-1))
  return(sigma2_h)
}


#' Iterative GARCH second moment formula
#'
#' @inheritParams e_sigma2_h
#' @param e_x4 Numeric vector, fourth moment of the residuals with length equal to `h`.
#' @param sigma4_t Numeric scalar, GARCH squared variance at time t.
#'
#' @name e_sigma4_h_mix
#' @rdname e_sigma4_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
e_sigma4_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  if (h == 0){
    return(c(`t+0` = sigma4_t))
  }
  sigma4_iter <- numeric(h + 1)
  sigma4_iter[1] <- sigma4_t
  # Iterative formula for E[sigma^2]
  sigma2_iter <- numeric(h + 1)
  sigma2_iter[1] <- sqrt(sigma4_t)
  for(i in 1:h){
    lambda <- alpha * e_x2[i] + beta
    gamma <- alpha^2 * e_x4[i] + beta * (2 * alpha * e_x2[i] + beta)
    sigma2_iter[i+1] <- omega + lambda * sigma2_iter[i]
    b_t <- omega^2 + 2 * omega * lambda * sigma2_iter[i]
    sigma4_iter[i+1] <- b_t + gamma * sigma4_iter[i]
  }
  sigma4_iter
}

#' Iterative GARCH variance formula
#'
#' @inheritParams e_sigma4_h_mix
#'
#' @name v_sigma2_h_mix
#' @rdname v_sigma2_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
v_sigma2_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Second moment GARCH variance
  e_sigma4 <- e_sigma4_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance
  e_sigma4 - e_sigma2^2
}

#' Iterative GARCH variance formula (approximated)
#'
#' @inheritParams e_sigma4_h_mix
#'
#' @name v_sigma_h_mix
#' @rdname v_sigma_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
v_sigma_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Expectation GARCH std. dev
  e_sigma <- e_sigma12_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Variance
  e_sigma2 - e_sigma^2
}

#' Conditional first moment GARCH std. dev (approximated)
#'
#' @inheritParams e_sigma4_h_mix
#'
#' @rdname e_sigma12_h_mix
#' @name e_sigma12_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
e_sigma12_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance GARCH variance
  v_sigma2 <- v_sigma2_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Moment to power 1/2 (Approximated)
  e_sigma2^(1/2) - (1/8) * v_sigma2 / sqrt(e_sigma2)^3
}

#' Conditional third moment GARCH std. dev (approximated)
#'
#' @inheritParams e_sigma4_h_mix
#'
#' @name e_sigma32_h_mix
#' @rdname e_sigma32_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
e_sigma32_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance variance
  v_sigma2 <- v_sigma2_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Moment to power 3/2 (Approximated)
  e_sigma2^(3/2) + (3/8) * v_sigma2 / sqrt(e_sigma2)
}

#' Covariance between two GARCH variances
#'
#' @inheritParams e_sigma4_h_mix
#'
#' @name cov_sigma2_h_mix
#' @rdname cov_sigma2_h_mix
#' @keywords GARCH
#' @note Version 1.0.0.
#' @export
#' @noRd
cov_sigma2_h_mix <- function(h, s, omega, alpha, beta, e_x2, e_x4, sigma4_t){
  # Compute the variances
  v_sigma2_s <- v_sigma2_h_mix(h-1, omega, alpha, beta, e_x2, e_x4, sigma4_t)[-c(1:(s))]
  # Extract the moments from the indexes
  m2 <- e_x2[s:(h-1)]
  # Compute lambda
  lambda <- (alpha * m2 + beta)
  # Compute cumulative product
  prod_lambda <- cumprod(lambda)[(h-s):1]
  # Compute covariances
  cv_sigma2_hs <- prod_lambda * v_sigma2_s
  # Assign correct names
  names(cv_sigma2_hs) <- paste0(names(v_sigma2_s),"|",paste0("t+",h))
  cv_sigma2_hs
}
