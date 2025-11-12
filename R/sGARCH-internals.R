#' GARCH_vector_b(0,2)
#' @keywords GARCH
#' @note Version 1.0.0.
#' @rdname GARCH_vector_b
#' @name GARCH_vector_b
#' @export
#' @noRd
GARCH_vector_b <- function(archOrder = 1, garchOrder = 0){
  # Arch order
  p <- archOrder
  # Basis vector
  e_p <- c(1)
  if (p > 0){
    e_p <- c(e_p, rep(0, p - 1))
  }
  # Garch order
  q <- garchOrder
  # Zero vector
  e_q <- c()
  if (q > 0){
    e_q <- rep(0, q)
  }
  matrix(c(e_q, e_p), ncol = 1)
}

#' Fast sGARCH filter (C)
#'
#' @param x numeric residuals/innovations.
#' @param omega Numeric scalar, strictly positive, intercept of the GARCH model.
#' @param alpha Numeric scalar, strictly positive, ARCH parameters.
#' @param beta  numeric length q (>=0).
#' @param eps0 optional numeric initial epsilons to prepend (length p+q).
#' @param sigma20 optional numeric initial variances to prepend (length p+q).
#' @return numeric vector of conditional variances sigma2, same length as the working series (prepended if initials are provided).
#' @examples
#' # Simple GARCH(1,1)
#' set.seed(1)
#' x <- rnorm(1000)
#' omega <- 0.1; alpha <- 0.05; beta <- 0.9
#' sGARCH_filter(x, omega, alpha, beta)  # no initials
#'
#' # With initials (k=max(p,q)=1)
#' sGARCH_filter(x, omega, alpha, beta, eps0 = 0, sigma20 = 0.2)
#'
#' @export
sGARCH_filter <- function(x, omega, alpha, beta, eps0 = NULL, sigma20 = NULL) {
  if (is.null(eps0) != is.null(sigma20)){
    stop("Either provide both eps0 and sigma20, or neither.")
  }
  .Call("sGARCH_filter_c",
        as.numeric(x),
        as.numeric(omega),
        as.numeric(alpha),
        as.numeric(beta),
        if (is.null(eps0)) NULL else as.numeric(eps0),
        if (is.null(sigma20)) NULL else as.numeric(sigma20))
}

#' Next-step value of an ARMA process
#'
#' @param eps Numeric, vector of new residuals.
#' @examples
#' # Companion matrix and vector b
#' omega <- 0.3
#' alpha <- c(0.2, 0.04)
#' beta <- c(0.1, 0.04, 0.02)
#' # Initial values
#' sigma20 <- c(`t-1` = 0.9, `t-2` = 0.2, `t-3` = 0.1)
#' eps0 <- c(`t-1` = 0.2, `t-2` = 0.4, `t-3` = 0.4)
#' # Next step forecast
#' sGARCH_next_step(omega, alpha, beta, sigma20, eps0)
#'
#' @keywords GARCG
#' @note Version 1.0.0.
#' @rdname sGARCH_next_step
#' @name sGARCH_next_step
#' @export
#' @noRd
sGARCH_next_step <- function(omega, alpha, beta, eps0 = 0, sigma20 = 1){
  eps0 <- eps0[1:length(alpha)]
  sigma20 <- sigma20[1:length(beta)]
  # Next step variance
  sigma2 <- omega + sum(alpha * eps0^2) + sum(beta * sigma20)
  return(c(sigma2 = sigma2[[1]]))
}


#' Simulate scenarios for an sGARCH(p, q) with 2-component Gaussian mixture shocks
#'
#' Generates \code{B} Monte Carlo paths of length \code{n} from an sGARCH process
#' with innovations \eqn{\varepsilon_t = \sqrt{\sigma_t^2}\,U_t}, where
#' \eqn{U_t} is a two-component Gaussian mixture with time-varying or constant
#' mixing probability.
#'
#' @param B Integer, number of scenarios (paths).
#' @param n Integer, path length.
#' @inheritParams sGARCH_filter
#' @param mu Numeric length-2, component means \code{c(mu1, mu2)}.
#' @param sd Numeric length-2, component standard deviations \code{c(sd1, sd2)}.
#' @param p_mix Numeric scalar or length-\code{n} vector in (0,1), mixture prob for component 1.
#'
#' @return A list of length \code{B}; each element is a tibble with columns \code{t, j, x, sigma2}.
#'
#' @examples
#' set.seed(123)
#' sim <- sGARCH_scenarios(B = 10, n = 500,
#' omega = 0.05, alpha = 0.05, beta = 0.9,
#' mu = c(0, 0), sd = c(1, 2), p_mix = 0.2)
#'
#' @export
sGARCH_scenarios <- function(B = 100, n = 1000, eps0 = 0, sigma2_0, omega, alpha, beta, mu, sd, probs){
  # ARCH order
  p <- length(alpha)
  # GARCH order
  q <- length(beta)
  # Maximum lag
  k <- max(c(p, q))
  # Simulate scenarios
  scenarios <- list()
  for(b in 1:B){
    Bt <- rbinom(n, 1, probs[1])
    # Standardized normal residuals
    z <- rnorm(n, 0, 1)
    # Mixture simulation
    u <- (mu[1]*Bt + mu[2]*(1-Bt)) + (sd[1]*Bt + sd[2]*(1-Bt)) * z
    # Initialization
    eps <- c(eps0)
    sigma2 <- c(sigma2_0)
    for(t in (k+1):n){
      sigma2[t] <- omega + sum(alpha * eps[(t-1):(t-1-p+1)]^2) + sum(beta * sigma2[(t-1):(t-1-q+1)])
      eps[t] <- sqrt(sigma2[t]) * u[t]
    }
    scenarios[[b]] <- dplyr::tibble(t = 1:n, j = b, x = eps, sigma2 = sigma2)
  }
  return(scenarios)
}

#' Map unconstrained parameters to constrained GARCH coefficients
#'
#' Transforms \code{zeta = (eta0, kappa_1..kappa_{p+q}, [xi|psi])} to
#' constrained \code{phi = (omega, alpha_1..alpha_p, beta_1..beta_q)} while
#' enforcing \eqn{\alpha_i \ge 0}, \eqn{\beta_j \ge 0}, and
#' \eqn{\sum \alpha_i + \sum \beta_j < 1}.
#'
#' @param zeta Numeric vector of unconstrained parameters: \code{eta0} (mass),
#'   \code{kappa_1..kappa_{p+q}} (softmax logits), and optionally \code{xi} or \code{psi}
#'   depending on \code{mode}.
#' @param p Integer, ARCH order.
#' @param q Integer, GARCH order.
#' @param mode Character, one of \code{"unitOmega"}, \code{"targetSigma2"}, \code{"freeOmega"}.
#' @param eps Numeric buffer in [0,1), enforces \eqn{\tau \le 1-\varepsilon}.
#'
#' @return Named numeric vector \code{c(omega, alpha_1..alpha_p, beta_1..beta_q)}.
#'
#' @examples
#' # eta0, kappa1..kappa2 for p+q=2
#' zeta <- c(0, rep(0, 2))
#' GARCH_params_to_phi(zeta, p = 1, q = 1)
#'
#' @export
sGARCH_params_to_phi <- function(zeta, p, q, mode = "unitOmega", eps = 0) {

  sigmoid <- function(x) 1/(1 + exp(-x))
  softmax <- function(kappa) {
    z <- kappa - max(kappa)
    ez <- exp(z)
    ez / sum(ez)
  }

  mode <- match.arg(mode, choices = c("unitOmega","targetSigma2","freeOmega"))
  m <- q + p
  stopifnot(length(zeta) >= 1 + m)

  eta0  <- zeta[1]
  kappa <- zeta[1 + seq_len(m)]
  s  <- sigmoid(eta0)                  # in (0,1)
  tau <- eps + (1 - eps) * s           # in (eps,1)
  w   <- softmax(kappa)                # length m, sums to 1

  alpha <- c(alpha1 = 0)
  if (p > 0) {
    alpha <- tau * w[seq_len(p)]
    names(alpha) <- paste0("alpha",  1:p)
  }
  beta <- c(beta1 = 0)
  if (q > 0) {
    beta  <- tau * w[p + seq_len(q)]
    names(beta) <- paste0("beta",  1:q)
  }

  if (mode == "unitOmega") {
    omega <- 1 - tau
  } else if (mode == "targetSigma2") {
    xi     <- zeta[1 + m + 1]
    sigma2 <- exp(xi)
    omega  <- sigma2 * (1 - tau)
  } else { # freeOmega
    psi   <- zeta[1 + m + 1]
    omega <- exp(psi)
  }

  c(omega = omega[[1]], alpha, beta)
}

#' Map constrained GARCH coefficients back to unconstrained parameters
#'
#' Inverts \code{phi = (omega, alpha_1..alpha_p, beta_1..beta_q)} into
#' \code{zeta = (eta0, kappa_1..kappa_{p+q}, [xi|psi])}, fixing the softmax gauge
#' by setting the last \code{kappa} to 0.
#'
#' @param phi Named numeric vector \code{c(omega, alpha_1..alpha_p, beta_1..beta_q)}.
#' @inheritParams sGARCH_params_to_phi
#'
#' @return Numeric vector \code{(eta0, kappa_1..kappa_{p+q}, [xi|psi])}.
#'
#' @examples
#' phi <- c(omega = 0.1, alpha1 = 0.05, beta1 = 0.9)
#' GARCH_params_to_zeta(phi, p = 1, q = 1)
#'
#' @export
sGARCH_params_to_zeta <- function(phi, p, q, mode = "unitOmega", eps = 0) {

  logit   <- function(p) log(p/(1-p))

  mode <- match.arg(mode, choices = c("unitOmega","targetSigma2","freeOmega"))
  m <- q + p
  stopifnot(length(phi) == 1 + m)

  omega <- phi[1]
  alpha <- phi[1 + seq_len(p)]
  beta  <- phi[1 + p + seq_len(q)]

  tau <- sum(alpha) + sum(beta)        # must be in (eps,1)
  if (!(tau > eps && tau < 1)) stop("Invalid tau: sum(alpha)+sum(beta) must be in (eps,1).")

  # Recover eta0 from tau = eps + (1-eps)*sigmoid(eta0)
  s    <- (tau - eps) / (1 - eps)
  if (!(s > 0 && s < 1)) stop("Invalid s from tau; check eps.")
  eta0 <- logit(s)

  # Weights and kappa (softmax inverse up to a constant -> fix last to 0)
  w <- c(alpha, beta) / tau
  if (any(w <= 0)) stop("All alpha/beta must be > 0 to invert softmax (or jitter).")
  kappa <- c(log(w[1:(m-1)] / w[m]), 0)

  if (mode == "unitOmega") {
    # Consistency check: omega should equal 1 - tau
    if (abs(omega - (1 - tau)) > 1e-10)
      stop("phi not consistent with unitOmega: omega != 1 - sum(alpha,beta).")
    zeta <- c(eta0, kappa)
    names(zeta) <- c("eta0", paste0("kappa", 1:m))
    return(zeta)
  }

  if (mode == "targetSigma2") {
    sigma2 <- omega / (1 - tau)        # target unconditional variance
    xi <- log(sigma2)
    zeta <- c(eta0, kappa, xi)
    names(zeta) <- c("eta0", paste0("kappa", 1:m), "xi")
    return(zeta)
  }

  # freeOmega
  psi  <- log(omega)
  zeta <- c(eta0, kappa, psi)
  names(zeta) <- c("eta0", paste0("kappa", 1:m), "psi")
  zeta
}

#' Jacobian of the (unconstrained -> constrained) GARCH reparametrization
#'
#' Builds the Jacobian matrix \eqn{\partial \phi / \partial \zeta} for the mapping
#' \code{zeta -> phi}, assuming the **reduced** softmax parameterization
#' (last \code{kappa} fixed to 0). Rows are ordered as
#' \code{(omega, alpha_1..alpha_p, beta_1..beta_q)}.
#'
#' @inheritParams sGARCH_params_to_phi
#'
#' @return A numeric matrix of size \code{(1+p+q) x (1+(p+q-1)+extra)},
#'   where \code{extra} is 0, 1 for \code{xi}, or 1 for \code{psi} depending on \code{mode}.
#'
#' @examples
#' zeta <- c(0, 0, 0)  # eta0, kappa1 (since last fixed), and psi/xi if needed
#' GARCH_params_to_zeta_jacobian(zeta, p = 1, q = 1)
#'
#' @export
sGARCH_params_to_zeta_jacobian <- function(zeta, p, q, mode = "unitOmega", eps = 0) {

  mode <- match.arg(mode, choices = c("unitOmega","targetSigma2","freeOmega"))
  m <- p + q

  eta0  <- zeta[1]
  kappa <- zeta[1 + seq_len(m - 1)]
  extra <- if (length(zeta) > 1 + (m - 1)) zeta[1 + (m - 1) + 1] else NULL

  s    <- 1/(1 + exp(-eta0))
  dtau <- (1 - eps) * s * (1 - s)
  tau  <- eps + (1 - eps) * s

  ez <- exp(kappa)
  denom <- 1 + sum(ez)
  w <- c(ez / denom, 1 / denom)  # length m

  if (m > 1) {
    colnames_vec <- c("eta0", paste0("kappa", 1:(m - 1)))
  } else {
    colnames_vec <- c("eta0")
  }

  if (mode == "targetSigma2") colnames_vec <- c(colnames_vec, "xi")
  if (mode == "freeOmega")    colnames_vec <- c(colnames_vec, "psi")

  J <- matrix(0, nrow = 1 + m, ncol = length(colnames_vec))
  # Row names
  row_names <- c("omega")
  if (p > 0){
    row_names <- c(row_names, paste0("alpha", 1:p))
  }
  if (q > 0){
    row_names <- c(row_names, paste0("beta", 1:q))
  }
  rownames(J) <- row_names
  colnames(J) <- colnames_vec

  dwdk <- function(i, r) {
    if (i < m) {
      wi <- w[i]; wr <- w[r]
      wi * ((i == r) - wr)
    } else {
      -w[m] * w[r]
    }
  }

  # alpha rows (1..p)
  if (p > 0) {
    for (i in 1:p) {
      J[1 + i, "eta0"] <- w[i] * dtau
      if (m > 1){
        for (r in 1:(m - 1)) J[1 + i, 1 + r] <- tau * dwdk(i, r)
      }
    }
  }

  # beta rows (p+1 .. p+q)
  if (q > 0) {
    for (j in 1:q) {
      idx <- p + j; row <- 1 + p + j
      J[row, "eta0"] <- w[idx] * dtau
      if (m > 1){
        for (r in 1:(m - 1)) J[row, 1 + r] <- tau * dwdk(idx, r)
      }
    }
  }

  if (mode == "unitOmega") {
    J["omega", "eta0"] <- -dtau
  } else if (mode == "targetSigma2") {
    sigma2 <- exp(extra); omega <- sigma2 * (1 - tau)
    J["omega", "eta0"] <- -sigma2 * dtau
    J["omega", "xi"]   <- omega
  } else {
    omega <- exp(extra)
    J["omega", "psi"]  <- omega
  }
  J
}

#' sGARCH Gaussian QMLE log-likelihood (parameterized by \eqn{\zeta})
#'
#' Computes the (per-observation or total) log-likelihood for an sGARCH(p, q)
#' under standard normal innovations, using the unconstrained parameter vector
#' \code{zeta} and the mapping in \code{GARCH_params_to_phi()}.
#'
#' @param zeta Unconstrained parameter vector (see \code{GARCH_params_to_phi()}).
#' @param y Numeric vector of observations or residuals.
#' @param archOrder Integer \eqn{p}.
#' @param garchOrder Integer \eqn{q}.
#' @param per_obs Logical; if \code{TRUE}, returns the log-likelihood by observation,
#'   else returns the scalar sum.
#' @inheritParams sGARCH_params_to_phi
#'
#' @return Numeric vector (if \code{per_obs=TRUE}) or numeric scalar (sum log-likelihood).
#'
#' @export
sGARCH_loglik <- function(y, weights, omega, alpha, beta, eps0 = NULL, sigma20 = NULL, per_obs = FALSE){
  if (missing(weights)){
    weights <- rep(1, length(y))
  }
  # GARCH filter
  sigma <- sqrt(sGARCH_filter(y, omega, alpha, beta, eps0, sigma20))
  # Standardized residuals
  z_hat <- y / sigma
  # Log-likelihood
  loglik <- log(dnorm(z_hat) / sigma) * weights
  if (!per_obs){
    loglik <- sum(loglik)
  }
  return(loglik)
}

#' sGARCH Gaussian QMLE
#'
#' Fits an sGARCH(p, q) by Gaussian QMLE using the unconstrained parameterization
#' \code{zeta} and the mapping in \code{GARCH_params_to_phi()}. Returns both
#' unconstrained and constrained estimates, and the Hessian-based covariance of
#' the unconstrained parameters. (If a Jacobian consistent with the chosen
#' parameterization is supplied, delta-method SEs for constrained params can be computed.)
#'
#' @param x Numeric vector to fit.
#' @param archOrder Integer \eqn{p}.
#' @param garchOrder Integer \eqn{q}.
#' @inheritParams sGARCH_params_to_phi
#'
#' @return A list with elements:
#' \describe{
#'   \item{par_unconstrained}{Estimated \code{zeta}.}
#'   \item{par_constrained}{Estimated \code{phi} (omega, alphas, betas).}
#'   \item{vcov_unconstrained}{Hessian-based covariance matrix of \code{zeta}.}
#'   \item{se_unconstrained}{Standard errors of \code{zeta}.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(500)
#' fit <- sGARCH_fit(x, archOrder = 1, garchOrder = 1)
#'
#' @export
sGARCH_fit <- function(y, weights, archOrder = 1, garchOrder = 1, mode = c("unitOmega","targetSigma2","freeOmega")){
  # Match argument
  mode <- match.arg(mode, choices = c("unitOmega","targetSigma2","freeOmega"))
  # Initial random parameters
  if (mode == "unitOmega") {
    # Total order of arch and GARCH parameters
    pq <- archOrder + garchOrder
    zeta <- runif(pq)
  } else {
    # Total order of arch and GARCH parameters
    pq <- archOrder + garchOrder + 1
    zeta <- runif(pq)
  }
  # ***************************************************************
  # Log likelihood function
  logLik_ <- function(x, per_obs = FALSE){
    # Convert bound parameters into unbounded parameters
    if (mode == "unitOmega") {
      zeta <- c(x, 0)
    } else {
      zeta <- c(x[-length(x)], 0, x[length(x)])
    }
    # Unconstraint parameters
    phi <-  sGARCH_params_to_phi(zeta, archOrder, garchOrder, mode)
    # Parameters names
    params_names <- names(phi)
    # GARCH filter
    omega <- phi[stringr::str_detect(params_names, "omega")]
    alpha <- phi[stringr::str_detect(params_names, "alpha")]
    beta <- phi[stringr::str_detect(params_names, "beta")]
    sGARCH_loglik(y, weights, omega, alpha, beta, eps0 = NULL, sigma20 = NULL, per_obs = per_obs)
  }

  # Optimal parameters
  if (pq == 1){
    opt <- optim(zeta, function(params) -logLik_(params, per_obs = FALSE),
                 lower = -10, upper = 10, method = "Brent")
  } else {
    opt <- optim(zeta, function(params) -logLik_(params, per_obs = FALSE))
  }

  # QMLE params (unconstraint)
  theta_star_qmle <- opt$par

  # Convert bound parameters into unbounded parameters
  if (mode == "unitOmega") {
    theta_star_qmle_full <- c(theta_star_qmle, 0)
  } else {
    theta_star_qmle_full <- c(theta_star_qmle[-length(theta_star_qmle)], 0, theta_star_qmle[length(theta_star_qmle)])
  }

  # QMLE params (constraint)
  theta_qmle <- sGARCH_params_to_phi(theta_star_qmle_full, archOrder, garchOrder, mode)
  # Jacobian
  J <- sGARCH_params_to_zeta_jacobian(theta_star_qmle, archOrder, garchOrder, mode)
  # Hessian and robust std. errors
  sandwitch <- logLik_std.errors(theta_star_qmle, logLik_, robust = TRUE)

  # Sandwitch var-cov matrix (constraint)
  V_orig <- J %*% sandwitch$V %*% t(J)
  std.errors <- sqrt(diag(V_orig))
  names(std.errors) <- names(theta_qmle)[theta_qmle != 0]

  # Output
  structure(
    list(
      coef_star = theta_star_qmle,
      std.errors_star = sandwitch$std.errors,
      coef = theta_qmle,
      std.errors = std.errors,
      log.likelihoods = logLik_(theta_star_qmle, TRUE)
    )
  )
}

#' sGARCH Gaussian QMLE with rugarch or manual routine
#'
#' Fits an sGARCH(p, q) by Gaussian QMLE using the unconstrained parameterization
#' \code{zeta} and the mapping in \code{GARCH_params_to_phi()}. Returns both
#' unconstrained and constrained estimates, and the Hessian-based covariance of
#' the unconstrained parameters. (If a Jacobian consistent with the chosen
#' parameterization is supplied, delta-method SEs for constrained params can be computed.)
#'
#' @param x Numeric vector to fit.
#' @param archOrder Integer \eqn{p}.
#' @param garchOrder Integer \eqn{q}.
#' @inheritParams sGARCH_params_to_phi
#'
#' @return A list with elements:
#' \describe{
#'   \item{par_unconstrained}{Estimated \code{zeta}.}
#'   \item{par_constrained}{Estimated \code{phi} (omega, alphas, betas).}
#'   \item{vcov_unconstrained}{Hessian-based covariance matrix of \code{zeta}.}
#'   \item{se_unconstrained}{Standard errors of \code{zeta}.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(500)
#' fit <- sGARCH_fit_rugarch(x, archOrder = 1, garchOrder = 1)
#'
#' @export
sGARCH_fit_rugarch <- function(y, archOrder = 1, garchOrder = 0, mode = "unitOmega") {

  # Rugarch specification
  if (mode == "unitOmega") {
    # Rugarch specification
    spec <- rugarch::ugarchspec(variance.model = list(garchOrder = c(archOrder, garchOrder), variance.targeting = 1),
                                mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
  } else {
    # Rugarch specification
    spec <- rugarch::ugarchspec(variance.model = list(garchOrder = c(archOrder, garchOrder)),
                                mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
  }
  # Safe GARCH fit
  safe_GARCH <- purrr::safely(rugarch::ugarchfit)
  # Fitted model
  model <- safe_GARCH(data = y, spec = spec)$result

  # Check convergence
  if (model@fit$convergence == 0) {
    # cli::cli_alert_success("GARCH routine converged!")
    params <- model@fit$ipars[,1][model@fit$ipars[,3]==1]
    # Extract ARCH parameters
    alpha <- params[stringr::str_detect(names(params), "alpha")]
    if (purrr::is_empty(alpha)){
      alpha <- c(alpha1 = 0)
    }
    # Extract GARCH parameters
    beta <- params[stringr::str_detect(names(params), "beta")]
    if (purrr::is_empty(beta)){
      beta <- c(beta1 = 0)
    }
    # Extract intercept
    omega <- params[names(params) == "omega"]
    # Store the log-likelihood
    loglik <- -model@fit$log.likelihoods
    # Robust standard errors
    std.errors <- model@fit$robust.se.coef
    if (archOrder == 1 & garchOrder == 0 & mode == "unitOmega"){
      names(model@fit$coef) <- c("alpha1", "omega")
    }
    names(std.errors) <- names(model@fit$coef)
    # Full std. errors
    std.errors_full <- c(omega = omega[[1]], alpha, beta)
    std.errors_full[1:length(std.errors_full)] <- NA
    std.errors_full[names(std.errors)] <- std.errors
  } else {
    cli::cli_alert_danger("GARCH routine do not converged!")
    return(NULL)
  }

  # Output
  structure(
    list(
      omega = omega,
      alpha = alpha,
      beta = beta,
      std.errors = std.errors_full,
      loglik = loglik
    )
  )
}


#' sGARCH Gaussian QMLE with rugarch or manual routine
#'
#' Fits an sGARCH(p, q) by Gaussian QMLE using the unconstrained parameterization
#' \code{zeta} and the mapping in \code{GARCH_params_to_phi()}. Returns both
#' unconstrained and constrained estimates, and the Hessian-based covariance of
#' the unconstrained parameters. (If a Jacobian consistent with the chosen
#' parameterization is supplied, delta-method SEs for constrained params can be computed.)
#'
#' @param x Numeric vector to fit.
#' @param archOrder Integer \eqn{p}.
#' @param garchOrder Integer \eqn{q}.
#' @inheritParams sGARCH_params_to_phi
#'
#' @return A list with elements:
#' \describe{
#'   \item{par_unconstrained}{Estimated \code{zeta}.}
#'   \item{par_constrained}{Estimated \code{phi} (omega, alphas, betas).}
#'   \item{vcov_unconstrained}{Hessian-based covariance matrix of \code{zeta}.}
#'   \item{se_unconstrained}{Standard errors of \code{zeta}.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(500)
#' fit <- sGARCH_robust_fit(x, archOrder = 1, garchOrder = 1)
#'
#' @export
sGARCH_robust_fit <- function(y, weights, archOrder = 1, garchOrder = 0, mode = "unitOmega", method = "rugarch") {

  res <- NULL
  if (missing(weights)){
    if (method == "rugarch") {
      res <- sGARCH_fit_rugarch(y, archOrder, garchOrder, mode)
    }
  }

  if (is.null(res) | method != "rugarch") {
    if (is.null(res) && method == "rugarch" & missing(weights)) cli::cli_alert_danger("GARCH routine do not converged!")
    model <- sGARCH_fit(y, weights, archOrder, garchOrder, mode)
    # Parameters
    params <- model$coef
    # Extract ARCH parameters
    alpha <- params[stringr::str_detect(names(params), "alpha")]
    if (purrr::is_empty(alpha)){
      alpha <- c(alpha1 = 0)
    }
    # Extract GARCH parameters
    beta <- params[stringr::str_detect(names(params), "beta")]
    if (purrr::is_empty(beta)){
      beta <- c(beta1 = 0)
    }
    # Extract intercept
    omega <- params[names(params) == "omega"]
    # Log-likelihoods
    loglik <- model$log.likelihoods
    # Robust standard errors
    std.errors <- model$std.errors
    names(std.errors) <- names(params)[params != 0]
    # Full std. errors
    std.errors_full <- c(omega = omega[[1]], alpha, beta)
    std.errors_full[1:length(std.errors_full)] <- NA
    std.errors_full[names(std.errors)] <- std.errors
    # Result
    res <- list(
      omega = omega,
      alpha = alpha,
      beta = beta,
      std.errors = std.errors_full,
      loglik = loglik
    )
  }
  return(res)
}

#' Forecast GARCH
#' @keywords GARCH
#' @note Version 1.0.0.
#' @rdname sGARCH_forecast_sigma2
#' @name sGARCH_forecast_sigma2
#' @export
#' @noRd
sGARCH_forecast_sigma2 <- function(h = 1, A, b, d, e1, S0, e_u2 = 1) {
  res <- .Call("sGARCH_forecast_sigma2_c",
               as.matrix(A),
               as.matrix(b),
               as.matrix(d),
               as.matrix(e1),
               as.matrix(S0),
               as.numeric(h),
               as.numeric(e_u2))
  names(res$E_S0) <- paste0("t+", 1:h)
  dplyr::tibble(h = 1:h, e_sigma2 = res$e_sigma2, E_S0 = res$E_S0)
}

#' Forecast GARCH
#' @keywords GARCH
#' @note Version 1.0.0.
#' @rdname sGARCH_forecast_sigma4
#' @name sGARCH_forecast_sigma4
#' @export
#' @noRd
sGARCH_forecast_sigma4 <- function(h = 1, A, b, d, e1, S0, E_S0, e_u2 = 1, e_u4 = 3) {
  res <- .Call("sGARCH_forecast_sigma4_c",
               as.matrix(A),
               as.matrix(b),
               as.matrix(d),
               as.matrix(e1),
               as.matrix(S0),
               as.list(E_S0),
               as.numeric(h),
               as.numeric(e_u2),
               as.numeric(e_u4))
  names(res$E2_S0) <- paste0("t+", 1:h)
  dplyr::tibble(h = 1:h, e_sigma4 = res$e_sigma4, E2_S0 = res$E2_S0)
}

#' Forecast GARCH
#' @keywords GARCH
#' @note Version 1.0.0.
#' @rdname sGARCH_forecast_covariance
#' @name sGARCH_forecast_covariance
#' @export
#' @noRd
sGARCH_forecast_covariance <- function(A, b, d, e1, E_S0, E2_S0, e_u2 = 1) {
  h <- length(E_S0)
  .Call("sGARCH_forecast_covariance_c",
        as.matrix(A),
        as.matrix(b),
        as.matrix(d),
        as.matrix(e1),
        as.list(E_S0),
        as.list(E2_S0),
        as.numeric(h),
        as.numeric(e_u2))
}

#' Forecast GARCH
#' @keywords GARCH
#' @note Version 1.0.0.
#' @rdname sGARCH_forecast
#' @name sGARCH_forecast
#' @export
#' @noRd
sGARCH_forecast <- function(h, A, b, d, e1, S0, e_u2, e_u4){
  # Forecast expectation
  E_S0 <- sGARCH_forecast_sigma2(h, A, b, d, e1, S0, e_u2)
  # Second moment
  E2_S0 <- sGARCH_forecast_sigma4(h, A, b, d, e1, S0, E_S0$E_S0, e_u2, e_u4)
  # Covariance matric
  Cv_S0 <- sGARCH_forecast_covariance(A, b, d, e1, E_S0$E_S0, E2_S0$E2_S0, e_u2)
  # Moments
  moments <- dplyr::tibble(h = 0:(h-1),
                    e_sigma2 = E_S0$e_sigma2,
                    e_sigma4 = E2_S0$e_sigma4,
                    v_sigma2 = diag(Cv_S0),
                    v_sigma = v_sigma2 / (4*e_sigma2),
                    e_sigma = e_sigma2^(1/2) - (1/8) * v_sigma2 / sqrt(e_sigma2)^3)
  Cv_S0[upper.tri(Cv_S0)] <- 0
  diag(Cv_S0) <- 0
  # Add covariances
  moments$cv_sigma2 <- 0
  moments$cv_sigma <- 0
  for(i in 1:nrow(Cv_S0)){
    moments$cv_sigma2[i] = list(Cv_S0[i,-i])
    moments$cv_sigma[i] = list(Cv_S0[i,-i] / (4 * sqrt(E_S0$e_sigma2[i] * E_S0$e_sigma2[-i])))
  }
  moments
}

