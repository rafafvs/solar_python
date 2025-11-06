#' Moments of a Gaussian Mixture
#'
#' Compute the first four moments, central moments, mean, variance, skewness and kurtosis
#' for a Gaussian Mixture with K components.
#'
#' @examples
#' means = c(-0.9, 0.5)
#' sd = c(0.4,1)
#' alpha = c(0.5, 0.5)
#' GM_moments(means, sd, alpha)
#' GM_moments(c(-0.8, 1.8), c(0.4,1), c(0.5, 0.5))
#'
#' @details
#' The non-central moments are denoted as `m1`, `m2`, `m3`, `m4`. The central moments as `mu2`, `mu3`, `mu4`.
#' @rdname GM_moments
#' @name GM_moments
#' @references https://en.wikipedia.org/wiki/Normal_distribution#Moments
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
#' @noRd
GM_moments <- function(means, sd, alpha){
  # First moment (Expectation)
  m1 <- sum(means * alpha)
  # Non-central moments
  # Second moment
  m2 <- sum((means^2 + sd^2) * alpha)
  # Third moment
  m3 <- sum((3 * means * sd^2 + means^3) * alpha)
  # Fourth moment
  m4 <- sum((3 * sd^4 + 6 * means^2 * sd^2 + means^4) * alpha)
  # Fifth moment
  # m5 <- sum((means^5 + 10 * means^3 * sd^2 + 15 * means * sd^4) * alpha)
  # Sixth moment
  # m6 <- sum((means^6 + 15 * means^4 * sd^2 + 45 * means^2 * sd^4 + 15 * sd^6)*alpha)
  # Central-moments
  delta_k <- means - m1
  # Second moment
  mu2 <- sum((delta_k^2 + sd^2) * alpha)
  # Third moment
  mu3 <- sum((3 * delta_k * sd^2 + delta_k^3) * alpha)
  # Fourth moment
  mu4 <- sum((3 * sd^4 + 6 * delta_k^2 * sd^2 + delta_k^4) * alpha)
  # Fifth moment
  # mu5 <- sum((delta_k^5 + 10 * delta_k^3 * sd^2 + 15 * delta_k * sd^4) * alpha)
  # Sixth moment
  # mu6 <- sum((delta_k^6 + 15 * delta_k^4 * sd^2 + 45 * delta_k^2 * sd^4 + 15 * sd^6)*alpha)

  # Derived statistics
  # Variance
  v_x <- mu2
  # Skewness
  sk_x <- mu3 / (v_x^(3/2))
  # Excess Kurtosis
  kt_x <- mu4 / (v_x^2) - 3

  dplyr::tibble(
    m1 = m1,
    m2 = m2,
    m3 = m3,
    m4 = m4,
    mu2 = mu2,
    mu3 = mu3,
    mu4 = mu4,
    mean = m1,
    variance = v_x,
    skewness = sk_x,
    kurtosis = kt_x
  )
}

#' Match the first three moments of a Gaussian Mixture
#'
#' @param d Numeric, distance between the two means.
#' @param m1 Numeric, first target moment.
#' @param m2 Numeric, second target moment.
#' @param m3 Numeric, third target moment.
#' @param p Numeric, probability.
#'
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
#' @noRd
GM_moments_match <- function(d, m1 = 0, m2 = 1, m3 = 0, p = 0.5){
  means <- c(mu1 = 0, mu2 = 0)
  sigma <- c(sd1 = 0, sd2 = 0)
  probs <- c(p1 = p, p2 = 1-p)
  # First moment
  means[1] <- m1 + (1-p) * d
  means[2] <- m1 - p * d
  # Second moment
  delta <- (m3 - p * (1 - p) * ((1 - p)^2 - p^2) * d^3) / (3 * p * (1-p) * d)
  sigma[2] <- m2 - p * (1 - p) * d^2 - p * delta
  sigma[1] <- sigma[2] + delta
  sigma <- sqrt(sigma)

  structure(
    list(
      means = means,
      sigma = sigma,
      probs = probs
    )
  )
}

#' Compute the log-likelihood of a Gaussian Mixture
#'
#' @param means description
#' @param sd description
#' @param alpha description
#' @examples
#' GM_loglik(c(-0.8, 0.8), c(0.4,1), c(0.5, 0.5), rnorm(100))
#'
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
#' @noRd
GM_loglik <- function(means, sd, alpha, x){
  # Log-likelihood
  if (!missing(x)) {
    lik <- dmixnorm(x, means, sd, alpha)
    loss <- sum(log(lik))
  } else {
    loss <- 0
  }
  return(loss)
}

#' Fit a Gaussian Mixture model
#'
#' @param x Numeric, vector on which the model will be fitted.
#' @param method Character, method use to fit the GM model, can be `mclust` or `mixtools`.
#' @param components Integer, number of components for the mixture.
#'
#' @examples
#' x <- rmixnorm(1000, c(0, 0), c(0.3, 2), c(0.4, 0.6))
#' # Fitted parameters without contraints
#' params <- GM_fit(x$X)
#' GM_moments(params$means, params$sd, params$p)
#'
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
#' @noRd
GM_fit <- function(x, method = c("mclust", "mixtools"), components = 2, maxit = 30000, maxrestarts = 500){
  clust <- NULL
  method <- match.arg(method, choices = c("mclust", "mixtools"))
  # Fitted parameters with mclust
  if (method == "mclust") {
    clust <- mclust::Mclust(x, G = components, modelNames = c("V"), verbose = FALSE)
    if (!is.null(clust)){
      # Estimated parameters
      means <- clust$parameters$mean
      sd <- sqrt(clust$parameters$variance$scale)
      p <- clust$parameters$pro
    }
  }
  # Fitted parameters with mixtools
  if (method == "mixtools" | is.null(clust)) {
    quiet_EM <- purrr::quietly(mixtools::normalmixEM)
    clust <- quiet_EM(x, maxit = maxit, k = components, maxrestarts = maxrestarts)$result
    # Estimated parameters
    means <- clust$mu
    sd <- clust$sigma
    p <- clust$lambda
  }
  idx <- order(means)
  means = means[idx]
  sd = sd[idx]
  p = p[idx]
  # Assign a name to the parameters
  names(means) <- paste0("mu", 1:components)
  names(sd) <- paste0("sd", 1:components)
  names(p) <- paste0("p", 1:components)

  structure(
    list(
      means = means,
      sd = sd,
      p = p
    )
  )
}

#' Fit a Gaussian Mixture with two components under moments constraints.
#'
#' @param x Numeric, vector on which the model will be fitted.
#' @param start_params Numeric vector of parameters
#' @param mu_target Numeric, scalar. Target mean of the mixture. If missing will be not applied any contraint.
#' @param var_target Numeric, scalar. Target variance of the mixture. If missing will be not applied any contraint.
#'
#' @examples
#' x <- rmixnorm(1000, c(0, 0), c(0.3, 2), c(0.4, 0.6))
#' # Fitted parameters without contraints
#' params <- GM_fit(x$X)
#' GM_moments(params$means, params$sd, params$p)
#'
#' # Fit the parameters with moments contraints
#' start_params <- unlist(purrr::flatten(params)[-6])
#' # Match a certain mean
#' match_params <- GM_fit_moments_match(x$X, start_params, mu_target = -0.09)
#' GM_moments(match_params$means, match_params$sd, match_params$p)
#' # Match a certain variance
#' match_params <- GM_fit_moments_match(x$X, start_params, var_target = 2.4)
#' GM_moments(match_params$means, match_params$sd, match_params$p)
#'
#' # Match a certain mean and variance
#' match_params <- GM_fit_moments_match(x$X, start_params, mu_target = -0.09, var_target = 2.4)
#' GM_moments(match_params$means, match_params$sd, match_params$p)
#'
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
#' @noRd
GM_fit_moments_match <- function(x, start_params, mu_target = NA, var_target = NA, maxit = 1000, abstol = 1e-5){
  # Objective: weighted log-likelihood
  objective <- function(par){
    mu1 <- par[1]
    mu2 <- par[2]
    sd1 <- par[3]
    sd2 <- par[4]
    p1 <- par[5]
    p2 <- 1 - par[5]
    # Log-likelihood
    loglik <- sum(log(p1 * dnorm(x, mu1, sd1) + p2 * dnorm(x, mu2, sd2)), na.rm = TRUE)
    return(-loglik)
  }

  # Equality constraints: match mean and variance
  constraints <- function(par){
    mu1 <- par[1]
    mu2 <- par[2]
    sd1 <- par[3]
    sd2 <- par[4]
    p1 <- par[5]
    p2 <- 1 - par[5]
    # Moments
    e_x <- p1 * mu1 + p2 * mu2
    e_x2 <- p1 * (sd1^2 + mu1^2) + p2 * (sd2^2 + mu2^2)
    v_x <- e_x2 - e_x^2

    if (is.na(mu_target) & !is.na(var_target)){
      loss <- c(0, v_x - var_target)
    } else if (!is.na(mu_target) & is.na(var_target)){
      loss <- c(e_x - mu_target, 0)
    } else {
      loss <- c(e_x - mu_target, v_x - var_target)
    }
    return(loss)
  }
  # Solver with constraints
  res <- nloptr::nloptr(
    x0 = start_params,
    eval_f = objective,
    eval_g_eq = constraints,
    # Lower bounds
    lb = c(-3, -3, 0.01, 0.01, 0.01),
    # Upper bounds
    ub = c(3, 3, 3, 3, 0.99),
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = abstol, maxeval = maxit)
  )
  structure(
    list(
      means = res$solution[1:2],
      sd = res$solution[3:4],
      p = c(res$solution[5], 1-res$solution[5])
    )
  )
}
