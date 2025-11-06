#' Density, distribution, quantile and random generator of Solar Radiation
#'
#' @param x,p Numeric vector of quantiles or probabilities.
#' @param Ct Numeric scalar, clear sky radiation
#' @param alpha Numeric scalar, parameter `alpha > 0`.
#' @param beta Numeric scalar, parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y Function, density of Y.
#' @param cdf_Y Function, distribution of Y.
#' @param log Logical, when `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p Logical, when `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail Logical, when `TRUE`, the default, the computed probabilities are \eqn{\mathbb{P}(X < x)}, otherwise \eqn{\mathbb{P}(X \ge x)}.
#'
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Y`. Then
#' the funtion `dsolarGHI` compute the density function of the following transformed random variable, i.e.
#' \deqn{R_t(y) = C(t) (1-\alpha-\beta \exp(-\exp(y)))}
#' where \eqn{R_t(y) \in [C(t)(1-\alpha-\beta), C(t)(1-\alpha)]}.
#' @examples
#' # Parameters
#' alpha = 0
#' beta = 0.9
#' Ct <- 7
#' # Grid of points
#' grid <- seq(Ct*(1-alpha-beta), Ct*(1-alpha), by = 0.01)
#'
#' # Density
#' dsolarGHI(5, Ct, alpha, beta, function(x) dnorm(x))
#' dsolarGHI(5, Ct, alpha, beta, function(x) dnorm(x, sd=2))
#' plot(grid, dsolarGHI(grid, Ct, alpha, beta, function(x) dnorm(x, mean = -1, sd = 0.9)), type="l")
#'
#' # Distribution
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) pnorm(x))
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) pnorm(x, sd=2))
#' plot(grid, psolarGHI(grid, Ct, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")
#'
#' # Quantile
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) pnorm(x))
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) pnorm(x, sd=2))
#'
#' # Random generator (I)
#' Ct <- Bologna$seasonal_data$Ct
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, alpha, beta, function(x) pnorm(x, sd=1.4)))
#' plot(1:366, GHI, type="l")
#'
#' # Random generator (II)
#' cdf <- function(x) pmixnorm(x, c(-0.8, 0.5), c(1.2, 0.7), c(0.3, 0.7))
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, 0.001, 0.9, cdf))
#' plot(1:366, GHI, type="l")
#' @rdname dsolarGHI
#' @aliases dsolarGHI
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsolarGHI <- function(x, Ct, alpha, beta, pdf_Y, log = FALSE, link = "invgumbel"){
  z_x <- (1 - x/Ct - alpha)/beta
  if (link == "invgumbel") {
    u_x <- suppressWarnings(log(-log(z_x)))
    den <- -Ct*beta*log(z_x^z_x)
  } else if (link == "gumbel") {
    u_x <- suppressWarnings(-log(-log(z_x)))
    den <- -Ct*beta*log(z_x^z_x)
  } else if (link == "logis") {
    u_x <- suppressWarnings(log(z_x/(1-z_x)))
    den <- Ct*beta*z_x*(1-z_x)
  } else if (link == "norm") {
    u_x <- suppressWarnings(qnorm(z_x))
    den <- Ct*beta*dnorm(qnorm(z_x))
  }
  # Numerator
  num <- pdf_Y(u_x)
  # Probability
  probs <- num/den
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  probs[is.nan(probs)] <- 0
  return(probs)
}

#' Distribution function for the GHI
#'
#' @rdname dsolarGHI
#' @export
psolarGHI  <- function(x, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE, link = "invgumbel"){
  z_x <- (1 - x/Ct - alpha)/beta
  if (link == "invgumbel") {
    u_x <- suppressWarnings(log(-log(z_x)))
    probs <- cdf_Y(u_x)
  } else if (link == "gumbel") {
    u_x <- suppressWarnings(-log(-log(z_x)))
    probs <- 1-cdf_Y(u_x)
  } else if (link == "logis") {
    u_x <- suppressWarnings(log(z_x/(1-z_x)))
    probs <- 1-cdf_Y(u_x)
  } else if (link == "norm") {
    u_x <- suppressWarnings(qnorm(z_x))
    probs <- 1-cdf_Y(u_x)
  }
  probs[x<=Ct*(1-alpha-beta)] <- 0
  probs[x>=Ct*(1-alpha)] <- 1
  # Lower tail
  if (!lower.tail) {
    probs <- 1 - probs
  }
  # Log-probability
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Quantile function for the GHI
#'
#' @rdname dsolarGHI
#' @export
qsolarGHI  <- function(p, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE, link = "invgumbel"){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for GHI
  interval = c(Ct*(1-alpha-beta), Ct*(1-alpha))
  # Density function
  cdf <- function(x) psolarGHI(x, Ct, alpha, beta, cdf_Y, link = link)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)

  return(x)
}

#' Random generator function for the GHI
#'
#' @rdname dsolarGHI
#' @export
rsolarGHI  <- function(n, Ct, alpha, beta, cdf_Y, link = "invgumbel"){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarGHI(u, Ct, alpha, beta, cdf_Y, link = link)
}
