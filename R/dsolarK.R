#' Clearness index random variable
#'
#' Clearness index density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param alpha parameter `alpha > 0`.
#' @param beta parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y density function of Y.
#' @param cdf_Y distribution function of Y.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Y`. Then
#' the funtion `dsolarK` compute the density function of the following transformed random variable, i.e.
#' \deqn{K(Y) = 1-\alpha-\beta \exp(-\exp(Y))}
#' where \eqn{K(Y) \in [1-\alpha-\beta, 1-\alpha]}.
#' @examples
#' # Parameters
#' alpha = 0.001
#' beta = 0.9
#' # Grid of points
#' grid <- seq(1-alpha-beta, 1-alpha, length.out = 50)[-50]
#'
#' # Density
#' dsolarK(0.4, alpha, beta, function(x) dnorm(x))
#' dsolarK(0.4, alpha, beta, function(x) dnorm(x, sd = 2))
#' plot(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = 0.2)), type="l")
#'
#' # Distribution
#' psolarK(0.493, alpha, beta, function(x) pnorm(x))
#' psolarK(0.493, alpha, beta, function(x) pnorm(x, sd = 2))
#' plot(grid, psolarK(grid, alpha, beta, function(x) pt(0.2*x, 3)), type="l")
#' plot(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")
#'
#' # Quantile
#' qsolarK(c(0.05, 0.95), alpha, beta, function(x) pnorm(x))
#' qsolarK(c(0.05, 0.95), alpha, beta, function(x) pnorm(x, sd = 2))
#'
#' # Random generator (I)
#' Kt <- rsolarK(366, alpha, beta, function(x) pnorm(x, sd = 1.3))
#' plot(1:366, Kt, type="l")
#'
#' # Random generator (II)
#' pdf <- function(x) pmixnorm(x, c(-1.8, 0.8), c(0.5, 0.7), c(0.6, 0.4))
#' Kt <- rsolarK(36, alpha, beta, pdf)
#' plot(1:36, Kt, type="l")
#' @rdname dsolarK
#' @aliases dsolarK
#' @aliases psolarK
#' @aliases qsolarK
#' @aliases rsolarK
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsolarK  <- function(x,  alpha, beta, pdf_Y, log = FALSE){
  z_x <- (1 - x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Y(u_x))/(beta*log(z_x^z_x))
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Distribution function for the Clearness index
#'
#' @rdname dsolarK
#' @export
psolarK  <- function(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  z_x <- (1 - x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- cdf_Y(u_x)
  probs[x<=(1-alpha-beta)] <- 0
  probs[x>=(1-alpha)] <- 1
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

#' Quantile function for the Clearness index
#'
#' @rdname dsolarK
#' @export
qsolarK  <- function(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for Clearness Index
  interval <- c(1-alpha-beta, 1-alpha)
  # Distribution function
  cdf <- function(x) psolarK(x, alpha, beta, cdf_Y)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)
  return(x)
}

#' Random generator function for the Clearness index
#'
#' @rdname dsolarK
#' @export
rsolarK  <- function(n, alpha, beta, cdf_Y){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarK(u, alpha, beta, cdf_Y)
}
