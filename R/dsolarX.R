#' Solar risk driver random variable
#'
#' Solar risk driver density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param alpha parameter `alpha > 0`.
#' @param beta parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y density of Y.
#' @param cdf_Y distribution function of Y.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Y`. Then
#' the funtion `dsolarX` compute the density function of the following transformed random variable, i.e.
#' \deqn{X(Y) = \alpha+\beta \exp(-\exp(Y))}
#' where \eqn{X(Y) \in [\alpha, \alpha+\beta]}.
#' @examples
#' # Parameters
#' alpha = 0.001
#' beta = 0.9
#' # Grid of points
#' grid <- seq(alpha, alpha+beta, length.out = 50)[-50]
#'
#' # Density
#' dsolarX(0.4, alpha, beta, function(x) dnorm(x))
#' dsolarX(0.4, alpha, beta, function(x) dnorm(x, sd = 2))
#' plot(grid, dsolarX(grid, alpha, beta, function(x) dnorm(x, sd = 0.2)), type="l")
#'
#' # Distribution
#' psolarX(0.493, alpha, beta, function(x) pnorm(x))
#' dsolarX(0.493, alpha, beta, function(x) pnorm(x, sd = 2))
#' plot(grid, psolarX(grid, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")
#'
#' # Quantile
#' qsolarX(c(0.05, 0.95), alpha, beta, function(x) pnorm(x))
#' qsolarX(c(0.05, 0.95), alpha, beta, function(x) pnorm(x, sd = 1.3))
#'
#' # Random generator (I)
#' set.seed(1)
#' Kt <- rsolarX(366, alpha, beta, function(x) pnorm(x, sd = 0.8))
#' plot(1:366, Kt, type="l")
#'
#' # Random generator (II)
#' cdf <- function(x) pmixnorm(x, c(-1.8, 0.9), c(0.5, 0.7), c(0.6, 0.4))
#' Kt <- rsolarX(366, alpha, beta, cdf)
#' plot(1:366, Kt, type="l")
#' @rdname dsolarX
#' @aliases dsolarX
#' @aliases psolarX
#' @aliases qsolarX
#' @aliases rsolarX
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsolarX  <- function(x, alpha, beta, pdf_Y, log = FALSE){
  z_x <- (x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Y(u_x))/(beta*log(z_x^z_x))
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Distribution function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
psolarX  <- function(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){

  z_x <- (x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- 1 - cdf_Y(u_x)
  probs[x<=alpha] <- 0
  probs[x>=(alpha+beta)] <- 1
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

#' Quantile function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
qsolarX  <- function(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for solar risk driver
  interval <- c(alpha, alpha + beta)
  # Density function
  cdf <- function(x) psolarX(x, alpha, beta, cdf_Y)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)
  return(x)
}

#' Random generator function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
rsolarX  <- function(n, alpha, beta, cdf_Y){
  # Simulated grades
  u <- runif(n, 0, 1)
  # Simulated random variable
  qsolarX(u, alpha, beta, cdf_Y)
}
