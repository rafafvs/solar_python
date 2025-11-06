#' Skewed Normal random variable
#'
#' Skewed Normal density, distribution, quantile and random generator.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape skewness parameter.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
# @references Skewed Normal Distribution [\href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{W}].
#'
#' @examples
#' # Grid of points
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function
#' # right tailed
#' plot(x, dsnorm(x, shape = 1.9), type = "l")
#' # left tailed
#' plot(x, dsnorm(x, shape = -1.9), type = "l")
#'
#' # Distribution function
#' plot(x, psnorm(x, shape = 4.9), type = "l")
#' plot(x, psnorm(x, shape = -4.9), type = "l")
#'
#' # Quantile function
#' dsnorm(0.1, shape = 4.9)
#' dsnorm(0.1, shape = -4.9)
#' psnorm(qsnorm(0.9, shape = 3), shape = 3)
#'
#' # Random generator
#' set.seed(1)
#' plot(rsnorm(100, shape = 5), type = "l")
#'
#' @name dsnorm
#' @rdname dsnorm
#' @aliases dsnorm
#' @aliases psnorm
#' @aliases qsnorm
#' @aliases rsnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsnorm <- function(x, location = 0, scale = 1, shape = 0, log = FALSE){
  # Standardized values
  z <- (x - location)/scale
  # Probabilities
  p <- 2*dnorm(z)*pnorm(shape*z)
  # Log-probabilities
  if (log) {
    p <- base::log(p)
  }
  return(p)
}


#' @rdname dsnorm
#' @export
psnorm <- function(q, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (q - location)/scale
  # Distribution function
  cdf <- CDF(dsnorm, location = location, scale = scale, shape = shape, lower=-Inf, log = FALSE)
  # Cumulated probabilities
  p <- cdf(z)
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Log-probabilities
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}


#' @rdname dsnorm
#' @export
qsnorm <- function(p, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE) {
  # Distribution function
  cdf <- function(x) psnorm(x, location = location, scale = scale, shape = shape)
  # Quantile function
  quantile_numeric <- Quantile(cdf, interval = c(-location - scale*10, location + scale*10))
  # Quantiles
  q <- quantile_numeric(p, log.p = log.p, lower.tail = lower.tail)
  return(q)
}


#' @rdname dsnorm
#' @export
rsnorm <- function(n, location = 0, scale = 1, shape = 0){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Quantiles
  q <- qsnorm(u, location, scale, shape, log.p = FALSE, lower.tail = TRUE)
  return(q)
}
