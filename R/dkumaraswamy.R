#' Kumaraswamy random variable
#'
#' Kumaraswamy density, distribution, quantile and random generator.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param a parameter `a > 0`.
#' @param b parameter `b > 0`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#'
# @references Kumaraswamy Distribution \href{https://en.wikipedia.org/wiki/Kumaraswamy_distribution}{W}.
#'
#' @examples
#' # Grid
#' x <- seq(0, 1, 0.01)
#'
#' # Density function
#' plot(x, dkumaraswamy(x, 0.2, 0.3), type = "l")
#' plot(x, dkumaraswamy(x, 2, 1.1), type = "l")
#'
#' # Distribution function
#' plot(x, pkumaraswamy(x, 2, 1.1), type = "l")
#'
#' # Quantile function
#' qkumaraswamy(0.2, 0.4, 1.4)
#' qkumaraswamy(pkumaraswamy(0.4, 2, 1.1),2, 1.1)
#'
#' # Random generator
#' rkumaraswamy(20, 0.4, 1.4)
#'
#' @name dkumaraswamy
#' @rdname dkumaraswamy
#' @aliases dkumaraswamy
#' @aliases pkumaraswamy
#' @aliases qkumaraswamy
#' @aliases rkumaraswamy
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dkumaraswamy <- function(x, a = 1, b = 1, log = FALSE){
  # Density function
  p <- a*b*x^(a - 1)*(1 - x^a)^(b - 1)
  # Ensure bounds
  p[x<0|x>1] <- 0
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname dkumaraswamy
pkumaraswamy <- function(q, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  # Distribution function
  p <- 1 - (1 - q^a)^b
  # Ensure bounds
  p[q<0|q>1] <- 0
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Log-probability
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname dkumaraswamy
qkumaraswamy <- function(p, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Ensure bounds
  p[p<=0] <- 0
  p[p>=1] <- 1
  # Quantiles
  q <- (1 - (1-p)^(1/b))^(1/a)
  return(q)
}

#' @export
#' @rdname dkumaraswamy
rkumaraswamy <- function(n, a = 1, b = 1){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated quantiles
  q <- qkumaraswamy(u, a = a, b = b, log.p = FALSE, lower.tail = TRUE)
  return(q)
}

