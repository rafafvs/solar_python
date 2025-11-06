#' Truncated Normal random variable
#'
#' Truncated Normal density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param a lower bound.
#' @param b upper bound.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @examples
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function
#' p <- dtnorm(x, mean = 0, sd = 1, a = -1)
#' plot(x, p, type = "l")
#'
#' # Distribution function
#' p <- ptnorm(x, mean = 0, sd = 1, b = 1)
#' plot(x, p, type = "l")
#'
#' # Quantile function
#' dtnorm(0.1)
#' ptnorm(qtnorm(0.1))
#'
#' # Random Numbers
#' rtnorm(1000)
#' plot(rtnorm(100, mean = 0, sd = 1, a = 0, b = 1), type = "l")
#'
#' @name dtnorm
#' @rdname dtnorm
#' @aliases dtnorm
#' @aliases ptnorm
#' @aliases qtnorm
#' @aliases rtnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dtnorm <- function(x, mean = 0, sd = 1, a = -3, b = 3, log = FALSE){
  x[x < a | x > b] <- NA
  z <- (x - mean)/sd
  p <- (1/sd)*(dnorm(z)/(pnorm((b - mean)/sd) - pnorm((a - mean)/sd)))
  p[which(is.na(p))] <- 0

  if (log){
    return(base::log(p))
  }
  return(p)
}


#' @rdname dtnorm
#' @export
ptnorm <- function(x, mean = 0, sd = 1, a = -3, b = 3, log.p = FALSE, lower.tail = TRUE){
  z <- (x - mean)/sd
  p <- c()
  for(i in 1:length(z)){
    p[i] <- integrate(dtnorm, lower=-Inf, upper = z[i], a = a, b = b)$value
  }

  if (!lower.tail) {
    p <- 1 - p
  }

  if (log.p) {
    return(base::log(p))
  }
  return(p)
}


#' @export
#' @rdname dtnorm
qtnorm <- function(p, mean = 0, sd = 1, a = -3, b = 3, log.p = FALSE, lower.tail = TRUE) {

  if (log.p) {
    p <- exp(p)
  }

  loss <- function(x, p) {
    p_hat <- ptnorm(x, mean, sd, a = a, b = b, lower.tail = lower.tail)
    (p_hat - p)^2
  }
  x <- c()
  for(i in 1:length(p)){
    x[i] <- suppressWarnings(optim(par = (a+b)/2, loss, p = p[i])$par)
  }
  return(x)
}


#' @export
#' @rdname dtnorm
rtnorm <- function(n, mean = 0, sd = 1, a = -100, b = 100){
  u <- runif(n, min = 0, max = 1)
  x <- qtnorm(u, mean = mean, sd = sd, a = a, b = b)
  x[x < a] = a
  x[x > b] = b
  x
}

