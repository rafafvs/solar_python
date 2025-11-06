#' Gumbel random variable
#'
#' Gumbel density, distribution, quantile and random generator.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#@references Gumbel distribution [\href{https://en.wikipedia.org/wiki/Gumbel_distribution}{W}].
#'
#' @examples
#' # Grid
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function
#' p <- dgumbel(x, location = 0, scale = 1)
#' plot(x, p, type = "l")
#'
#' # Distribution function
#' p <- pgumbel(x, location = 0, scale = 1)
#' plot(x, p, type = "l")
#'
#' # Quantile function
#' qgumbel(0.1)
#' pgumbel(qgumbel(0.1))
#'
#' # Random Numbers
#' rgumbel(1000)
#' plot(rgumbel(1000), type = "l")
#'
#' @name dgumbel
#' @rdname dgumbel
#' @aliases dgumbel
#' @aliases pgumbel
#' @aliases qgumbel
#' @aliases rgumbel
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dgumbel <- function(x, location = 0, scale = 1, log = FALSE){
  # Standardized values
  z <- (x-location)/scale
  # Density
  p <- (1/scale)*exp(-(z + exp(-z)))
  # Log probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}


#' @export
#' @rdname dgumbel
pgumbel <- function(q, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (q - location)*scale
  # Distribution
  p <- exp(-exp(-z))
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
#' @rdname dgumbel
qgumbel <- function(p, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE) {
  # Log probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Quantiles
  q <- location - scale*log(-log(p))
  return(q)
}


#' @export
#' @rdname dgumbel
rgumbel <- function(n, location = 0, scale = 1){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated values
  q <- qgumbel(u, location, scale, log.p = FALSE, lower.tail = FALSE)
  return(q)
}

