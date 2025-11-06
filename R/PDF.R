#' Density, distribution and quantile function
#'
#' Return a function of `x` given the specification of a function of `x`.
#'
#' @param .f density function
#' @param cdf cumulative distribution function.
#' @param lower lower bound for integration (CDF).
#' @param interval lower and upper bounds for unit root (Quantile).
#' @param ... other parameters to be passed to `.f`.
#'
#' @examples
#' # Density
#' pdf <- PDF(dnorm, mean = 0.3, sd = 1.3)
#' pdf(3)
#' dnorm(3, mean = 0.3, sd = 1.3)
#' # Distribution
#' cdf <- CDF(dnorm, mean = 0.3, sd = 1.3)
#' cdf(3)
#' pnorm(3, mean = 0.3, sd = 1.3)
#' # Numeric quantile function
#' pnorm(Quantile(pnorm)(0.9))
#' @name PDF
#' @rdname PDF
#' @aliases PDF
#' @aliases CDF
#' @aliases Quantile
#' @export
PDF <- function(.f, ...){
  function(x, log = FALSE){
    probs <- .f(x, ...)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' @rdname PDF
#' @export
CDF <- function(.f, lower = -Inf, ...){
  # Density
  pdf <- PDF(.f, ...)
  # Distribution
  cdf <- function(x, pdf) integrate(pdf, lower = lower, upper = x)$value
  # Distribution function
  function(x,  lower.tail = TRUE, log.p = FALSE) {
    probs <- purrr::map_dbl(x, ~cdf(.x, pdf))
    # Lower tail
    if (!lower.tail) {
      probs <- 1 - probs
    }
    # Log-probabilities
    if (log.p) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' @rdname PDF
#' @export
Quantile <- function(cdf, interval = c(-100, 100)){

  # Find the quantile numerically
  quantile_root <- function(p, cdf, interval){
    uniroot(function(x) cdf(x) - p,
            interval = interval,
            tol = 10^{-16})$root
  }
  # Quantile function
  quantile_numeric <- function(p) purrr::map_dbl(p, ~quantile_root(.x, cdf, interval))

  function(p, log.p = FALSE, lower.tail = TRUE){
    probs <- p
    # Log probability
    if (log.p) {
      probs <- exp(probs)
    }
    # Lower tail
    if (!lower.tail) {
      probs <- 1 - probs
    }
    # Quantiles
    x <- quantile_numeric(probs)
    return(x)
  }
}





