#' Sugeno-distorted density and distribution
#'
#' Compute the Sugeno-distorted distribution for a given distribution or density.
#'
#' @param cdf Function, distribution function.
#' @param pdf Function, density function.
#' @param lambda Numeric, distortion parameter.
#' @examples
#' # Distribution and density
#' cdf <- function(x) pnorm(x)
#' pdf <- function(x) dnorm(x)
#' x <- seq(-4, 4, 0.01)
#' plot(x, psugeno(cdf, lambda = -0.2)(x))
#' plot(x, dsugeno(pdf, cdf, lambda = -0.2)(x))
#'
#' @rdname dsugeno
#' @name dsugeno
#' @aliases dsugeno
#' @aliases psugeno
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsugeno <- function(pdf, cdf, lambda = 0){
  function(x) {
    den <- pdf(x) * (1 + lambda * (1-cdf(x))) - cdf(x) * (-lambda * pdf(x))
    num <- (1 + lambda *(1-cdf(x)))^2
    den/num
  }
}

#' @rdname dsugeno
#' @export
psugeno <- function(cdf, lambda = 0){
  function(x) {
    cdf(x) / (1 + lambda * (1 - cdf(x)))
  }
}

#' Sugeno upper and lower parameters.
#'
#' @param lambda Numeric, distortion parameter.
#' @rdname sugeno_bounds
#' @name sugeno_bounds
#' @export
sugeno_bounds <- function(lambda){
  # Identify the positive and negative parameters
  par <- list(pos = 0, bar = 0, neg = 0)
  if (lambda >= 0) {
    par[["pos"]] <- lambda
    par[["neg"]] <- -lambda / (1 + lambda)
  } else {
    par[["pos"]] <- -lambda / (1 + lambda)
    par[["neg"]] <- lambda
  }
  # Average parameter
  par[["bar"]] <- 0.5*(par[["pos"]] + par[["neg"]])
  return(par)
}
