#' Esscher-distorted density and distribution of a Gaussian Mixture
#'
#' @inheritParams dmixnorm
#' @param theta Numeric, distortion parameter.
#'
#' @examples
#' library(ggplot2)
#' grid <- seq(-5, 5, 0.01)
#' # Density
#' pdf_1 <- desscherMixture(mean = c(-3, 3), theta = 0)(grid)
#' pdf_2 <- desscherMixture(mean = c(-3, 3), theta = -0.5)(grid)
#' pdf_3 <- desscherMixture(mean = c(-3, 3), theta = 0.5)(grid)
#' ggplot()+
#'  geom_line(aes(grid, pdf_1), color = "black")+
#'  geom_line(aes(grid, pdf_2), color = "green")+
#'  geom_line(aes(grid, pdf_3), color = "red")
#' # Distribution
#' cdf_1 <- pesscherMixture(mean = c(-3, 3), theta = 0)(grid)
#' cdf_2 <- pesscherMixture(mean = c(-3, 3), theta = -0.2)(grid)
#' cdf_3 <- pesscherMixture(mean = c(-3, 3), theta = 0.2)(grid)
#' ggplot()+
#'   geom_line(aes(grid, cdf_1), color = "black")+
#'   geom_line(aes(grid, cdf_2), color = "green")+
#'   geom_line(aes(grid, cdf_3), color = "red")
#'
#' @rdname desscherMixture
#' @aliases desscherMixture
#' @aliases pesscherMixture
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
desscherMixture <- function(mean = c(0,0), sd = c(1,1), alpha = c(0.5, 0.5), theta = 0){
  # Distorted parameters
  params <- solarEsscher_GM_mixture(mean, sd, alpha, theta)
  # Mixture pdf
  pdf <- function(x) dmixnorm(x, params$mean, params$sd, params$alpha)
  # Esscher pdf
  function(x, log = FALSE){
    probs <- pdf(x)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' @rdname desscherMixture
#' @export
pesscherMixture <- function(mean = c(0,0), sd = c(1,1), alpha = c(0.5, 0.5), theta = 0){
  # Distorted parameters
  params <- solarEsscher_GM_mixture(mean, sd, alpha, theta)
  # Mixture pdf
  cdf <- function(x) pmixnorm(x, params$mean, params$sd, params$alpha)
  # Esscher pdf
  function(x, log = FALSE){
    probs <- cdf(x)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}
