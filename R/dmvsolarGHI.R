#' Bivariate PDF GHI
#' @param x vector of quantiles.
#' @param Ct clear sky radiation
#' @param alpha parameters `alpha > 0`.
#' @param beta parameters `beta > 0` and `alpha + beta < 1`.
#' @param joint_pdf_Yt joint density of Y1_t, Y2_t.
#'
#' @name dmvsolarGHI
#' @rdname dmvsolarGHI
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmvsolarGHI <- function(x, Ct, alpha, beta, joint_pdf_Yt){
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  z <- x
  z[,1] <- (1 - x[,1]/Ct[1] - alpha[1])/beta[1]
  z[,2] <- (1 - x[,2]/Ct[2] - alpha[2])/beta[2]
  u <- log(-log(z))
  # Compute denominator
  z_prod <- apply(z, 1, prod)
  # Denominator
  den <- prod(Ct*beta)*apply(z, 1, prod)*apply(log(z), 1, prod)
  # Mixture probabilities
  probs <- (1/den)*joint_pdf_Yt(u)
  probs[is.nan(probs)] <- 0
  return(probs)
}
