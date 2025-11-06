#' Riccati Root
#'
#' Compute the square root of a symmetric matrix.
#'
#' @param x squared and symmetric matrix.
#'
#' @examples
#' cv <- matrix(c(1, 0.3, 0.3, 1), nrow = 2, byrow = TRUE)
#' riccati_root(cv)
#'
#' @rdname riccati_root
#' @name riccati_root
#' @export
riccati_root <- function(x){
  dec <- eigen(x)
  e <- dec$vectors
  lam <- dec$values
  e %*% diag(sqrt(lam)) %*% t(e)
}

#' Make a matrix positive semi-definite
#'
#' @param x matrix, squared and symmetric.
#' @param neg_values numeric, the eigenvalues lower the zero will be substituted with this value.
#'
#' @examples
#' m <- matrix(c(2, 2.99, 1.99, 2), nrow = 2, byrow = TRUE)
#' makeSemiPositive(m)
#'
#' @rdname makeSemiPositive
#' @name makeSemiPositive
#' @export
makeSemiPositive <- function(x, neg_values = 1e-10){
  # Spectral decomposition
  dec <- eigen(x)
  # Eigenvectors
  e <- dec$vectors
  # Eigenvalues
  lam <- dec$values
  # Substitute negative values
  if (any(lam < 0)){
    lam[lam < 0] <- neg_values
    mat <- e %*% diag(lam) %*% t(e)
  }
  # Store informations
  attr(x, "index_neg_values") <- which(dec$values < 0)
  attr(x, "original_values") <- dec$values[dec$values < 0]
  return(x)
}



