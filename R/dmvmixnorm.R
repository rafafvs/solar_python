#' Multivariate Gaussian mixture random variable
#'
#' Multivariate Gaussian mixture density, distribution, quantile and random generator.
#'
#' @examples
#' # Means components
#' mean_1 = c(-1.8,-0.4)
#' mean_2 = c(0.6, 0.5)
#' # Dimension of the random variable
#' j = length(mean_1)
#' # Matrix of means
#' means = matrix(c(mean_1, mean_2), j,j, byrow = TRUE)
#'
#' # Variance components
#' var_1 = c(1,1.4)
#' var_2 = c(1.3, 1.2)
#' # Matrix of variances
#' sigma2 = matrix(c(var_1, var_2), j,j, byrow = TRUE)
#'
#' # Correlations
#' rho <- c(rho_1 = 0.2, rho_2 = 0.3)
#'
#' # Probability for each component
#' p <- c(0.4, 0.6)
#'
#' x <- matrix(c(0.1,-0.1), nrow = 1)
#' dmvmixnorm(x, means, sigma2, p, rho)
#' pmvmixnorm(x, means, sigma2, p, rho)
#' qmvmixnorm(0.35, means, sigma2, p, rho)
#' @name dmvmixnorm
#' @rdname dmvmixnorm
#' @aliases dmvmixnorm
#' @aliases pmvmixnorm
#' @aliases qmvmixnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), log = FALSE){
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  # Number of observations
  n <- nrow(x)
  # Number of components
  k <- nrow(means)
  # Number of variables
  j <- ncol(means)
  # Covariance matrix
  cv_ <- list()
  for(s in 1:k){
    cv_k <- diag(sigma2[s,])
    cv_k[upper.tri(cv_k)] <- cv_k[lower.tri(cv_k)] <- rho[1]*prod(sigma2[s,])
    cv_ [[s]] <- cv_k
  }
  # Density
  probs <- rep(0, n)
  for(i in 1:n){
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*mvtnorm::dmvnorm(x[i,], mean = means[s,], sigma  = cv_[[s]])
    }
  }
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmvmixnorm
#' @export
pmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), lower = -Inf, log.p = FALSE){
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  # Number of observations
  n <- nrow(x)
  # Number of components
  k <- nrow(means)
  # Number of variables
  j <- ncol(means)
  # Covariance matrix
  cv_ <- list()
  for(s in 1:k){
    cv_k <- diag(sigma2[s,])
    cv_k[upper.tri(cv_k)] <- cv_k[lower.tri(cv_k)] <- rho[1]*prod(sigma2[s,])
    cv_ [[s]] <- cv_k
  }
  # Density
  probs <- rep(0, n)
  for(i in 1:n){
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*mvtnorm::pmvnorm(lower = lower, upper = x[i,], mean = means[s,], sigma = cv_[[s]])
    }
  }
  # Log-probability
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmvmixnorm
#' @export
qmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), log.p = FALSE){
  # Log-probabilities
  probs <- x
  if (log.p) {
    probs <- base::exp(probs)
  }
  # Distribution function
  cdf <- function(x) pmvmixnorm(x, means = means, sigma2 = sigma2, p = p, rho = rho)

  # Empirical quantile optimizer
  quantile_ <- function(x, prob){(cdf(x) - prob)^2}
  # Quantiles
  q <- matrix(0, nrow = length(probs), ncol = ncol(means))
  for(i in 1:length(probs)){
    q[i,] <- suppressWarnings(optim(par = rep(0, ncol(means)), quantile_, prob = probs[i])$par)
  }
  return(q)
}
