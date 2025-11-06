#' Gaussian mixture random variable
#'
#' Gaussian mixture density, distribution, quantile and random generator.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mean vector of means parameters.
#' @param sd vector of std. deviation parameters.
#' @param alpha vector of probability parameters for each component.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{\mathbb{P}(X < x)}, otherwise \eqn{\mathbb{P}(X \ge x)}.
#'
# @references Mixture Models [\href{https://en.wikipedia.org/wiki/Mixture_model}{W}].
#'
#' @examples
#' # Parameters
#' mean = c(-3,0,3)
#' sd = rep(1, 3)
#' alpha = c(0.2, 0.3, 0.5)
#' # Density function
#' dmixnorm(3, mean, sd, alpha)
#' # Distribution function
#' dmixnorm(c(1.2, -3), mean, sd, alpha)
#' # Quantile function
#' qmixnorm(0.2, mean, sd, alpha)
#' # Random generator
#' rmixnorm(1000, mean, sd, alpha)
#'
#' @rdname dmixnorm
#' @name dmixnorm
#' @aliases dmixnorm
#' @aliases pmixnorm
#' @aliases qmixnorm
#' @aliases rmixnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmixnorm <- function(x, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), log = FALSE){
  # Density
  p <- c()
  for(i in 1:length(x)){
    p[i] <- 0
    for(s in 1:length(mean)){
      p[i] <- p[i] + alpha[s]*dnorm(x[i], mean = mean[s], sd = sd[s])
    }
  }
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' @rdname dmixnorm
#' @export
pmixnorm <- function(q, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), lower.tail = TRUE, log.p = FALSE){
  # Distribution
  p <- c()
  for(i in 1:length(q)){
    p[i] <- 0
    for(s in 1:length(mean)){
      p[i] <- p[i] + alpha[s]*pnorm(q[i], mean = mean[s], sd = sd[s], lower.tail = lower.tail)
    }
  }
  # Log-probability
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}

#' @rdname dmixnorm
#' @export
qmixnorm <- function(p, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), lower.tail = TRUE, log.p = FALSE) {
  # Log-probabilities
  if (log.p) {
    p <- base::exp(p)
  }
  # Distribution function
  cdf <- function(x) pmixnorm(x, mean, sd, alpha, lower.tail = lower.tail)
  # Empirical quantile
  quantile_numeric <- Quantile(cdf, interval = c(min(mean) - max(sd)*10, max(mean) + max(sd)*10))
  # Quantiles
  q <- quantile_numeric(p)
  #q <- p
  #q[p <= 0] <- -Inf
  #q[p >= 1] <- Inf
  #q[p > 0 & p < 1] <- quantile_numeric(p[p > 0 & p < 1])
  return(q)
}

#' @rdname dmixnorm
#' @export
rmixnorm <- function(n, mean = rep(0, 3), sd = rep(1, 3), alpha = rep(1/3, 3)){
  # Number of components
  k <- length(mean)
  X <- matrix(NA, nrow = n, ncol = k)
  B <- matrix(0, nrow = n, ncol = k)
  index <- 1:n
  for(s in 1:k){
    # Simulated bernoulli
    if (s == k) {
      B[index,][,s] <- 1
    } else {
      B[index,][,s] <- rbinom(n, 1, alpha[s]/sum(alpha[s:k]))
    }
    # Simulated component
    X[B[,s] == 1, s] <- rnorm(sum(B[, s]), mean = mean[s], sd = sd[s])
    # Update number of remaining elements
    n <- n - sum(B[,s])
    # Update the remaining indexes
    index <- index[!(index %in% which(!is.na(X[,s])) )]
    # Substitue NA values with 0
    X[,s] <- ifelse(is.na(X[,s]), 0, X[,s])
  }
  colnames(X) <- paste0("X", 1:k)
  colnames(B) <- paste0("B", 1:k)
  sim <- dplyr::bind_cols(t = 1:nrow(X), X, B, X = rowSums(X))
  return(sim)
}
