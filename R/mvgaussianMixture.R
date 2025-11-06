#' Multivariate gaussian mixture
#'
#' @rdname mvgaussianMixture
#' @name mvgaussianMixture
#' @export
mvgaussianMixture <- function(x, means, sd, p, components = 2, maxit = 100, abstol = 10e-15, na.rm = FALSE){

  # Ensure that there are not NAs or NaN observations
  idx_NA <- is.na(x)
  if (any(idx_NA)) {
    x <- na.omit(x)
    wrn <- paste0("Removed ", sum(idx_NA), " NA observations!")
    warning(wrn)
  }
  # Number of observations
  n_w <- nrow(x)
  # Number of variables
  j_w <- ncol(x)
  # Empirical moments
  e_x_hat <- colMeans(x, na.rm = na.rm)
  v_x_hat <- apply(x, 2, var)

  # Default starting means
  if (missing(means) || any(is.na(means))){
    # Initialize a matrix for the components
    means <- matrix(0, nrow = components, ncol = j_w, dimnames = list(1:components, dimnames(x)[[2]]))
    probs <- seq(0.8, 0.2, length.out = components)
    for(k in 1:components){
      means[k,] <- apply(x, 2, quantile, probs = probs[k])
    }
  }

  # Default std. deviations
  if (missing(sd) || any(is.na(sd))){
    sd <- list()
    for(k in 1:components){
      sd[[k]] <- diag(v_x_hat)
    }
  }
  # Default probabilities
  if (missing(p) || any(is.na(p))) {
    p <- rep(1/components, components)
  }

  # Routine
  # 0. Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  prev_responsibilities <- matrix(0, nrow = n_w, ncol = components)
  previous_params <- list(mean = means, sd = sd, p = p)
  iteration <- 1
  # EM Algorithm
  for (iteration in 1:maxit) {
    # E-step: posterior probabilities
    responsibilities <- prev_responsibilities
    for (i in 1:n_w) {
      for(k in 1:components){
        responsibilities[i, k] <- previous_params$p[k]*mvtnorm::dmvnorm(x[i,], mean = previous_params$mean[k,], sigma = previous_params$sd[[k]])
      }
      # Normalize the posterior probabilities
      responsibilities[i,] <- (responsibilities[i,])/sum(responsibilities[i,], na.rm = TRUE)
      responsibilities[i,][is.na(responsibilities[i,])] <- 0
    }

    # Optimal parameters
    k <- 1
    params <- previous_params
    # M-step: Update the parameters
    for(k in 1:components){
      # Normalizing factor for each component
      n_k <- sum(responsibilities[, k], na.rm = na.rm)
      # Mean parameters k-component
      params$mean[k,] <- apply(responsibilities[, k]*x, 2, sum)/n_k
      # Covariance matrix k-component
      params$sd[[k]] <- diag(apply(x^2*responsibilities[, k], 2, sum)/n_k - params$mean[k,]^2)
      params$sd[[k]][1,2] <- params$sd[[k]][2,1] <- sum(x[,1]*x[,2]*responsibilities[, k])/n_k - params$mean[k,][1]*params$mean[k,][2]
      # Probability k-component
      params$p[k] <- n_k/n_w
    }

    if(any(params$p > 0.9)){
      warning("Probability greater than 0.9 Break!")
      params <- previous_params
      break
    }

    # Calculate the log-likelihood
    log_likelihood <- 0
    for(i in 1:n_w) {
      ll <- 0
      for(k in 1:components){
        ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
      }
      log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
    }

    # Check for convergence
    stop_condition <- abs(log_likelihood - previous_log_likelihood) < abstol
    if (stop_condition) {
      break
    } else {
      # Update log-likelihood
      previous_log_likelihood <- log_likelihood
      # Update parameters
      previous_params <- params
    }
    print(log_likelihood)
    if (iteration == maxit) {
      message("Max iteration reached (", iteration, ")")
    }
  }

  # Final classification of each component
  B_hat <- matrix(0, nrow = n_w, ncol = components)
  for(i in 1:n_w) {
    ll <- c()
    for(k in 1:components){
      ll[k] <- sum(responsibilities[i,k])
    }
    B_hat[i, which.max(ll)] <- 1
  }
  colnames(B_hat) <- paste0("B", 1:components)
  B_hat <- dplyr::as_tibble(B_hat)

  # ML-parameters
  params <- previous_params
  # Reorder the components by decreasing means
  colnames(responsibilities) <- paste0("B", 1:components)
  responsibilities <- dplyr::as_tibble(responsibilities)

  # Log-likelihood on fitted parameters
  # Calculate the log-likelihood
  # Calculate the log-likelihood
  log_likelihood <- 0
  for(i in 1:n_w) {
    ll <- 0
    for(k in 1:components){
      ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
    }
    log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
  }

  upd_params <- list()
  upd_params$means <- params$mean
  upd_params$sigma2 <- params$mean
  upd_params$rho <- c(0,0)
  upd_params$p <- params$p
  i <- 1
  for(i in 1:length(params$sd)){
    upd_params$sigma2[i,] <- diag(params$sd[[i]])
    upd_params$rho[i] <-  params$sd[[i]][upper.tri(params$sd[[i]])]/prod(sqrt(diag(params$sd[[i]])))
  }
  structure(
    list(
      B_hat = B_hat,
      iteration = iteration,
      params = upd_params,
      responsibilities = responsibilities,
      log_lik = log_likelihood
    ),
    class = c("mvgaussianMixture")
  )
}



