#' Monthly Gaussian mixture with two components
#'
#' @examples
#' # Model fit
#' model <- solarModel$new(spec)
#' model$fit()
#' # Inputs
#' x <- model$data$u_tilde
#' w <- model$data$weights
#' date <- model$data$date
#' # Solar Mixture object
#' sm <- solarMixture$new()
#' sm$fit(x, date, w)
#' params <- sm$parameters
#' sm$std.errors
#' # params[1,]$mu1 <- params[1,]$mu1*0.9
#' # sm$update(means = params[,c(2,3)])
#' @rdname solarMixture
#' @name solarMixture
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
solarMixture <-  R6::R6Class("solarMixture",
                             public = list(
                               #' @field maxit Integer, maximum number of iterations.
                               maxit = 100,
                               #' @field maxrestarts Integer, maximum number of restarts.
                               maxrestarts = 500,
                               #' @field abstol Numeric, absolute level for convergence.
                               abstol = 10e-5,
                               #' @field components Integer, number of components.
                               components = 2,
                               #' @field mu1 Function, see \code{\link{monthlyParams}}.
                               mu1 = NA,
                               #' @field mu2 Function, see \code{\link{monthlyParams}}.
                               mu2 = NA,
                               #' @field sd1 Function, see \code{\link{monthlyParams}}.
                               sd1 = NA,
                               #' @field sd2 Function, see \code{\link{monthlyParams}}.
                               sd2 = NA,
                               #' @field prob Function, see \code{\link{monthlyParams}}.
                               prob = NA,
                               #' @description
                               #' Initialize a `solarMixture` object
                               #' @param components (`integer(1)`), number of components.
                               #' @param maxit (`integer(1)`) Numeric, maximum number of iterations.
                               #' @param abstol (`numeric(1)`) Numeric, absolute level for convergence.
                               #' @param maxrestarts (`integer(1)`) Numeric, maximum number of restarts.
                               initialize = function(components = 2, maxit = 5000, maxrestarts = 500, abstol = 1e-8){
                                 self$components <- components
                                 self$maxit <- maxit
                                 self$maxrestarts <- maxrestarts
                                 self$abstol <- abstol
                               },
                               #' @description
                               #' Fit the parameters with mclust package
                               #' @param x vector
                               #' @param date date vector
                               #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                               #' When `missing` all the available observations will be used.
                               #' @param method Character, package used to fit the parameters. Can be `mclust` or `mixtools`.
                               #' @param mu_target Numeric vector with length 12, target mean of the mixture to match.
                               #' @param var_target Numeric vector with length 12, target variance of the mixture to match.
                               fit = function(x, date, weights, method = "mixtools", mu_target = rep(NA, 12), var_target = rep(NA, 12)){
                                 # Add data
                                 private$x <- x
                                 private$date <- date
                                 # Weights
                                 if (missing(weights)){
                                   w <- rep(1, n)
                                 } else {
                                   w <- ifelse(weights == 0, 0, 1)
                                 }
                                 private$w <- w
                                 data <- self$data
                                 # Gaussian Mixture parameters
                                 for(m in 1:12){
                                   data_months <- dplyr::filter(data, Month == m)
                                   w <- data_months$w
                                   # Monthly data
                                   eps <- data_months$x
                                   private$date_month[[m]] <- data_months$date
                                   # Fitted model
                                   GM_model <- gaussianMixture$new(maxit = self$maxit, maxrestarts = self$maxrestarts,
                                                                   abstol = self$abstol, components = self$components)
                                   GM_model$fit(eps, w, method = method, mu_target = mu_target[m], var_target = var_target[m])
                                   # Add model
                                   private$..model[[m]] <- GM_model$clone(deep = TRUE)
                                 }
                                 names(private$..model) <- lubridate::month(1:12, label = TRUE)
                                 names(private$date_month) <- lubridate::month(1:12, label = TRUE)
                                 # Initialize monthly function for Mixture parameters
                                 self$mu1 <- monthlyParams$new(self$coefficients$mu1)
                                 self$mu2 <- monthlyParams$new(self$coefficients$mu2)
                                 self$sd1 <- monthlyParams$new(self$coefficients$sd1)
                                 self$sd2 <- monthlyParams$new(self$coefficients$sd2)
                                 self$prob <- monthlyParams$new(self$coefficients$p1)
                               },
                               #' @description
                               #' Update means, sd, p .
                               #' @param means Numeric matrix of means parameters.
                               #' @param sd Numeric matrix of std. deviation parameters.
                               #' @param p Numeric matrix of probability parameters.
                               update = function(means, sd, p){
                                 # Mean parameters
                                 if (missing(means)) {
                                   means <- self$means
                                 }
                                 # Std. deviations parameters
                                 if (missing(sd)) {
                                   sd <- self$sd
                                 }
                                 # Probability parameters
                                 if (missing(p)) {
                                   p <- self$p
                                 }
                                 M <- length(private$..model)
                                 # Update parameters
                                 for(m in 1:M){
                                   private$..model[[m]]$update(means = means[m,], sd = sd[m,], p = p[m,])
                                 }
                                 # Update monthly function for Mixture parameters
                                 self$mu1 <- monthlyParams$new(self$coefficients$mu1)
                                 self$mu2 <- monthlyParams$new(self$coefficients$mu2)
                                 self$sd1 <- monthlyParams$new(self$coefficients$sd1)
                                 self$sd2 <- monthlyParams$new(self$coefficients$sd2)
                                 self$prob <- monthlyParams$new(self$coefficients$p1)
                               },
                               #' @description
                               #' Apply the `$update_logLik()` method to all the `gaussianMixture` models
                               update_logLik = function(){
                                 M <- length(private$..model)
                                 # Update log-likelihood
                                 for(m in 1:M){
                                   private$..model[[m]]$update_logLik()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Apply the `$update_empiric_parameters()` method to all the `gaussianMixture` models
                               update_empiric_parameters = function(){
                                 M <- length(private$..model)
                                 # Update empiric parameters
                                 for(m in 1:M){
                                   private$..model[[m]]$update_empiric_parameters()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Apply the `$filter()` method to all the `gaussianMixture` models
                               filter = function(){
                                 M <- length(private$..model)
                                 # Update parameters, fitted data and log-likelihood
                                 for(m in 1:M){
                                   private$..model[[m]]$filter()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Apply the `$Hessian()` method to all the `gaussianMixture` models
                               Hessian = function(){
                                 M <- length(private$..model)
                                 # Update the Hessian and std. errors
                                 for(m in 1:M){
                                   private$..model[[m]]$Hessian()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Substitute the empiric parameters with EM parameters. If evaluated again
                               #' the EM parameters will be substituted back.
                               use_empiric_parameters = function(){
                                 M <- length(private$..model)
                                 for(i in 1:M){
                                   private$..model[[i]]$use_empiric_parameters()
                                 }
                                 private$..use_empiric <- private$..model[[1]]$use_empiric
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Log-likelihood
                               #' @param x vector
                               #' @param date dates
                               logLik = function(x, date){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date), x = x)
                                 log_lik <- 0
                                 for(m in 1:M){
                                   data_month <- dplyr::filter(df, Month == m)
                                   log_lik <- log_lik + self$model[[m]]$logLik(filter(df, Month == m)$x)
                                 }
                                 return(log_lik)
                               },
                               #' @description
                               #' Compute the grades
                               #' @param x vector
                               #' @param date dates
                               grades = function(x, date){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date), x = x, u = NA)
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly cdf
                                   cdf_M <- function(x) pmixnorm(x, means, sd, p)
                                   # Grades
                                   df$u[idx_month] <- cdf_M(x[idx_month])
                                 }
                                 return(df$u)
                               },
                               #' @description
                               #' Compute the VaR with certain confidence levels
                               #' @param x vector
                               #' @param date dates
                               #' @param alpha confidence levels for the VaR
                               VaR = function(date, alpha = 0.05){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date))
                                 # Compute VaR
                                 VaR_alpha <- matrix(0, nrow = nrow(df), ncol = length(alpha))
                                 # Assign standard names
                                 colnames(VaR_alpha) <- paste0("VaR_", alpha)
                                 alpha <- alpha[order(alpha)]
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly quantile
                                   Q_M <- function(x) qmixnorm(x, means, sd, p)
                                   # VaR
                                   VaR_alpha_M <- Q_M(alpha)
                                   for(j in 1:length(idx_month)){
                                     VaR_alpha[idx_month[j],] <- VaR_alpha_M
                                   }
                                 }
                                 return(VaR_alpha)
                               },
                               #' @description
                               #' Compute the VaR with certain confidence levels
                               #' @param x vector
                               #' @param date dates
                               #' @param alpha confidence levels for the VaR
                               ES = function(date, alpha = 0.05){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date))
                                 # Compute VaR
                                 ES_alpha <- matrix(0, nrow = nrow(df), ncol = length(alpha))
                                 # Assign standard names
                                 colnames(ES_alpha) <- paste0("ES_", alpha)
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly quantile
                                   Q_M <- function(x) qmixnorm(x, means, sd, p)
                                   VaR_alpha <- Q_M(alpha)
                                   # Expected shortfalls
                                   # Density function
                                   f_R <- function(x) dmixnorm(x, means, sd, p)
                                   # Expected Shortfall
                                   ES <- function(VaR, alpha) integrate(function(x) x*f_R(x), lower = -Inf,
                                                                        upper = VaR, subdivisions = 50L, stop.on.error = FALSE)$value/alpha
                                   ES_alpha_M <- purrr::map2_dbl(VaR_alpha, alpha, ~ES(.x, .y))
                                   for(j in 1:length(idx_month)){
                                     ES_alpha[idx_month[j],] <- ES_alpha_M
                                   }
                                 }
                                 return(ES_alpha)
                               },
                               #' @description
                               #' Print method for `solarMixture` class.
                               print = function(){
                                 M <- length(private$..model)
                                 for(m in 1:M){
                                   nmonth <- as.character(lubridate::month(m, label = TRUE, abbr = FALSE))
                                   self$model[[m]]$print(nmonth)
                                 }
                               }
                             ),
                             private = list(
                               version = "1.0.0",
                               x = NA,
                               w = NA,
                               date = NA,
                               date_month = list(),
                               ..model = list(),
                               ..use_empiric = FALSE,
                               deep_clone = function(name, value){
                                 if (name == "..model") {
                                   # Clonazione profonda degli oggetti kernelRegression all'interno della lista ..models
                                   return(lapply(value, function(model) {
                                     if (!is.null(model)) {
                                       model$clone(deep = TRUE)  # Clonazione profonda per ogni oggetto kernelRegression
                                     } else {
                                       NULL
                                     }
                                   }))
                                 } else {
                                   # Per altri campi, usa il comportamento di clonazione predefinito
                                   return(value)
                                 }
                               }
                             ),
                             active = list(
                               #' @field data A tibble with the following columns:
                               #' \describe{
                               #'  \item{date}{Time series of dates.}
                               #'  \item{Month}{Vector of Month.}
                               #'  \item{x}{Time series used for fitting.}
                               #'  \item{w}{Time series of weights.}}
                               data = function(){
                                 dplyr::tibble(date = private$date, Month = lubridate::month(date), x = private$x, w = private$w)
                               },
                               #' @field means Matrix of means where a row represents a month and a column a mixture component.
                               means = function(){
                                 M <- length(private$..model)
                                 means <- matrix(0, nrow = M, ncol = self$components)
                                 for(m in 1:M){
                                   means[m,] <- unlist(self$model[[m]]$means)
                                 }
                                 colnames(means) <- names(self$model[[1]]$means)
                                 rownames(means) <- lubridate::month(1:M, label = TRUE)
                                 return(means)
                               },
                               #' @field sd Matrix of std. deviations where a row represents a month and a column a mixture component.
                               sd = function(){
                                 M <- length(private$..model)
                                 sd <- matrix(0, nrow = M, ncol = self$components)
                                 for(m in 1:M){
                                   sd[m,] <- unlist(self$model[[m]]$sd)
                                 }
                                 colnames(sd) <- names(self$model[[1]]$sd)
                                 rownames(sd) <- lubridate::month(1:M, label = TRUE)
                                 return(sd)
                               },
                               #' @field p Matrix of probabilities where a row represents a month and a column a mixture component.
                               p = function(){
                                 M <- length(private$..model)
                                 p <- matrix(0, nrow = M, ncol = self$components)
                                 for(m in 1:M){
                                   p[m,] <- unlist(self$model[[m]]$p)
                                 }
                                 colnames(p) <- names(self$model[[1]]$p)
                                 rownames(p) <- lubridate::month(1:M, label = TRUE)
                                 return(p)
                               },
                               #' @field model Named List with 12 \code{\link{gaussianMixture}} objects.
                               model = function(){
                                 private$..model
                               },
                               #' @field use_empiric logical to denote if empiric parameters are currently used
                               use_empiric = function(){
                                 private$..use_empiric
                               },
                               #' @field loglik Numeric, total log-likelihood.
                               loglik = function(){
                                 sum(purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field fitted A `tibble` with the classified series
                               fitted = function(){
                                 df_fitted <- purrr::map2_df(self$model, private$date_month, ~ dplyr::bind_cols(date = .y, dplyr::select(.x$fitted, B = "B1", uncertanty)))
                                 dplyr::arrange(df_fitted, date)
                               },
                               #' @field moments A `tibble` with the theoric moments. It contains:
                               #' \describe{
                               #'  \item{Month}{Month of the year.}
                               #'  \item{mean}{Theoric monthly expected value of the mixture model.}
                               #'  \item{variance}{Theoric monthly variance of the mixture model.}
                               #'  \item{skewness}{Theoric monthly skewness.}
                               #'  \item{kurtosis}{Theoric monthly kurtosis.}
                               #'  \item{nobs}{Number of observations used for fitting.}
                               #'  \item{loglik}{Monthly log-likelihood.}}
                               moments = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~.x$moments), loglik = purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field coefficients A `tibble` with the fitted parameters.
                               coefficients = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~.x$model))
                               },
                               #' @field std.errors A `tibble` with the fitted std.errors
                               std.errors = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~dplyr::bind_cols(purrr::map(.x$std.errors, ~dplyr::bind_rows(as.list(.x))))))
                               },
                               #' @field summary A `tibble` with the fitted std.errors
                               summary = function(){
                                 M <- length(private$..model)
                                 purrr::map2_df(1:M, self$model, ~dplyr::bind_cols(Month = .x, .y$summary))
                               }
                             )
                            )

#' Correct the moments to ensure moments matching
#'
#' @keywords solarMixture internal
#' @export
#' @noRd
solarMixture_moments_match <- function(coefficients, e_target, v_target, sk_target, kt_target, x){
  # Loss function
  loss_function <- function(params, prob, e_target, v_target, sk_target, kt_target, x){
    # Extract the parameters
    means <- params[stringr::str_detect(names(params), "mu")]
    sigma <- params[stringr::str_detect(names(params), "sd")]
    probs <- c(p1 = prob[[1]], p2 = 1 - prob[[1]])
    # Compute the moments
    mom <- GM_moments(means, sigma, probs)
    loss <- 0
    # Compute the loss from target moments
    if (!is.na(e_target)){
      loss <- loss + abs(mom$mean - e_target)
    }
    if (!is.na(v_target)){
      loss <- loss + abs(mom$variance - v_target)
    }
    if (!is.na(sk_target)){
      loss <- loss + abs(mom$skewness - sk_target)
    }
    if (!is.na(kt_target)){
      loss <- loss + abs(mom$kurtosis - kt_target)
    }
    # Compute log-likelihood
    loss <- 10000*loss^2 #- GM_loglik(means, sigma, probs, x)
    return(loss)
  }
  # Monthly optimization
  opt_coefficients <- coefficients
  for(nmonth in 1:nrow(coefficients)) {
    # Initial parameters
    params <- unlist(purrr::flatten(coefficients[nmonth,]))
    # Optimal parameters
    opt <- optim(par = params[-5], fn = loss_function, x = x[[nmonth]], prob = params[5],
                 e_target = e_target[nmonth], v_target = v_target[nmonth], sk_target = sk_target[nmonth],
                 kt_target = kt_target[nmonth])
    opt_coefficients[nmonth, 1:4] <- dplyr::bind_rows(opt$par)
  }
  opt_coefficients$p2 <- 1 - opt_coefficients$p1
  return(opt_coefficients)
}


# Monthly multivariate Gaussian mixture with two components
solarModel_mvmixture <- function(model_Ct, model_GHI){

  # Extract a bivariate dataset
  data_GHI <- dplyr::select(model_GHI$data, date, Month, GHI = "u_tilde")
  data_Ct <- dplyr::select(model_Ct$data, date, Month, Ct = "u_tilde")
  data <- dplyr::inner_join(data_GHI, data_Ct, by = c("date", "Month"))
  # Remove outliers
  outliers_date <- c(model_Ct$.__enclos_env__$private$outliers$date, model_GHI$.__enclos_env__$private$outliers$date)
  data <- dplyr::filter(data, !(date %in% outliers_date))

  # Extract Gaussian mixture parameters
  NM_model_GHI <- model_GHI$.__enclos_env__$private$..NM_model
  NM_model_GHI$rho_up <- NM_model_GHI$rho_dw <- 0
  NM_model_Ct <- model_Ct$.__enclos_env__$private$..NM_model
  NM_model_Ct$rho_up <- NM_model_Ct$rho_dw <- 0

  # Gaussian Mixture parameters
  m <- 1
  # Monthly data
  data_months <- list()
  for(m in 1:12){
    data_months[[m]] <- dplyr::filter(data, Month == m)
    # Monthly data
    eps <- data_months[[m]][,c(3,4)]
    # Fitted model
    gm <- mvgaussianMixture(eps, components = 2, na.rm = TRUE)
    # Update Gaussian mixture parameters (GHI)
    NM_model_GHI$mu_up[m] <- gm$params$means[1,2]
    NM_model_GHI$mu_dw[m] <- gm$params$means[2,2]
    NM_model_GHI$sd_up[m] <- sqrt(gm$params$sigma2[1,2])
    NM_model_GHI$sd_dw[m] <- sqrt(gm$params$sigma2[2,2])
    NM_model_GHI$p_up[m] <- gm$params$p[1]
    NM_model_GHI$p_dw[m] <- gm$params$p[2]
    NM_model_GHI$rho_up[m] <- gm$params$rho[1]
    NM_model_GHI$rho_dw[m] <- gm$params$rho[2]
    # Update Gaussian mixture parameters (Ct)
    NM_model_Ct$mu_up[m] <- gm$params$means[1,1]
    NM_model_Ct$mu_dw[m] <- gm$params$means[2,1]
    NM_model_Ct$sd_up[m] <- sqrt(gm$params$sigma2[1,1])
    NM_model_Ct$sd_dw[m] <- sqrt(gm$params$sigma2[2,1])
    NM_model_Ct$p_up[m] <- gm$params$p[1]
    NM_model_Ct$p_dw[m] <- gm$params$p[2]
    NM_model_Ct$rho_up[m] <- gm$params$rho[1]
    NM_model_Ct$rho_dw[m] <- gm$params$rho[2]
    # Add fitted Bernoulli series
    data_months[[m]]$B <- gm$B_hat$B1
  }
  # Fitted series of bernoulli
  data_months <- dplyr::select(dplyr::bind_rows(data_months), date, B)

  # Update data (Ct)
  data_Ct <- model_Ct$.__enclos_env__$private$..data
  data_Ct <- dplyr::left_join(dplyr::select(data_Ct, -B), data_months, by = "date")
  data_Ct$B[is.na(data_Ct$B)] <- 0
  model_Ct$.__enclos_env__$private$..data <- data_Ct
  model_Ct$.__enclos_env__$private$..NM_model <- NM_model_Ct
  # Update data (GHI)
  data_GHI <- model_GHI$.__enclos_env__$private$..data
  data_GHI <- dplyr::left_join(dplyr::select(data_GHI, -B), data_months, by = "date")
  data_GHI$B[is.na(data_GHI$B)] <- 0
  model_GHI$.__enclos_env__$private$..data <- data_GHI
  model_GHI$.__enclos_env__$private$..NM_model <- NM_model_GHI

  structure(
    list(
      model_Ct = model_Ct,
      model_GHI = model_GHI
    ),
    class = c("solarModelMixture", "list")
  )
}





