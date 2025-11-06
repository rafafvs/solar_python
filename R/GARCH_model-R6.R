#' Implementation of `rugarch` methods for a GARCH(p,q) as R6 class
#'
#' @rdname GARCH_modelR6
#' @name GARCH_modelR6
#' @keywords GARCH
#' @note Version 1.0.0
#' @export
GARCH_modelR6 <- R6::R6Class("GARCH_modelR6",
                             public = list(
                               #' @description
                               #' Initialize a GARCH model with `rugarch` specification
                               #' @param spec GARCH specification from `ugarchspec`.
                               #' @param x Numeric, vector. Time series to be fitted.
                               #' @param weights Numeric, vector. Optional custom weights.
                               #' @param sigma20 Numeric scalar. Target unconditional variance.
                               initialize = function(spec, x, weights, sigma20){
                                 # Store the specification
                                 private[["..spec"]] <- spec
                                 # Time series
                                 private[["..x"]] <- x
                                 # Flexible weights (useful if model do not converge)
                                 if (missing(weights)) {
                                   private[["..w"]] <- rep(1, length(x))
                                 } else {
                                   private[["..w"]] <- ifelse(weights == 0, 0, 1)
                                 }
                                 # Unconditional variance constraint
                                 private[["..sigma20"]] <- ifelse(missing(sigma20), NA_integer_, sigma20)
                               },
                               #' @description
                               #' Fit the GARCH model with `rugarch` function.
                               fit = function(){
                                 # Specification
                                 spec <- private[["..spec"]]
                                 # Safe GARCH fit
                                 safe_GARCH <- purrr::safely(rugarch::ugarchfit)
                                 # Fitted model
                                 model <- safe_GARCH(data = private[["..x"]], spec = spec)$result
                                 if (model@fit$convergence == 0) {
                                   # cli::cli_alert_success("GARCH routine converged!")
                                   params <- model@fit$ipars[,1][model@fit$ipars[,3]==1]
                                   # Extract ARCH parameters
                                   params_arch <- params[stringr::str_detect(names(params), "alpha")]
                                   if (!purrr::is_empty(params_arch)){
                                     private[["..alpha"]] <- params_arch
                                   }
                                   # Extract GARCH parameters
                                   params_garch <- params[stringr::str_detect(names(params), "beta")]
                                   if (!purrr::is_empty(params_garch)){
                                     private[["..beta"]] <- params_garch
                                   }
                                   # Extract intercept
                                   private[["..omega"]] <- params[names(params) == "omega"]
                                   # Store the log-likelihood
                                   private[["..loglik"]] <- -sum(model@fit$log.likelihoods)
                                   # Store the hessian matrix
                                   private[["..hessian"]] <- model@fit$hessian
                                   names_coefs <- names(self$coefficients)
                                   names_coefs <- names_coefs[self$coefficients!=0]
                                   if (!is.na(private$..sigma20)){
                                     names_coefs <- names_coefs[-1]
                                   }
                                   colnames(private[["..hessian"]]) <- names_coefs
                                 } else {
                                   cli::cli_alert_danger("GARCH routine do not converged!")
                                   private$safe_fit(x = private$..x, weights = private$..w)
                                 }
                                 self$update_std.errors()
                               },
                               #' @description
                               #' Filter method from `rugarch` package to compute GARCH variance, residuals and log-likelihoods.
                               #' @param x Numeric, vector. Time series to be filtered.
                               #' @param coefficients Numeric, named vector. Model's coefficients. When missing will be used the fitted parameters.
                               #' @param ... Other arguments passed to `ugarchfilter` function.
                               filter = function(x, coefficients, ...){
                                 # Default time series
                                 if (missing(x)) {
                                   x <- private[["..x"]]
                                 }
                                 # GARCH specification
                                 spec <- self$spec
                                 # Default parameters
                                 if (missing(coefficients)) {
                                   coefficients <- self$coefficients
                                 }
                                 # Remove zero coefficients
                                 coefficients <- coefficients[coefficients != 0]
                                 # Remove intercept if unconditional variance is constrained
                                 if (!is.na(private$..sigma20)) {
                                   coefficients <- coefficients[!stringr::str_detect(names(coefficients), "omega")]
                                 }
                                 # Remove arch component if not present
                                 if (self$order[1] == 0){
                                   coefficients <- coefficients[!stringr::str_detect(names(coefficients), "alpha")]
                                 }
                                 # Remove garch component if not present
                                 if (self$order[2] == 0){
                                   coefficients <- coefficients[!stringr::str_detect(names(coefficients), "beta")]
                                 }
                                 # Update parameters in the specification
                                 spec@model$fixed.pars <- coefficients
                                 # GARCH filter
                                 rugarch::ugarchfilter(spec = spec, data = x)@filter
                               },
                               #' @description
                               #' Log-likelihoods function
                               #' @param coefficients Numeric, named vector. Model's coefficients. When missing will be used the fitted parameters.
                               #' @param x Numeric, vector. Time series used to compute log-likelihoods.
                               #' @param weights Numeric, vector. Optional custom weights.
                               #' @param update Logical. When true the internal log-likelihood will be updated.
                               #' @param ... Other arguments passed to `ugarchfilter` function.
                               logLik = function(coefficients, x, weights, update = FALSE, ...){
                                 # Default values when missing
                                 if (missing(weights) | missing(x)){
                                   x <- private[["..x"]]
                                   weights <- private[["..w"]]
                                 } else if (missing(weights) & !missing(x)){
                                   weights <- rep(1, length(x))
                                 }
                                 # Log-likelihood
                                 log.likelihoods <- self$filter(x, coefficients, ...)$log.likelihoods
                                 log.likelihoods <- sum(log.likelihoods * weights, na.rm = TRUE)
                                 # Update log-likelihood
                                 if (update) {
                                   private[["..loglik"]] <- log.likelihoods
                                 }
                                 return(log.likelihoods)
                               },
                               #' @description
                               #' Update the coefficients of the model
                               #' @param coefficients Numeric, named vector. Model's coefficients.
                               update = function(coefficients){
                                 # Update arch parameters
                                 params_arch <- c()
                                 if (self$order[1] > 0){
                                   params_arch <- coefficients[stringr::str_detect(names(coefficients), "alpha")]
                                   private[["..alpha"]] <- params_arch
                                 }
                                 params_garch <- c()
                                 # Update garch parameters
                                 if (self$order[2] > 0){
                                   params_garch <- coefficients[stringr::str_detect(names(coefficients), "beta")]
                                   private[["..beta"]] <- params_garch
                                 }
                                 # Parameters
                                 params <- c(params_arch, params_garch)
                                 # Unconditional variance constraint
                                 sigma20 <- private$..sigma20
                                 if (!is.na(sigma20)) {
                                   private[["..omega"]] <- c(omega = sigma20 * (1-sum(params)))
                                 } else {
                                   private[["..omega"]] <- coefficients["omega"]
                                   params <- c(private[["..omega"]], params)
                                 }
                                 # Set std. errors to NA
                                 private$..std.errors <- rep(NA, length(self$coefficients))
                               },
                               #' @description
                               #' Numerical computation of the Hessian matrix.
                               #' @param coefficients Numeric, named vector. Model's coefficients. When missing will be used the fitted parameters.
                               #' @param logLik Function, log-likelihood function depending on the parameters and x.
                               #' @param ... Other arguments passed to `logLik` function.
                               update_hessian = function(coefficients, logLik, ...){
                                 # Default coefficients
                                 if (missing(coefficients)) {
                                   coefficients <- self$coefficients
                                 }
                                 # Remove coefficients exactly zero
                                 coefficients <- coefficients[coefficients != 0]
                                 # Default log-likelihood
                                 if (missing(logLik)) {
                                   logLik <- self$logLik
                                 }
                                 # Numerical computation of the Hessian matrix
                                 if (!is.na(private$..sigma20)){
                                   coefficients <- coefficients[!stringr::str_detect(names(coefficients), "omega")]
                                 }
                                 private$..hessian <- numDeriv::hessian(func = logLik, x = coefficients, ...)
                                 colnames(private$..hessian) <- names(coefficients)
                               },
                               #' @description
                               #' Numerical computation of the std. errors of the parameters.
                               #' @param std.errors Numeric std. errors.
                               update_std.errors = function(std.errors){
                                 # Update the standard errors only if hessian is not NA
                                 if (missing(std.errors) & sum(is.na(private$..hessian))==0) {
                                   # Compute std. errors by inverting the Hessian
                                   H_inv <- diag(solve(private$..hessian))
                                   std.errors <- sqrt(H_inv)
                                   names(std.errors) <- colnames(private$..hessian)
                                 } else {
                                   private$..hessian[,] <- NA
                                 }
                                 if (!missing(std.errors)) {
                                   # Extract the coefficients to match the names
                                   std.errors_updated <- self$coefficients
                                   # Remove std. errors that does not match the params names
                                   std.errors <- std.errors[(names(std.errors) %in% names(std.errors_updated))]
                                   # Update std. errors that are not included with NA
                                   std.errors_updated[!(names(std.errors_updated) %in% names(std.errors))] <- NA_integer_
                                   # Updated std. errors included
                                   std.errors_updated[names(std.errors_updated) %in% names(std.errors)] <- std.errors
                                   #  ================ Private ================
                                   private$..std.errors <- std.errors_updated
                                 }
                               },
                               #' @description
                               #' Next step GARCH std. deviation forecast
                               #' @param x Numeric, vector. Past residuals.
                               #' @param sigma Numeric, vector. Past garch std. deviations.
                               #' @param n.ahead Numeric, scalar. Number of steps ahead.
                               next_step = function(x = 1, sigma = 1, n.ahead = 1){
                                 x_t <- x
                                 sigma_t <- sigma
                                 prev_sigma_t <- sigma_t[1]*10
                                 nstep <- 0
                                 while(prev_sigma_t - sigma_t[1] != 0 & nstep < n.ahead){
                                   # Initialize GARCH variance
                                   sigma2_t_plus_1 <- self$omega
                                   # ARCH(p) component
                                   if (self$order[1] > 0) {
                                     sigma2_t_plus_1 <- sigma2_t_plus_1 + sum(self$alpha * x_t^2)
                                   }
                                   # GARCH(q) component
                                   if (self$order[2] > 0) {
                                     sigma2_t_plus_1 <- sigma2_t_plus_1 + sum(self$beta * sigma_t^2)
                                   }
                                   # Store previous std. dev
                                   prev_sigma_t <- sigma_t[1]
                                   # Update residuals
                                   x_t <- c(sqrt(sigma2_t_plus_1), lag(x_t, 1)[-1])
                                   # Update std. dev
                                   sigma_t <- c(sqrt(sigma2_t_plus_1), lag(sigma_t, 1)[-1])
                                   # Update number of step
                                   nstep <- nstep + 1
                                 }
                                 return(sigma_t[1])
                               },
                               #' @description
                               #' Print method for `GARCH_modelR6` class.
                               print = function(){
                                 # Format parameters with std. errors
                                 green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                                 red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                                 msg_col <- function(x) ifelse(x, green(x), red(x))
                                 # Intercept
                                 par.intercept <- format(self$omega, digits = 4)
                                 std.error.intercept <- private$..std.errors[names(private$..std.errors) == names(par.intercept)]
                                 format.intercept <- paste0(par.intercept, " (", format(std.error.intercept, digits = 3), ")")
                                 # ARCH parameters
                                 par.arch <- format(self$alpha, digits = 4)
                                 std.error.arch <- private$..std.errors[names(private$..std.errors) %in% names(par.arch)]
                                 format.arch <- paste0(par.arch, " (", format(std.error.arch, digits = 3), ")")
                                 # GARCH parameters
                                 par.garch <- format(self$beta, digits = 4)
                                 std.error.garch <- private$..std.errors[names(private$..std.errors) %in% names(par.garch)]
                                 format.garch <- paste0(par.garch, " (", format(std.error.garch, digits = 3), ")")
                                 # Model name
                                 model_name <- paste0("GARCH", " (", self$order[1], ", ", self$order[2], ")")
                                 # ********************************************************
                                 msg_0 <- paste0("--------------------- ", model_name, "--------------------- \n")
                                 msg_1 <- paste0("ARCH: ", msg_col(self$order[1] != 0), "\n")
                                 msg_2 <- paste0("GARCH: ", msg_col(self$order[2] != 0), "\n")
                                 msg_3 <- paste0("Version: ", private$version, "\n")
                                 msg_4 <- paste0("Unconditional Variance: ", private$..sigma20, "\n")
                                 msg_5 <- "------------------------------------------------\n"
                                 msg_6 <- paste0("Intercept: ", format.intercept, " \n")
                                 msg_7 <- paste0("Alpha: ", format.arch, " \n")
                                 msg_8 <- paste0("Beta: ", format.garch, " \n")
                                 msg_9 <- paste0("Log-lik: ", format(self$loglik, digits = 3), "\n")
                                 cat(c(msg_0, msg_1, msg_2, msg_3, msg_4,
                                       msg_5, msg_6, msg_7, msg_8, msg_9))
                               }
                             ),
                             private = list(
                               version = "1.0.0",
                               ..x = NA,
                               ..w = NA,
                               ..sigma20 = NA,
                               ..spec = NA,
                               ..omega = c(omega = 1),
                               ..alpha = c(alpha1 = 0),
                               ..beta = c(beta1 = 0),
                               ..loglik = NA,
                               ..hessian = matrix(NA, 3, 3),
                               ..std.errors = c(NA, NA, NA),
                               #' Manual fit of the GARCH model
                               safe_fit = function(params, x, weights){
                                 # GARCH specification
                                 spec <- private[["..spec"]]
                                 # Address failed convergence with initial random parameters
                                 if (missing(params)) {
                                   # Detect parameters to be estimated
                                   fixed.pars <- spec@model$pars[,4][spec@model$pars[,4] == 1]
                                   fixed.pars[1:length(fixed.pars)] <- c(runif(length(fixed.pars), 0, 0.2))
                                   params <- fixed.pars
                                 } else {
                                   # Ensure omega is removed
                                   if (!is.na(private[["..sigma20"]])) {
                                     params <- params[which(names(params) != "omega")]
                                   }
                                 }
                                 # Define a custom loss function
                                 loss_function = function(params, x, weights){
                                   # Extract arch parameters
                                   index_arch <- which(stringr::str_detect(names(params), "alpha"))
                                   alpha <- c(alpha1 = 0)
                                   if (!purrr::is_empty(index_arch)) {
                                     alpha <- params[index_arch]
                                   }
                                   # Extract garch parameters
                                   index_garch <- which(stringr::str_detect(names(params), "beta"))
                                   beta <- c(beta1 = 0)
                                   if (!purrr::is_empty(index_garch)) {
                                     beta <- params[index_garch]
                                   }
                                   # Check stationarity condition
                                   stationary_condition <- sum(beta < 0) == 0 & sum(alpha < 0) == 0 & sum(alpha + beta) < 1
                                   if (!stationary_condition) {return(NA)}
                                   # Compute log likelihoods
                                   loglik <- self$filter(x, params)$log.likelihoods
                                   sum(-loglik * weights)
                                 }
                                 # Optimization function
                                 opt <- optim(params, self$logLik, x = x, weights = weights)
                                 # Extract the parameters
                                 params <- opt$par
                                 # Extract ARCH parameters
                                 private[["..alpha"]] <- params[stringr::str_detect(names(params), "alpha")]
                                 # Extract GARCH parameters
                                 private[["..beta"]] <- params[stringr::str_detect(names(params), "beta")]
                                 # Update log-likelihood
                                 private[["..loglik"]] <- opt$value
                                 # Update the Hessian
                                 self$update_hessian()
                               }
                             ),
                             active = list(
                               #' @field spec A `rugarch` object containing the model's specification.
                               spec = function(){
                                 private$..spec
                               },
                               #' @field omega Numeric scalar, intercept parameter.
                               omega = function(){
                                 private$..omega
                               },
                               #' @field alpha Numeric vector, ARCH parameters.
                               alpha = function(){
                                 private$..alpha
                               },
                               #' @field beta Numeric vector, GARCH parameters.
                               beta = function(){
                                 private$..beta
                               },
                               #' @field order Numeric scalar, model's order.
                               order = function(){
                                 c(ARCH = length(self$alpha[self$alpha != 0]), GARCH = length(self$beta[self$beta != 0]))
                               },
                               #' @field coefficients Numeric vector, model's coefficients.
                               coefficients = function(){
                                 c(self$omega, self$alpha, self$beta)
                               },
                               #' @field vol Numeric scalar, long-term unconditional std. deviation.
                               vol = function(){
                                 sqrt(self$omega / (1 - sum(self$alpha+self$beta)))
                               },
                               #' @field loglik model log-likelihood
                               loglik = function(){
                                 private$..loglik
                               },
                               #' @field tidy Tibble with estimated parameters and relative std. errors.
                               tidy = function(){
                                 dplyr::tibble(
                                   term = names(self$coefficients),
                                   estimate = self$coefficients,
                                   std.error = private$..std.errors,
                                 )
                               }
                             )
                            )
