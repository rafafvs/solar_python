#' R6 class for ARMA(p, q) model
#'
#' @rdname ARMA_modelR6
#' @name ARMA_modelR6
#' @keywords ARMA
#' @note Version 1.0.1
#' @seealso [stats::arima()] which is wrapped in the method `fit`.
#' @export
ARMA_modelR6 <- R6::R6Class("ARMA_modelR6",
                          public = list(
                            #' @description
                            #' Initialize an ARMA model
                            #' @param arOrder Numeric scalar, order for Autoregressive component.
                            #' @param maOrder Numeric scalar, order for Moving-Average component.
                            #' @param include.intercept Logical. When `TRUE` the intercept will be included. The default is `FALSE`.
                            initialize = function(arOrder = 1, maOrder = 1, include.intercept = FALSE){
                              # Logical value for the intercept parameter
                              private$include.intercept <- include.intercept
                              # Store AR order
                              private$..arOrder <- arOrder
                              # Store MA order
                              private$..maOrder <- maOrder
                              # Pre-compute vector b
                              private$..b <- ARMA_vector_b(arOrder, maOrder)
                            },
                            #' @description
                            #' Fit the ARMA model with `arima` function.
                            #' @param x Numeric vector, time series to fit.
                            fit = function(x){
                              # Fitted model
                              ARMA_model <- arima(x, order = c(self$order[1], 0, self$order[2]),
                                                  include.mean = private$include.intercept, method = "CSS")
                              # Standardize parameters names
                              params <- ARMA_model$coef
                              # Extract the std.errors
                              std.errors <- sqrt(diag(ARMA_model$var.coef))
                              names(std.errors) <- names(params)
                              # Extract intercept
                              intercept <- c(intercept = 0)
                              std.errors_intercept <- c(intercept = NA)
                              if (private$include.intercept) {
                                index <- which(names(params) == "intercept")
                                intercept <- c(intercept = params[[index]])
                                std.errors_intercept <- c(intercept = std.errors[[index]])
                                std.errors <- std.errors[-c(index)]
                                params <- params[-c(index)]
                              }
                              # AR coefficients without intercept
                              phi <- c()
                              std.errors_ar <- c()
                              if (self$order[1] > 0){
                                index <- stringr::str_detect(names(params), "ar")
                                phi <- params[index]
                                std.errors_ar <- std.errors[index]
                                names(phi) <- names(std.errors_ar) <- paste0("phi_", 1:self$order[1])
                              }
                              # MA coefficients without intercept
                              theta <- c()
                              std.errors_ma <- c()
                              if (self$order[2] > 0){
                                index <- stringr::str_detect(names(params), "ma")
                                theta <- params[index]
                                std.errors_ma <- std.errors[index]
                                names(theta) <- names(std.errors_ma) <- paste0("theta_", 1:self$order[2])
                              }
                              # *********************************************
                              # Store coefficients
                              private$..intercept <- intercept
                              private$..phi <- phi
                              private$..theta <- theta
                              # Compute companion matrix
                              private$..A <- ARMA_companion_matrix(phi, theta)
                              # Store the fitted model
                              private$..model <- ARMA_model
                              # Store the std. errors
                              private$..std.errors <- c(std.errors_intercept, std.errors_ar, std.errors_ma)
                              # Update fitted variance
                              private$..sigma2 <- sqrt(ARMA_model$sigma2)
                            },
                            #' @description
                            #' Filter the time-series and compute fitted values and residuals.
                            #' @param x Numeric vector, time series to filter.
                            filter = function(x){
                              ARMA_filter(x, self$A, self$b, self$intercept)
                            },
                            #' @description
                            #' Next step function
                            #' @param x Numeric vector, state vector with past observations and residuals.
                            #' @param n.ahead Numeric scalar, forecasted steps ahead.
                            #' @param eps Numeric vector, optional realized residuals.
                            next_step = function(x, n.ahead = 1, eps = 0){
                              ARMA_next_step(n.ahead, x, self$A, self$b, self$intercept, eps)
                            },
                            #' @description
                            #' Forecast expected value
                            #' @param h Numeric scalar, number of steps ahead.
                            #' @param X0 Numeric vector with length `p + q`, state vector of past values.
                            expectation = function(h = 1, X0){
                              if (missing(X0)){
                                X0 <- rep(0, sum(self$order))
                              }
                              ARMA_expectation(h, X0, self$A, self$b, self$intercept)
                            },
                            #' @description
                            #' Forecast variance
                            #' @param h Numeric scalar, number of steps ahead.
                            #' @param sigma2 Numeric scalar, std. deviation of the residuals.
                            variance = function(h = 1, sigma2 = 1){
                              ARMA_variance(h, self$A, self$b, sigma2)
                            },
                            #' @description
                            #' Update the model's parameters
                            #' @param coefficients Numeric named vector, model's coefficients. If missing nothing will be updated.
                            update = function(coefficients){
                              if (missing(coefficients)) {
                                return(invisible(NULL))
                              }
                              # Extract old coefficients
                              new_coefs <- self$coefficients
                              # Extract names
                              names_old <- names(new_coefs)
                              names_new <- names(coefficients)
                              # Warning
                              if (length(names_new) != length(names_old)) {
                                cli::cli_alert_warning("In ARMA_model$update(): The lenght of new `coefficients` do not match the length of the old coefficients.")
                              }
                              # Update only if they are present
                              for(i in 1:length(coefficients)){
                                if (names_new[i] %in% names_old) {
                                  new_coefs[names_new[i]] <- coefficients[i]
                                  private$..std.errors[names_new[i]] <- NA_integer_
                                }
                              }
                              # Update intercept
                              if (private$include.intercept){
                                private$..intercept <- new_coefs["intercept"]
                              }
                              # Update AR parameters
                              if (self$order[1] > 0) {
                                private$..phi <- new_coefs[stringr::str_detect(names_old, "phi")]
                              }
                              # Update MA parameters
                              if (self$order[2] > 0) {
                                private$..theta <- new_coefs[stringr::str_detect(names_old, "theta")]
                              }
                              if (!private$include.intercept){
                                new_coefs <- new_coefs[-which(names_old == "intercept")]
                              } else {
                                index <- which(names_old == "intercept")
                                new_coefs <- c(new_coefs[-index], new_coefs[index])
                              }
                              # Update the parameters inside the ARMA model
                              names(new_coefs) <- names(private$..model$coef)
                              private$..model$coef <- new_coefs
                              # Update companion matrix
                              private$..A <- ARMA_companion_matrix(self$phi, self$theta)
                            },
                            #' @description
                            #' Update the standard errors of the parameters.
                            #' @param std.errors Numeric named vector, parameters' standard errors. If missing nothing will be updated.
                            update_std.errors = function(std.errors){
                              if (missing(std.errors) || purrr::is_empty(std.errors)) {
                                return(invisible(NULL))
                              }
                              # Extract old coefficients
                              new_std.errors <- private$..std.errors
                              # Extract names
                              names_old <- names(new_std.errors)
                              names_new <- names(std.errors)
                              # Warning
                              if (length(names_new) != length(names_old)) {
                                cli::cli_alert_warning("In ARMA_model$update_std.errors(): The lenght of new `std.errors` do not match the length of the old std. errors!")
                              }
                              # Update only if they are present
                              for(i in 1:length(std.errors)){
                                if (names_new[i] %in% names_old) {
                                  new_std.errors[names_new[i]] <- std.errors[i]
                                }
                              }
                              # Update private std. errors
                              private$..std.errors <- new_std.errors
                            },
                            #' @description
                            #' Update the variance of the residuals.
                            #' @param sigma2 Numeric scalar, variance of the residuals.
                            update_sigma2 = function(sigma2){
                              if (!missing(sigma2)){
                                private$..sigma2 <- sigma2
                              }
                            },
                            #' @description
                            #' Print method for `AR_modelR6` class.
                            print = function(){
                              green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                              red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                              msg_col <- function(x) ifelse(x, green(x), red(x))
                              # Intercept
                              par.intercept <- format(self$intercept, digits = 4)
                              std.error.intercept <- private$..std.errors[names(private$..std.errors) == names(par.intercept)]
                              format.intercept <- paste0(par.intercept, " (", format(std.error.intercept, digits = 3), ")")
                              # AR parameters
                              par.ar <- format(self$phi, digits = 4)
                              std.error.ar <- private$..std.errors[names(private$..std.errors) %in% names(par.ar)]
                              format.ar <- paste0(par.ar, " (", format(std.error.ar, digits = 3), ")")
                              # MA parameters
                              par.ma <- format(self$theta, digits = 4)
                              std.error.ma <- private$..std.errors[names(private$..std.errors) %in% names(par.ma)]
                              format.ma <- paste0(par.ma, " (", format(std.error.ma, digits = 3), ")")
                              # Model name
                              model_name <- paste0("ARMA", "(", self$order[1], ", ", self$order[2], ")")
                              # ********************************************************
                              msg_0 <- paste0("--------------------- ", model_name, "--------------------- \n")
                              msg_1 <- paste0("Include Intercept: ", msg_col(private$include.intercept), "\n")
                              msg_2 <- paste0("AR: ", msg_col(self$order[1] != 0), "\n")
                              msg_3 <- paste0("MA: ", msg_col(self$order[2] != 0), "\n")
                              msg_4 <- paste0("Version: ", private$version, "\n")
                              msg_5 <- "------------------------------------------------\n"
                              msg_6 <- paste0("Intercept: ", format.intercept, " \n")
                              msg_7 <- paste0("AR parameters: ", format.ar, " \n")
                              msg_8 <- paste0("MA parameters: ", format.ma, " \n")
                              cat(c(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5, msg_6, msg_7, msg_8))
                            }
                          ),
                          private = list(
                            version = "1.0.1",
                            ..model = NA,
                            ..arOrder = 0,
                            ..maOrder = 0,
                            ..intercept = 0,
                            ..phi = c(),
                            ..theta = c(),
                            ..b = c(),
                            ..A = c(),
                            ..sigma2 = 1,
                            ..std.errors = NA,
                            include.intercept = FALSE
                          ),
                          active = list(
                            #' @field model An object with the fitted ARMA model from the function [stats::arima()].
                            model = function(){
                              private$..model
                            },
                            #' @field arOrder Numeric scalar, Autoregressive order.
                            arOrder = function(){
                              private$..arOrder
                            },
                            #' @field maOrder Numeric scalar, Moving Average order.
                            maOrder = function(){
                              private$..maOrder
                            },
                            #' @field order Numeric named vector, orders of the ARMA model. The first element is the AR order, while the second the MA order.
                            order = function(){
                              c(AR = self$arOrder, MA = self$maOrder)
                            },
                            #' @field intercept Numeric named scalar, intercept of the model.
                            intercept = function(){
                              private$..intercept
                            },
                            #' @field phi Numeric named vector, AR parameters.
                            phi = function(){
                              private$..phi
                            },
                            #' @field theta Numeric named vector, MA parameters.
                            theta = function(){
                              private$..theta
                            },
                            #' @field coefficients Numeric named vector, intercept and ARMA parameters.
                            coefficients = function(){
                              c(self$intercept, self$phi, self$theta)
                            },
                            #' @field std.errors Numeric named vector, std.errors of the intercept and ARMA parameters.
                            std.errors = function(){
                              private$..std.errors
                            },
                            #' @field sigma2 Numeric scalar, std.errors of the residuals.
                            sigma2 = function(){
                              private$..sigma2
                            },
                            #' @field A Numeric matrix, companion matrix to govern the transition between two time steps.  See the function [ARMA_companion_matrix()].
                            A = function(){
                              private$..A
                            },
                            #' @field b Numeric vector, unitary vector for the residuals. See the function [ARMA_vector_b()].
                            b = function(){
                              private$..b
                            },
                            #' @field tidy Tibble with estimated parameters and relative std. errors.
                            tidy = function(){
                              dplyr::tibble(
                                term = names(self$coefficients),
                                estimate = self$coefficients,
                                std.error = self$std.errors
                              )
                            }
                          ))

