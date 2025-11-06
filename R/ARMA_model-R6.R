#' R6 class for ARMA(p, q) model
#'
#' @rdname ARMA_modelR6
#' @name ARMA_modelR6
#' @keywords ARMA
#' @note Version 1.0.0
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
                              private$include.intercept <- include.intercept
                              private$..arOrder <- arOrder
                              private$..maOrder <- maOrder
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
                              # Extract intercept
                              intercept <- c(intercept = 0)
                              if (private$include.intercept) {
                                index <- which(names(params) == "intercept")
                                intercept <- c(intercept = params[[index]])
                                params <- params[-c(index)]
                              }
                              # AR coefficients without intercept
                              phi <- c()
                              if (self$order[1] > 0){
                                phi <- params[stringr::str_detect(names(params), "ar")]
                                names(phi) <- paste0("phi_", 1:self$order[1])
                              }
                              # MA coefficients without intercept
                              theta <- c()
                              if (self$order[2] > 0){
                                theta <- params[stringr::str_detect(names(params), "ma")]
                                names(theta) <- paste0("theta_", 1:self$order[2])
                              }
                              # Extract the std.errors
                              std.errors <- broom::tidy(ARMA_model)$std.error
                              names(std.errors) <- names(params)
                              # *********************************************
                              # Store coefficients
                              private$..intercept <- intercept
                              private$..phi <- phi
                              private$..theta <- theta
                              # Store the fitted model
                              private$..model <- ARMA_model
                              # Store the std. errors
                              private$..std.errors <- std.errors
                            },
                            #' @description
                            #' Filter the time-series and compute fitted values and residuals.
                            #' @param x Numeric vector, time series to filter.
                            #' @param eps0 Numeric vector, initial residuals of the same length of the MA order.
                            filter = function(x, eps0){
                              # Maximum order
                              k <- max(self$order)
                              # Length of the time series
                              n <- length(x)
                              # Vector to store the fitted residuals
                              e_hat <- eps0
                              # Vector to store the fitted time series
                              x_hat <- x
                              # Initialize the state vector
                              x_t <- c(x[k:(k-self$order[1]+1)], eps0[self$order[2]:1])
                              # i <- k + 1
                              b <- self$b
                              A <- self$Phi
                              for(i in (k + 1):n){
                                x_t <- ARMA_next_step(x_t, A, b, 1, self$intercept)
                                # Fitted series
                                x_hat[i] <- x_t[1,1]
                                # Fitted residuals
                                e_hat[i] <- x[i] - x_hat[i]
                                # Update state vector
                                x_t <- x_t + b * e_hat[i]
                              }
                              dplyr::tibble(fitted = x_hat, residuals = e_hat)
                            },
                            #' @description
                            #' Next step function
                            #' @param x Numeric vector, state vector with past observations and residuals.
                            #' @param n.ahead Numeric scalar, forecasted steps ahead.
                            #' @param eps Numeric vector, optional realized residuals.
                            next_step = function(x, n.ahead = 1, eps = 0){
                              ARMA_next_step(x, self$Phi, self$b, n.ahead, self$intercept, eps)
                            },
                            #' @description
                            #' Update the model's parameters
                            #' @param coefficients Numeric named vector, model's coefficients. If missing nothing will be updated.
                            update = function(coefficients){
                              # 1) Update the parameters
                              if (!missing(coefficients)) {
                                names(coefficients) <- names(self$coefficients)
                                # Intercept
                                if(private$include.intercept){
                                  # Update intercept
                                  private$..intercept <- c(intercept = coefficients[[1]])
                                } else {
                                  # Set intercept equal to zero
                                  private$..intercept <- c(intercept = 0)
                                  # Remove intercept from coefficients
                                  coefficients <- coefficients[-c(1)]
                                }
                                # Update the parameters inside the ARMA model
                                private$..model$coef <- coefficients
                                # Update AR parameters
                                if (self$order[1] > 0) {
                                  private$..phi <- coefficients[stringr::str_detect(names(coefficients), "phi")]
                                }
                                # Update MA parameters
                                if (self$order[2] > 0) {
                                  private$..theta <- coefficients[stringr::str_detect(names(coefficients), "theta")]
                                }
                                # Set std. errors equal to NA
                                private$..std.errors <- rep(NA, length(coefficients))
                              }
                            },
                            #' @description
                            #' Update the standard errors of the parameters.
                            #' @param std.errors Numeric named vector, parameters' standard errors. If missing nothing will be updated.
                            update_std.errors = function(std.errors){
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
                              msg_7 <- paste0("Phi: ", format.ar, " \n")
                              msg_8 <- paste0("Theta: ", format.ma, " \n")
                              cat(c(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5, msg_6, msg_7, msg_8))
                            }
                          ),
                          private = list(
                            version = "1.0.0",
                            ..model = NA,
                            ..arOrder = 0,
                            ..maOrder = 0,
                            ..intercept = 0,
                            ..phi = c(),
                            ..theta = c(),
                            ..std.errors = NA,
                            include.intercept = FALSE
                          ),
                          active = list(
                            #' @field model An object with the fitted ARMA model from the function [stats::arima()].
                            model = function(){
                              private$..model
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
                            #' @field order Numeric named vector, ARMA order. The first element is the AR order, while the second the MA order.
                            order = function(){
                              c(AR = private$..arOrder, MA = private$..maOrder)
                            },
                            #' @field coefficients Numeric named vector, intercept and ARMA parameters.
                            coefficients = function(){
                              c(self$intercept, self$phi, self$theta)
                            },
                            #' @field mean Numeric scalar, long term expectation.
                            mean = function(){
                              self$intercept / (1 - sum(self$phi))
                            },
                            #' @field variance Numeric scalar, long term variance.
                            variance = function(){
                              ARMA_variance(self$Phi, self$b, 1, 1000)
                            },
                            #' @field Phi Numeric matrix, companion matrix to govern the transition between two time steps.
                            Phi = function(){
                              ARMA_companion_matrix(self$phi, self$theta)
                            },
                            #' @field b Numeric vector, unitary vector for the residuals.
                            b = function(){
                              ARMA_vector_b(self$order[1], self$order[2])
                            },
                            #' @field tidy Tibble with estimated parameters and relative std. errors.
                            tidy = function(){
                              dplyr::tibble(
                                term = names(self$coefficients),
                                estimate = self$coefficients,
                                std.error = private$..std.errors
                              )
                            }
                          ))

