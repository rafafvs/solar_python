#' R6 implementation for a clear sky seasonal model
#'
#' @rdname seasonalClearsky
#' @name seasonalClearsky
#' @keywords clearsky
#' @note Version 1.0.1
#' @export
seasonalClearsky <- R6::R6Class("seasonalClearsky",
                                inherit = seasonalModel,
                                public = list(
                                  #' @field lat Numeric, scalar, latitude of the location considered.
                                  lat = NA_integer_,
                                  #' @description
                                  #' Initialize a `seasonalClearsky` object.
                                  #' @param control Named list, control parameters. See the function \code{\link{control_seasonalClearsky}} for more details.
                                  initialize = function(control = control_seasonalClearsky()){
                                    # Store control parameters
                                    private$..control <- control
                                    # Update order
                                    private$..order <- control$order
                                    # Update period
                                    private$..period <- control$period
                                  },
                                  #' @description
                                  #' Fit the seasonal model for clear sky radiation.
                                  #' @param x Numeric vector, time series of solar radiation.
                                  #' @param date Character or Date vector, time series of dates.
                                  #' @param lat Numeric scalar, reference latitude.
                                  #' @param clearsky Numeric vector, time series of target clear sky radiation.
                                  fit = function(x, date, lat, clearsky){
                                    # Self arguments
                                    control = self$control
                                    # Control parameters
                                    include.intercept = control$include.intercept
                                    include.trend = control$include.trend
                                    # Ensure clearsky is specified
                                    if (missing(clearsky)) {
                                      stop('`clearsky` time series must be specified.')
                                    }
                                    # Add the function to compute extraterrestrial radiation
                                    private$..ssf <- seasonalSolarFunctions$new("spencer")
                                    # Store reference latitude
                                    self$lat <- lat[1]
                                    # Initialize the dataset
                                    data <- dplyr::tibble(date = as.Date(date))
                                    data <- dplyr::mutate(data,
                                                          Year = lubridate::year(date),
                                                          Month = lubridate::month(date),
                                                          Day = lubridate::day(date),
                                                          t = Year - max(Year),
                                                          n = number_of_day(date),
                                                          Rt = x,
                                                          H0 = self$ssf$Hon(n, self$lat),
                                                          clearsky = clearsky)
                                    # Method: Estimate with Extraterrestrial and clear sky radiation
                                    # ========================================================================
                                    # 1. Daily maximum clearsky: Ct_max ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
                                    # ========================================================================
                                    # Compute maximum clear sky for each day
                                    base_formula <- "clearsky ~ H0"
                                    if (control$order_H0 > 1){
                                      for(i in 2:control$order_H0){
                                        data[[paste0("H0_", i)]] <- data$H0^i
                                        base_formula <- paste0(base_formula, " + ", paste0("H0_", i))
                                      }
                                    }
                                    base_formula <- ifelse(include.trend, paste0(base_formula, " + t"), base_formula)
                                    base_formula <- ifelse(include.intercept, base_formula, paste0(base_formula, "-1"))
                                    # Fit the coefficients of the clear sky max model
                                    super$fit(formula = base_formula, data = data)
                                    # Initial fit average clear sky
                                    data$Ct_hat <- self$predict(newdata = data)
                                    # Optimize the fit
                                    data <- dplyr::select(data, n, t, H0, Rt, Ct_hat)
                                    # ========================================================================
                                    # 2. Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
                                    # ========================================================================
                                    # Optimize the fit
                                    delta <- clearsky_delta_optimizer(data$Rt, data$Ct_hat*control$delta0, control$lower, control$upper, control$by, control$ntol)
                                    # Standard names for coefficients
                                    coefs_names <- c()
                                    orig_names <- super$.__enclos_env__$private$..model$coefficients_names
                                    if (include.intercept) {
                                      coefs_names[1] <- "delta_0"
                                      orig_names <- orig_names[-c(1)]
                                      for(i in 1:control$order_H0){
                                        coefs_names[i+1] <- paste0("delta_extra", i)
                                        orig_names <- orig_names[-c(1)]
                                      }
                                    } else {
                                      for(i in 1:control$order_H0){
                                        coefs_names[i] <- paste0("delta_extra", i)
                                        orig_names <- orig_names[-c(1)]
                                      }
                                    }
                                    if (include.trend) {
                                      coefs_names <- c(coefs_names, "t")
                                      orig_names <- orig_names[-c(1)]
                                    }
                                    if (self$order > 0) {
                                      coefs_names <- c(coefs_names, paste0("delta_", orig_names))
                                    }

                                    # Store original coefficients
                                    private$coefficients_orig <- self$coefficients
                                    # Store delta parameter
                                    private$delta <- delta * control$delta0
                                    # Update coefficients values and names
                                    super$update(super$.__enclos_env__$private$..model$coefficients * private$delta)
                                    super$.__enclos_env__$private$..model$coefficients_names <- coefs_names
                                    # Update std errors values and names
                                    super$update_std.errors(super$.__enclos_env__$private$..std.errors * private$delta)
                                    names(super$.__enclos_env__$private$..std.errors) <- coefs_names
                                  },
                                  #' @description
                                  #' Predict method for `seasonalClearsky` object.
                                  #' @param n Integer, scalar or vector. number of day of the year.
                                  #' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
                                  predict = function(n, newdata){
                                    if (missing(newdata)) {
                                      if (missing(n)) {
                                        predict.lm(private$..model)
                                      } else {
                                        H0 <- self$ssf$Hon(n, self$lat)
                                        newdata <- data.frame(n = n, H0 = H0)
                                        if (self$control$order_H0 > 1){
                                          for(i in 2:self$control$order_H0){
                                            newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                          }
                                        }
                                        predict.lm(private$..model, newdata = newdata)
                                      }
                                    } else {
                                      newdata$H0 <- self$ssf$Hon(newdata$n, self$lat)
                                      if (self$control$order_H0 > 1){
                                        for(i in 2:self$control$order_H0){
                                          newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                        }
                                      }
                                      predict.lm(private$..model, newdata = newdata)
                                    }
                                  },
                                  #' @description
                                  #' Differential method for `seasonalClearsky` object.
                                  #' @param n Integer, scalar or vector. number of day of the year.
                                  #' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
                                  differential = function(n, newdata){
                                    if (missing(newdata)) {
                                      if (missing(n)) {
                                        predict.lm(private$..dmodel)
                                      } else {
                                        H0 <- self$ssf$Hon(n, self$lat, deriv = TRUE)
                                        newdata <- data.frame(n = n, H0 = H0)
                                        if (self$control$order_H0 > 1){
                                          for(i in 2:self$control$order_H0){
                                            newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                          }
                                        }
                                        predict.lm(private$..dmodel, newdata = newdata)
                                      }
                                    } else {
                                      newdata$H0 <- self$ssf$Hon(newdata$n, self$lat, deriv = TRUE)
                                      if (self$control$order_H0 > 1){
                                        for(i in 2:self$control$order_H0){
                                          newdata[[paste0("H0_", i)]] <- i * newdata$H0^(i-1)
                                        }
                                      }
                                      predict.lm(private$..dmodel, newdata = newdata)
                                    }
                                  },
                                  #' @description
                                  #' Print method for `seasonalClearsky` object.
                                  print = function(){
                                    cat(paste0("----------------------- seasonalClearsky ----------------------- \n"))
                                    msg_1 <- paste0(" - Order: ", self$order, "\n - Period: ", self$period, "\n")
                                    msg_2 <- paste0("- External regressors: 1 (H0) \n")
                                    msg_3 <- paste0("- Version: ", private$version, "\n")
                                    cat(c(msg_1, msg_2, msg_3))
                                    cat(paste0("--------------------------------------------------------------\n"))
                                    print(self$model)
                                  }
                                ),
                                private = list(
                                  version = "1.0.1",
                                  coefficients_orig = NA,
                                  delta = NA,
                                  ..ssf = NA,
                                  ..control = list()
                                ),
                                active = list(
                                  #' @field control Named list, control parameters.
                                  control = function(){
                                    private$..control
                                  },
                                  #' @field ssf Solar Seasonal Functions
                                  ssf = function(){
                                    private$..ssf
                                  }
                                )
                              )
