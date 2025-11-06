#' Kernel regression
#'
#' Fit a kernel regression.
#'
#' @rdname kernelRegression
#' @name kernelRegression
#' @note Version 1.0.0
#' @export
kernelRegression <- R6::R6Class("kernelRegression",
                                public = list(
                                  #' @description
                                  #' Fit a `kernelRegression` class
                                  #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                                  #' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                                  #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                                  #' @param ... other parameters to be passed to the function [np::npreg()].
                                  fit = function(formula, data, ...){
                                    # Model formula
                                    formula <- as.formula(formula)
                                    # Fit a kernel regression
                                    private$..model <- np::npreg(formula, data = data, ...)
                                  },
                                  #' @description
                                  #' Predict method for `kernelRegression` class
                                  #' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
                                  predict = function(newdata){
                                    if (missing(newdata)) {
                                      newdata <- private$..model$eval
                                    }
                                    np:::predict.npregression(private$..model, newdata = newdata)
                                  }
                                ),
                                private = list(
                                  version = "1.0.0",
                                  ..model = NA
                                ),
                                active = list(
                                  #' @field model an object of the class `npreg`.
                                  model = function(){
                                    private$..model
                                  }
                                ))

