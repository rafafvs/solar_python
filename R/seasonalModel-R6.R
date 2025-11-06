#' Seasonal Model
#'
#' @description
#' The `seasonalModel` class implements a seasonal regression model as a linear combination of sine and cosine functions.
#' This model is designed to capture periodic effects in time series data, particularly for applications involving seasonal trends.
#'
#' @details
#' The seasonal model is fitted using a specified formula, which allows for the inclusion of external regressors along with sine and cosine terms
#' to model seasonal variations. The periodicity can be customized, and the model can be updated with new coefficients after the initial fit.
#'
#' @examples
#' sm <- seasonalModel$new(1, 365)
#' formula <- "Yt ~ 1"
#' data = data.frame(Yt = rnorm(1000), n = 1:1000)
#' sm$fit(formula, data = data)
#' sm
#' sm$coefficients
#' sm$update(sm$coefficients*3)
#' sm$predict(20)
#'
#' @rdname seasonalModel
#' @name seasonalModel
#' @keywords seasonalModel
#' @note Version 1.0.0
#' @export
seasonalModel <- R6::R6Class("seasonalModel",
                             public = list(
                               #' @field extra_params List to contain custom extra parameters.
                               extra_params = list(),
                               #' @description
                               #' Initialize a `seasonalModel` object.
                               #' @param order Integer, number of combinations of sines and cosines.
                               #' @param period Integer, seasonality period. The default is 365.
                               initialize = function(order = 1, period = 365){
                                 # Store period and order
                                 private$..period = period
                                 private$..order = order
                               },
                               #' @description
                               #' Fit a seasonal model as a linear combination of sine and cosine functions and
                               #' eventual external regressors specified in the formula. The external regressors used should
                               #' have the same periodicity, i.e. not stochastic regressors are allowed.
                               #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                               #' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
                               #' @param data  an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param ... other parameters to be passed to the function `lm`.
                               fit = function(formula, data, ...){
                                 # Formula with standard names
                                 base_formula <- formula
                                 if (length(self$order) == 1) {
                                   for(period in self$period) {
                                     for (order in self$order){
                                       for(i in 1:order){
                                         base_formula <- formula_fourier(base_formula, order = i, period = period)
                                       }
                                     }
                                   }
                                 } else if (length(self$order) == length(self$period)) {
                                   for(j in 1:length(self$period)) {
                                     period <- self$period[j]
                                     for (i in 1:self$order[j]){
                                       base_formula <- formula_fourier(base_formula, order = i, period = period)
                                     }
                                   }
                                 }
                                 # Store the formula
                                 private$mformula <- as.formula(base_formula)
                                 attr(private$mformula, "coef_names") <- attr(base_formula, "coef_names")
                                 # Fit seasonal model
                                 private$..model <- lm(private$mformula, data = data)#, ...)
                                 # Detect number of regressors
                                 n_regressors <- length(private$..model$coefficients)
                                 # Extract regressors from the formula excluding target variable
                                 regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 # Seasonal regressors
                                 idx_seasonal_regressors <- which(stringr::str_detect(regressors, "sin|cos"))
                                 n_seasonal_reg <- length(idx_seasonal_regressors)
                                 # External regressors
                                 idx_external_regressors <- which(!stringr::str_detect(regressors, "sin|cos"))
                                 private$external_regressors <- regressors[idx_external_regressors]
                                 n_external_regressors <- length(idx_external_regressors)
                                 # Standard names
                                 coefs_names <- regressors
                                 coefs_names[idx_seasonal_regressors] <- attr(base_formula, "coef_names")
                                 coefs_names[idx_external_regressors] <- regressors[idx_external_regressors]
                                 # Intercept
                                 if (n_regressors - n_seasonal_reg - n_external_regressors == 1) {
                                   coefs_names <- c("Intercept", coefs_names)
                                 }
                                 names(coefs_names) <- NULL
                                 # Store coefficients names
                                 private$..model$coefficients_names <- coefs_names
                                 # Update std.errors
                                 private$..std.errors <- broom::tidy(private$..model)$std.error
                                 names(private$..std.errors) <- names(self$coefficients)
                               },
                               #' @description
                               #' Fit the differential of the sinusoidal function.
                               #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                               #' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
                               #' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param ... other parameters to be passed to the function `lm`.
                               fit_differential = function(formula, data, ...){
                                 # Extract seasonal regressors
                                 regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 seasonal_regressors <- which(stringr::str_detect(regressors, "sin|cos")) + 1
                                 external_regressors <- which(!stringr::str_detect(regressors, "sin|cos")) + 1
                                 # Reparametrize the coefficients
                                 coefs <- private$..model$coefficients[c(1, seasonal_regressors)]
                                 index <- seq(1, length(coefs)-1, by = 2)
                                 param <- c(private$..model$coefficients[c(1, external_regressors)])
                                 names(param) <- c("A", names(private$..model$coefficients[external_regressors]))
                                 k <- 1
                                 i <- 1
                                 for(i in index){
                                   coefs_names <- names(param)
                                   param <- c(param, sqrt(coefs[i + 1]^2 + coefs[i + 2]^2))
                                   param <- c(param, atan(coefs[i + 2]/coefs[i + 1]) + ifelse(coefs[i + 1] < 0, -base::pi, 0))
                                   names(param) <- c(coefs_names, paste0(c("B", "phi"), k))
                                   k <- k + 1
                                 }
                                 # Reparametrize the coefficients
                                 coefs <- param
                                 coef_phi <- coefs[stringr::str_detect(names(coefs), "phi")]
                                 coef_B <- coefs[stringr::str_detect(names(coefs), "B")]
                                 # Formula with standard names
                                 base_formula <- paste0(formula.tools::lhs(private$mformula), " ~ - 1")
                                 k <- 1
                                 if (length(self$order) == 1) {
                                   for(period in self$period) {
                                     for (order in self$order){
                                       for(i in 1:order){
                                         base_formula <- formula_fourier(base_formula, order = i, reparam = TRUE, period = period, phi = coef_phi[k])
                                         k <- k + 1
                                       }
                                     }
                                   }
                                 } else if (length(self$order) == length(self$period)) {
                                   for(j in 1:length(self$period)) {
                                     period <- self$period[j]
                                     for (i in 1:self$order[j]){
                                       base_formula <- formula_fourier(base_formula, order = i, reparam = TRUE, period = period, phi = coef_phi[k])
                                       k <- k + 1
                                     }
                                   }
                                 }
                                 # Store the formula
                                 private$dformula <- as.formula(paste0(base_formula, " - 1"))
                                 attr(private$mformula, "coef_names") <- paste0("B_", attr(base_formula, "coef_names"))
                                 # Fit seasonal model
                                 private$..dmodel <- lm(private$dformula, data = data)
                                 names(coef_B) <- names(private$..dmodel$coefficients)
                                 private$..dmodel$coefficients <- coef_B
                               },
                               #' @description
                               #' Predict method for the class `seasonalModel`.
                               #' @param n Integer vector, numbers of day of the year.
                               #' @param newdata an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param dt Numeric, time step.
                               predict = function(n, newdata, dt = 1){
                                 if (missing(newdata)) {
                                   if (missing(n)) {
                                     predict.lm(private$..model)
                                   } else {
                                     predict.lm(private$..model, newdata = data.frame(n = n))
                                   }
                                 } else {
                                   predict.lm(private$..model, newdata = newdata)
                                 }
                               },
                               #' @description
                               #' Compute the differential of the sinusoidal function.
                               #' @param n Integer, number of day of the year.
                               #' @param newdata an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param dt Numeric, time step.
                               differential = function(n, newdata, dt = 1){
                                 if (missing(newdata)) {
                                   if (missing(n)) {
                                     predict.lm(private$..dmodel)
                                   } else {
                                     predict.lm(private$..dmodel, newdata = data.frame(n = n))
                                   }
                                 } else {
                                   predict.lm(private$..dmodel, newdata = newdata)
                                 }
                               },
                               #' @description
                               #' Update the model's parameters.
                               #' @param coefficients Named vector, new parameters.
                               update = function(coefficients){
                                 old_coef <- self$coefficients
                                 # Check length
                                 condition <- length(old_coef) == length(coefficients)
                                 if (!condition) {
                                   cli::cli_alert_danger("In seasonalModel$update(): The lenght of new `coefficients` do not match the length of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Check names
                                 condition <- names(old_coef) %in% names(coefficients)
                                 condition <- sum(condition) == length(coefficients)
                                 if (!condition) {
                                   cli::cli_alert_warning("In seasonalModel$update(): The names of new `coefficients` do not match the names of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Ensure that the parameters are correctly ordered
                                 new_coef <- coefficients[names(old_coef)]
                                 # Set std.errors of the parameters that are changed equal to NA
                                 idx_new_coefs <- which(new_coef != old_coef)
                                 if (!purrr::is_empty(idx_new_coefs)) {
                                   private$..std.errors[idx_new_coefs] <- NA
                                 }
                                 # Set the names equal to the original one
                                 names(new_coef) <- names(private$..model$coefficients)
                                 # Update the parameters inside the lm object
                                 private$..model$coefficients <- new_coef
                                 # Set std.errors of the parameters that are changed equal to NA
                                 idx_new_coefs <- which(new_coef != old_coef)
                                 if (!purrr::is_empty(idx_new_coefs)) {
                                   private$..std.errors[idx_new_coefs] <- NA
                                 }
                               },
                               #' @description
                               #' Update the parameter's std. errors.
                               #' @param std.errors Named vector, new standard errors of the parameters.
                               update_std.errors = function(std.errors){
                                 if (!missing(std.errors)) {
                                   # Update the vector of std. errors
                                   std.errors_updated <- self$coefficients
                                   std.errors_updated[!(names(std.errors_updated) %in% names(std.errors))] <- NA_integer_
                                   std.errors_updated[names(std.errors_updated) %in% names(std.errors)] <- std.errors
                                   # Update private std. errors
                                   private$..std.errors <- std.errors_updated
                                 }
                               },
                               #' @description
                               #' Print method for the class `seasonalModel`.
                               print = function(){
                                 cat(paste0("----------------------- seasonalModel ----------------------- \n"))
                                 msg_1 <- paste0(" - Order: ", self$order, "\n - Period: ", self$period, "\n")
                                 if (any(is.na(private$external_regressors))) {
                                   msg_2 <- paste0("- External regressors: 0 \n")
                                 } else {
                                   n_external_regressors <- length(private$external_regressors)
                                   msg_2 <- paste0("- External regressors: ", n_external_regressors, " (", private$external_regressors, ")\n")
                                 }
                                 msg_3 <- paste0("- Version: ", private$version, "\n")
                                 cat(c(msg_1, msg_2, msg_3))
                                 cat("--------------------------------------------------------------\n")
                                 print(self$model)
                               }
                             ),
                             private = list(
                               version = "1.0.0",
                               ..model = NA,
                               ..dmodel = NA,
                               mformula = NA,
                               dformula = NA,
                               ..period = 1,
                               ..order = 365,
                               ..std.errors = c(),
                               ..coefficients2 = NA,
                               external_regressors = NA
                             ),
                             active = list(
                               #' @field coefficients Named vector, fitted coefficients.
                               coefficients = function(){
                                 coefs <- private$..model$coefficients
                                 names(coefs) <- private$..model$coefficients_names
                                 return(coefs)
                               },
                               #' @field coefficients2 Named vector, reparametrized coefficients into a linear combination of shifted sine functions.
                               coefficients2 = function(){
                                 private$..coefficients2
                               },
                               #' @field model A slot with the fitted `lm` object.
                               model = function(){
                                 private$..model
                               },
                               #' @field period Integer scalar, periodicity of the seasonality.
                               period = function(){
                                 private$..period
                               },
                               #' @field order Integer scalar, number of combinations of sines and cosines.
                               order = function(){
                                 private$..order
                               },
                               #' @field std.errors Named vector, with the parameters' std. errors.
                               std.errors = function(){
                                 private$..std.errors
                               },
                               #' @field tidy A tibble with estimated parameters and std. errors.
                               tidy = function(){
                                 dplyr::tibble(
                                   term = names(self$coefficients),
                                   estimate = self$coefficients,
                                   std.error = private$..std.errors
                                 )
                               }
                             )
)

#' Create a custom fourier formula
#'
#' @keywords seasonalModel
#' @noRd
formula_fourier <- function(formula, reparam = FALSE, order = 1, period = 365, phi = 0, sin = TRUE, cos = TRUE, t_idx = "n"){
  if (order == 0) {
    return(formula)
  }
  coefs_names <- c()
  base_formula <- formula
  if (!reparam) {
    if (sin) {
      # Standard name
      new_coef_name <- paste0("sin_", order, "_", period)
      # Check if it was not included
      if (!new_coef_name %in% attr(base_formula, "coef_names")){
        base_formula <- paste0(base_formula, " + ", "I(sin(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
        coefs_names <- c(coefs_names, new_coef_name)
      }
    }
    if (cos) {
      # Standard name
      new_coef_name <- paste0("cos_", order, "_", period)
      # Check if it was not included
      if (!new_coef_name %in% attr(base_formula, "coef_names")){
        base_formula <- paste0(base_formula, " + ", "I(cos(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
        coefs_names <- c(coefs_names, new_coef_name)
      }
    }
  } else {
    # Standard name
    new_coef_name <- paste0("cos_", order, "_", period)
    # Check if it was not included
    if (!new_coef_name %in% attr(base_formula, "coef_names")){
      base_formula <- paste0(base_formula, " + ", "I(2 * base::pi / ", eval(period),
                             "*cos(2 * base::pi / ", eval(period), " * ", t_idx, " * ", eval(order), " + ", eval(phi),"))")
      coefs_names <- c(coefs_names, new_coef_name)
    }
  }
  attr(base_formula, "coef_names") <- c(attr(formula, "coef_names"), coefs_names)
  return(base_formula)
}
