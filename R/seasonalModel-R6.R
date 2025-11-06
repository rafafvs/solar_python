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
#' @note Version 1.0.1
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
                                 base_formula <- base_formula_dt <- formula
                                 if (length(self$order) == 1) {
                                   for(period in self$period) {
                                     for (order in self$order){
                                       for(i in 1:order){
                                         base_formula <- seasonalModel_formula(base_formula, order = i, period = period)
                                         base_formula_dt <- seasonalModel_formula_dt(base_formula_dt, order = i, period = period)
                                       }
                                     }
                                   }
                                 } else if (length(self$order) == length(self$period)) {
                                   for(j in 1:length(self$period)) {
                                     period <- self$period[j]
                                     for (i in 1:self$order[j]){
                                       base_formula <- seasonalModel_formula(base_formula, order = i, period = period)
                                       base_formula_dt <- seasonalModel_formula_dt(base_formula_dt, order = i, period = period)
                                     }
                                   }
                                 }
                                 # Store the main formula
                                 private$mformula <- as.formula(base_formula)
                                 attr(private$mformula, "coef_names") <- attr(base_formula, "coef_names")
                                 # Fit seasonal model
                                 private$..model <- lm(private$mformula, data = data)

                                 # Fit the differential w.r.t. t
                                 private$dformula <- as.formula(base_formula_dt)
                                 attr(private$dformula, "coef_names") <- attr(base_formula_dt, "coef_names")
                                 private$..dmodel <- lm(private$dformula, data = data)
                                 # Update coefficients
                                 dcoefs <- private$..model$coefficients
                                 names(dcoefs) <- names(private$..dmodel$coefficients)
                                 dcoefs[stringr::str_detect(names(dcoefs), "sin")] <- -dcoefs[stringr::str_detect(names(dcoefs), "sin")]
                                 dcoefs[stringr::str_detect(names(dcoefs), "(Intercept)")] <- 0
                                 private$..dmodel$coefficients <- dcoefs

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
                               #' Predict method for the class `seasonalModel`.
                               #' @param n Integer vector, numbers of day of the year.
                               #' @param newdata an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               predict = function(n, newdata){
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
                               differential = function(n, newdata){
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
                                 # Extract old coefficients
                                 new_coefs <- self$coefficients
                                 # Extract names
                                 names_old <- names(new_coefs)
                                 names_new <- names(coefficients)
                                 # Warning
                                 if (length(names_new) != length(names_old)) {
                                   cli::cli_alert_warning("In seasonalModel$update(): The lenght of new `coefficients` do not match the length of the old coefficients.")
                                 }
                                 # Update only if they are present
                                 for(i in 1:length(coefficients)){
                                   if (names_new[i] %in% names_old) {
                                     new_coefs[names_new[i]] <- coefficients[i]
                                     private$..std.errors[names_new[i]] <- NA_integer_
                                   }
                                 }
                                 # Set the names equal to the original one
                                 names(new_coefs) <- names(private$..model$coefficients)
                                 # Update the parameters inside the lm object
                                 private$..model$coefficients <- new_coefs

                                 # Update coefficients of the differential
                                 dcoefs <- new_coefs
                                 names(dcoefs) <- names(private$..dmodel$coefficients)
                                 dcoefs[stringr::str_detect(names(dcoefs), "sin")] <- -dcoefs[stringr::str_detect(names(dcoefs), "sin")]
                                 dcoefs[stringr::str_detect(names(dcoefs), "(Intercept)")] <- 0
                                 private$..dmodel$coefficients <- dcoefs
                               },
                               #' @description
                               #' Update the parameter's std. errors.
                               #' @param std.errors Named vector, new standard errors of the parameters.
                               update_std.errors = function(std.errors){
                                 # Extract old coefficients
                                 new_std.errors <- self$std.errors
                                 # Extract names
                                 names_old <- names(new_std.errors)
                                 names_new <- names(std.errors)
                                 # Warning
                                 if (length(names_new) != length(names_old)) {
                                   cli::cli_alert_warning("In seasonalModel$update_std.errors(): The lenght of new `std.errors` do not match the length of the old std. errors!")
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
                               version = "1.0.1",
                               ..model = NA,
                               ..dmodel = NA,
                               mformula = NA,
                               dformula = NA,
                               ..period = 1,
                               ..order = 365,
                               ..std.errors = c(),
                               external_regressors = NA
                             ),
                             active = list(
                               #' @field coefficients Named vector, fitted coefficients.
                               coefficients = function(){
                                 coefs <- private$..model$coefficients
                                 names(coefs) <- private$..model$coefficients_names
                                 return(coefs)
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
                               #' @field omega Integer, periodicity.
                               omega = function(){
                                 2 * base::pi / private$..period
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
