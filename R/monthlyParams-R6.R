#' Create a function of time for monthly parameters
#'
#' @param params Vector of length 12 with the monthly parameters.
#'
#' @examples
#' set.seed(1)
#' params <- runif(12)
#' mp <- monthlyParams$new(params)
#' t_now <- as.Date("2022-01-01")
#' t_hor <- as.Date("2024-12-31")
#' dates <- seq.Date(t_now, t_hor, by = "1 day")
#' plot(mp$predict(dates), type = "l")
#' @note Version 1.0.0
#' @export
monthlyParams <- R6::R6Class("monthlyParams",
                             public = list(
                               #' @description
                               #' Initialize a `monthlyParams` object
                               #' @param params numeric vector of parameters with length 12.
                               initialize = function(params){
                                 if (length(params) != 12) {
                                   stop("The length of the vector of parameters must be 12!")
                                 }
                                 private$..parameters <- params
                               },
                               #' @description
                               #' Predict the monthly paramete
                               #' @param x date as character or month as numeric.
                               predict = function(x){
                                 nmonths <- lubridate::month(x)
                                 par <- c()
                                 for(i in 1:length(nmonths)){
                                   par[i] <- self$parameters[nmonths[i]]
                                 }
                                 return(par)
                               },
                               #' @description
                               #' Update the monthly parameters
                               #' @param params numeric vector of parameters with length 12.
                               update = function(params){
                                 if (length(params) != 12) {
                                   stop("The length of the vector of parameters must be 12!")
                                 }
                                 private$..parameters <- params
                               }
                             ),
                             private = list(
                               version = "1.0.0",
                               ..parameters = NA
                             ),
                             active = list(
                               #' @field parameters vector of parameters with length 12.
                               parameters = function(){
                                 private$..parameters
                               }
                             ))
