#' Create a SoRad / SoREd contract specification
#'
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption <- R6::R6Class("solarOption",
                           public = list(
                             #' @field ticker description
                             ticker = "",
                             #' @field strike Strike price for solar radiation.
                             strike = NA,
                             #' @field t_pricing Character, pricing date.
                             t_pricing = NA,
                             #' @field t_now Character, today date.
                             t_now = NA,
                             #' @field t_init Character, inception date.
                             t_init = NA,
                             #' @field t_hor Character, maturity date.
                             t_hor = NA,
                             #' @field tick Numeric, monetary conversion tick.
                             tick = 1,
                             #' @field contract_type Character, maturity date.
                             contract_type = "SoRad",
                             #' @description
                             #' Initialize the contract
                             #' @param contract_type Character, contract type "SoRad" or "SoREd"
                             initialize = function(contract_type = "SoRad"){
                               self$contract_type <- match.arg(contract_type, choices = c("SoRad", "SoREd"))
                             },
                             #' @description
                             #' Initialize the contract
                             #' @param t_pricing Character, pricing date.
                             #' @param t_init Character, inception date.
                             #' @param t_hor Character, maturity date.
                             #' @param strike Numeric, strike price.
                             #' @param tick Numeric monetary tick.
                             set_contract = function(t_pricing, t_init, t_hor, strike, tick = 1){
                               # Conversion in dates
                               self$t_pricing <- self$t_now <- as.Date(t_pricing)
                               self$t_init <- as.Date(t_init)
                               self$t_hor <- as.Date(t_hor)
                               self$strike <- strike
                               self$tick <- tick
                               # Standard ticker name
                               self$ticker <- paste0(self$contract_type, "_", stringr::str_replace_all(self$t_hor, "-", ""))
                             },
                             #' @description
                             #' Store a list of custom control parameters
                             #' @param control List, control parameters.
                             set_control = function(control){
                               private$..control <- control
                             },
                             #' @description
                             #' Print method
                             print = function(){
                               cat("----------- Solar Option  ----------- \n")
                               cat("- Ticker: ", self$ticker, " \n")
                               cat("- Strike: ", self$strike, " \n")
                               cat("- Tick: ", self$tick, " \n")
                               cat("- Today: ", as.character(self$t_now), " \n")
                               cat("- Pricing: ", as.character(self$t_pricing), " \n")
                               cat("- Inception: ", as.character(self$t_init), " \n")
                               cat("- Maturity: ", as.character(self$t_hor), " \n")
                               cat("- Version: ", as.character(private$version), " \n")
                             }
                           ),
                           private = list(
                             version = "1.0.0",
                             ..control = list()
                           ),
                           active = list(
                             #' @field control control parameters
                             control = function(){
                               private$..control
                             },
                             #' @field tau Numeric, scalar. Time from `t_now` till `t_hor` in days.
                             tau = function(){
                               as.numeric(difftime(self$t_hor, self$t_now, units = "days"))
                             },
                             #' @field tau_accrued Numeric, scalar. Time from `t_pricing` till `t_hor` in days.
                             tau_accrued = function(){
                               as.numeric(difftime(self$t_now, self$t_pricing, units = "days"))
                             }
                           ))


#' Create a Solar Option Index portfolio
#'
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
SoRadPorfolio <- function(model, t_now, t_init, t_hor){
  # Pricing date
  t_now <- as.Date(t_now)
  # Inception date
  t_init <- as.Date(t_init)
  # Horizon date
  t_hor <- as.Date(t_hor)
  # Sequence of dates
  t_seq <- seq.Date(t_init, t_hor, 1)
  # Portfolio of contracts
  portfolio <- list()
  for(t in t_seq){
    t_hor <- as.Date(t)
    df_t <- dplyr::filter(model$data, date == t_hor)
    solar_option <- solarOption$new()
    solar_option$set_contract(t_now, t_init, t_hor, df_t$GHI_bar)
    portfolio <- append(portfolio, solar_option)
  }
  names(portfolio) <- purrr::map_chr(portfolio, ~.x$ticker)
  return(portfolio)
}

