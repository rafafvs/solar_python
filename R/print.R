#' Print method for the class `solarOptionPayoff`
#'
#' @param object an object of the class `solarOptionPayoff`.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarOptionPayoff <- function(object){
  msg_1 <- paste0("------------------------ \033[1;35mSolar Option Payoffs\033[0m (", object$payoff_year$side, ") ------------------------ \n")
  msg_2 <- paste0("Yearly payoff: \033[1;31m", format(object$payoff_year$premium, digits = 5), "\033[0m\n")

  # Monthly premiums
  premiums <- object$payoff_month$premium
  # Count the number of integers on the left
  n_integers <- purrr::map_dbl(as.integer(premiums), ~length(stringr::str_split(.x, "")[[1]]))
  n_integers <- ifelse(n_integers == 1, 2, 1)
  # Format the number accordingly
  premiums <- purrr::map2_chr(premiums, n_integers, ~format(.x, digits = 3, nsmall = .y))

  msg_3 <- paste0("Monthly payoffs: \n ",
                  " Jan: \033[1;32m", premiums[1], "\033[0m",
                  "   Feb: \033[1;32m", premiums[2], "\033[0m",
                  "   Mar: \033[1;32m", premiums[3], "\033[0m\n",
                  "  Apr: \033[1;32m", premiums[4], "\033[0m",
                  "   May: \033[1;32m", premiums[5], "\033[0m",
                  "   Jun: \033[1;32m", premiums[6], "\033[0m\n",
                  "  Jul: \033[1;32m", premiums[7], "\033[0m",
                  "   Ago: \033[1;32m", premiums[8], "\033[0m",
                  "   Sep: \033[1;32m", premiums[9], "\033[0m\n",
                  "  Oct: \033[1;32m", premiums[10], "\033[0m",
                  "   Nov: \033[1;32m", premiums[11], "\033[0m",
                  "   Dec: \033[1;32m", premiums[12], "\033[0m\n")
  cat(paste0(msg_1, msg_2, msg_3))
}

#' Print method for the class `solarOptionChoquet`
#'
#' @param object an object of the class  \code{\link{solarOption_choquet}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarOptionChoquet <- function(object){
  msg_1 <- paste0("------------------------ \033[1;35mSolar Option Bootstrap\033[0m (", object$payoff_year$side, ") ------------------------ \n")
  msg_2 <- paste0("Yearly payoff: \033[1;31m", format(object$payoff_year$premium, digits = 5), "\033[0m\n")
  msg_3 <- paste0("Quantile ", "(", names(object$payoff_year$premium_dw), ")",
                  "\033[1;31m ", format(object$payoff_year$premium_dw, digits = 5), "\033[0m", "\n",
                  "Quantile ", "(", names(object$payoff_year$premium_up), ")",
                  "\033[1;31m ", format(object$payoff_year$premium_up, digits = 5), "\033[0m \n")

  format_premiums <- function(premiums){
    # Count the number of integers on the left
    n_integers <- purrr::map_dbl(as.integer(premiums), ~length(stringr::str_split(.x, "")[[1]]))
    n_integers <- ifelse(n_integers == 1, 2, 1)
    # Format the number accordingly
    premiums <- purrr::map2_chr(premiums, n_integers, ~format(.x, digits = 3, nsmall = .y))
    premiums
  }

  # Monthly premiums
  premiums <- format_premiums(object$payoff_month$premium)
  # Monthly premiums (up)
  premiums_up <- format_premiums(object$payoff_month$premium_up)
  # Monthly premiums (dw)
  premiums_dw <- format_premiums(object$payoff_month$premium_dw)
  msg_line <- paste0(paste0(rep("-", 78), collapse = ""), "\n")
  msg_4 <- paste0("Monthly payoffs: \n ",
                  " Jan: \033[1;32m", premiums[1], "\033[0m", paste0(" (\033[1;31m", premiums_dw[1], "\033[0m", " - ", "\033[1;31m", premiums_up[1], "\033[0m)"),
                  "   Feb: \033[1;32m", premiums[2], "\033[0m", paste0(" (\033[1;31m", premiums_dw[2], "\033[0m", " - ", "\033[1;31m", premiums_up[2], "\033[0m)"),
                  "   Mar: \033[1;32m", premiums[3], "\033[0m", paste0(" (\033[1;31m", premiums_dw[3], "\033[0m", " - ", "\033[1;31m", premiums_up[3], "\033[0m) \n"),
                  "  Apr: \033[1;32m", premiums[4], "\033[0m", paste0(" (\033[1;31m", premiums_dw[4], "\033[0m", " - ", "\033[1;31m", premiums_up[4], "\033[0m)"),
                  "   May: \033[1;32m", premiums[5], "\033[0m", paste0(" (\033[1;31m", premiums_dw[5], "\033[0m", " - ", "\033[1;31m", premiums_up[5], "\033[0m)"),
                  "   Jun: \033[1;32m", premiums[6], "\033[0m", paste0(" (\033[1;31m", premiums_dw[6], "\033[0m", " - ", "\033[1;31m", premiums_up[6], "\033[0m) \n"),
                  "  Jul: \033[1;32m", premiums[7], "\033[0m", paste0(" (\033[1;31m", premiums_dw[7], "\033[0m", " - ", "\033[1;31m", premiums_up[7], "\033[0m)"),
                  "   Ago: \033[1;32m", premiums[8], "\033[0m", paste0(" (\033[1;31m", premiums_dw[8], "\033[0m", " - ", "\033[1;31m", premiums_up[8], "\033[0m)"),
                  "   Sep: \033[1;32m", premiums[9], "\033[0m", paste0(" (\033[1;31m", premiums_dw[9], "\033[0m", " - ", "\033[1;31m", premiums_up[9], "\033[0m) \n"),
                  "  Oct: \033[1;32m", premiums[10], "\033[0m", paste0(" (\033[1;31m", premiums_dw[10], "\033[0m", " - ", "\033[1;31m", premiums_up[10], "\033[0m)"),
                  "   Nov: \033[1;32m", premiums[11], "\033[0m", paste0(" (\033[1;31m", premiums_dw[11], "\033[0m", " - ", "\033[1;31m", premiums_up[11], "\033[0m)"),
                  "   Dec: \033[1;32m", premiums[12], "\033[0m", paste0(" (\033[1;31m", premiums_dw[12], "\033[0m", " - ", "\033[1;31m", premiums_up[12], "\033[0m)"))
  cat(paste0(msg_1, msg_2, msg_3, msg_line, msg_4))
}


#' Print method for the class `spatialModel`
#'
#' @param object an object of the class \code{\link{spatialModel}}.
#'
#' @keywords internal
#' @noRd
print_spatialModel <- function(object){
  range_lat <- range(object$locations$lat)
  range_lon <- range(object$locations$lon)
  n_points <- nrow(object$locations)
  n_models <- length(object$models)
  n_params <- length(object$params_models$models)

  msg_1 <- paste0("Spatial model: ", "Points", " (", n_points, ")",
                  " - Models ", "(", n_models, ")",
                  " - Parameters", " (", n_params, ") \n")
  msg_2 <- paste0("    Latitude: ", range_lat[1], " - ", range_lat[2], " - ",
                  "Longitude: ", range_lon[1], " - ", range_lon[2], " \n")
  cat(paste0(msg_1, msg_2))
}

#' Print method for the class `spatialScenarioSpec`
#'
#' @param object an object of the class \code{\link{spatialScenario_spec}}.
#'
#' @keywords internal
#' @noRd
#' @export
print_spatialScenarioSpec <- function(object){
  cat("------------------------- spatialScenarioSpec -------------------------\n")
  msg_1 <- paste0(" Number of models: ", length(object$spec), "\n")
  msg_2 <- paste0("Dates: ", object$from, " - ", object$to, "\n")
  msg_3 <- paste0("Number of simulations: ", object$nsim, "\n",
                  " - Residuals: (", object$residuals, ")", "\n",
                  " - Filter: (", object$filter, ")")
  cat(msg_1, msg_2, msg_3)
}


