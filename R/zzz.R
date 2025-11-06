#' Detect the season
#'
#' Detect the season from a vector of dates
#'
#' @param x vector of dates in the format `YYYY-MM-DD`.
#' @param invert logica, when `TRUE` the seasons will be inverted.
#'
#' @examples
#' detect_season("2040-01-31")
#' detect_season(c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
#' @name detect_season
#' @rdname detect_season
#'
#' @return a character vector containing the correspondent season. Can be `spring`, `summer`, `autumn`, `winter`.
#' @export
detect_season <- function(x, invert = FALSE){
  # first day of the year
  date_start <- as.Date("2000-01-01")
  # start dates for the seasons
  date_spring <- as.Date("2000-03-23")
  date_summer <- as.Date("2000-06-23")
  date_autumn <- as.Date("2000-09-23")
  date_winter <- as.Date("2000-12-23")
  # last day of the year
  date_end <- as.Date("2000-12-31")

  seasons <- dplyr::tibble(
    date = as.Date(x),
    Month = lubridate::month(date),
    Day = lubridate::day(date)
  )
  # Standardized year
  seasons <- dplyr::mutate(seasons, date = as.Date(paste0("2000-", Month, "-", Day)))
  seasons <- dplyr::mutate(seasons,
                           season = dplyr::case_when(
                             date < date_spring  & date >= date_start ~ "Winter",
                             date >= date_spring & date < date_summer ~ "Spring",
                             date >= date_summer & date < date_autumn ~ "Summer",
                             date >= date_autumn & date < date_winter ~ "Autumn",
                             date >= date_winter & date <= date_end ~ "Winter"))
  if (invert) {
    seasons <- dplyr::mutate(seasons,
                             season = dplyr::case_when(
                               season == "Winter" & invert ~ "Summer",
                               season == "Summer" & invert ~ "Winter",
                               season == "Spring" & invert ~ "Autumn",
                               season == "Autumn" & invert ~ "Spring",
                               TRUE ~ season))
  }
  return(seasons$season)
}

#' Is leap year?
#'
#' Check if a given year is leap (366 days) or not (365 days).
#'
#' @param x numeric value or dates vector in the format `YYYY-MM-DD`.
#'
#' @examples
#' is_leap_year("2024-02-01")
#' is_leap_year(c(2023:2030))
#' is_leap_year(c("2024-10-01", "2025-10-01"))
#' is_leap_year("2029-02-01")
#' @name is_leap_year
#' @rdname is_leap_year
#' @return Boolean. `TRUE` if it is a leap year, `FALSE` otherwise.
#' @export
is_leap_year <- function(x){
  if (is.numeric(x)) {
    date_year <- x
  } else {
    date_year <- lubridate::year(as.Date(x))
  }
  return(date_year %% 4 == 0)
}

#' Number of day
#'
#' Compute the number of day of the year given a vector of dates.
#'
#' @param x dates vector in the format `YYYY-MM-DD`.
#'
#' @examples
#' number_of_day("2040-01-31")
#' number_of_day(c("2015-03-31", "2016-03-31", "2017-03-31"))
#' number_of_day(c("2015-02-28", "2016-02-28", "2017-02-28"))
#' number_of_day(c("2015-03-01", "2016-03-01", "2017-03-01"))
#' @name number_of_day
#' @rdname number_of_day
#' @return Numeric vector with the number of the day during the year.
#' @export
number_of_day <- function(x){
  if (any(is.numeric(x))){
    return(x)
  }
  n_of_day <- lubridate::yday(x)
  #n_of_day[is_leap_year(x)] <- n_of_day[is_leap_year(x)] - 1
  #n_of_day[lubridate::day(x) == 29 & lubridate::month(x) == 2] <- 59.5
  return(n_of_day)
}

#' Power of a matrix
#'
#' Compute the power of a matrix
#' @param x Matrix
#' @param n power, if zero will return the identity matrix.
#' @param eigen_dec Logical, when `TRUE`, the default, the power of the matrix is
#' computed with an eigenvalue decomposition.
#' @rdname pow_matrix
#' @name pow_matrix
#' @useDynLib solarr, .registration = TRUE
#' @export
#' @noRd
pow_matrix <- function(x, n = 0, eigen_dec = FALSE){
  .Call("pow_matrix_c", as.matrix(x), as.integer(n))
}

#' Format the pvalue
#'
#' @param x pvalue
#' @rdname format_pval
#' @name format_pval
#' @export
#' @noRd
format_pval <- function(x){
  ifelse(x > 0.1, "", ifelse(x > 0.05, "*", ifelse(x > 0.01, "**", "***")))
}

