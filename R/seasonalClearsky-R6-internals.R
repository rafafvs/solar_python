#' Control parameters for a `seasonalClearsky` object
#'
#' @param include.intercept Logical, when `TRUE`, the default, the intercept will be included in the clear sky model.
#' @param order Integer scalar, number of combinations of sines and cosines.
#' @param period Integer scalar, seasonality period. The default is 365.
#' @param delta0 Numeric scalar, initial value for the optimization. The estimated clear sky is increased by delta0.
#' @param quiet Logical, when `FALSE`, the default, the functions displays warning or messages.
#' @inheritParams clearsky_optimizer
#' @details The parameters `ntol`, `lower`, `upper` and `by` are used exclusively in \code{\link{clearsky_optimizer}}.
#' @examples
#' control = control_seasonalClearsky()
#' @return Named list of control parameters.
#' @keywords control
#' @note Version 1.0.0
#' @rdname control_seasonalClearsky
#' @export
control_seasonalClearsky <- function(order = 1, order_H0 = 1, period = 365, include.intercept = TRUE, include.trend = FALSE,
                                     delta0 = 1.4, lower = 0, upper = 3, by = 0.001, ntol = 0, quiet = FALSE){
  structure(
    list(
      order = order,
      order_H0 = order_H0,
      period = period,
      include.intercept = include.intercept,
      include.trend = include.trend,
      delta0 = delta0,
      lower = lower,
      upper = upper,
      by = by,
      ntol = ntol,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param x Numeric vector, time series of solar radiation.
#' @param Ct Numeric vector, time series of  clear sky radiation.
#' @param lower Numeric scalar, lower bound for grid of delta parameters used for optimization. Default is `0`.
#' @param upper Numeric scalar, upper bound for grid of delta parameters used for optimization. Default is `3`.
#' @param by Numeric scalar, step for the grid. Default is `0.01`.
#' @param ntol Integer scalar, Tolerance for the maximum number of violations admitted of the condition `Ct > x`. The default is `0`.
#' @details
#' Detect the best parameter `delta` such that the condition `x > delta * Ct` is satisfied for all the observations except for `ntol` observations.
#'
#' @name clearsky_delta_optimizer
#' @rdname clearsky_delta_optimizer
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
#' @noRd
clearsky_delta_optimizer <- function(x, Ct, lower = 0, upper = 3, by = 0.01, ntol = 0){
  # Grid of points
  grid <- seq(lower, upper, by = by)
  # Loss
  opt <- data.frame(delta = grid, loss = purrr::map_dbl(grid, ~sum(.x*Ct - x < 0)))
  # Return minimum delta satisfying the constraint
  delta <- opt[which(opt$loss <= ntol)[1],]$delta
  return(delta)
}

#' Optimizer for clear sky model with restricted least squares (RLS).
#'
#' @param seasonal_model_Ct An object of the class `seasonalClearsky`. See the function \code{\link{seasonalClearsky}} for more details.
#' @param newdata A data frame to input to the method `seasonal_model_Ct$predict()`.
#' @param ntol Integer scalar. Tolerance for the maximum number of violations admitted of the condition `clearsky > GHI`. Default is `0`.
#' @name clearsky_optimizer
#' @rdname clearsky_optimizer
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_optimizer <- function(seasonal_model_Ct, newdata, ntol = 0){
  # Clone the model
  sm <- seasonal_model_Ct$clone(TRUE)
  # Loss function
  loss_function <- function(params, sm, data){
    sm$update(params)
    # Prediction
    pred <- sm$predict(newdata = data)
    # Violations
    violations <- ifelse(data$GHI - pred > 0, 1, 0)
    # Check number of violations lower than ntol
    mse <- sum((data$GHI - pred)^2) + 1000000*(sum(violations) - ntol)
    return(mse)
  }
  # Optimal parameters
  opt <- optim(sm$coefficients, loss_function, sm = sm, data = newdata)
  # Update the parameters
  sm$update(opt$par)
  return(sm)
}

#' Impute clear sky outliers
#'
#' Detect and impute outliers with respect to a maximum level of radiation (Ct)
#'
#' @param x Numeric vector, time series of solar radiation.
#' @param Ct Numeric vector, time series of fitted clear sky radiation.
#' @param date Character or Date vector, time series of dates used to precisely impute solar radiation according to the realized values in the same day of the year.
#' @param threshold Numeric, scalar, threshold value used for imputation. Default is `0.0001`.
#' @param quiet Logical.
#' @examples
#' clearsky_outliers(c(1,2,3), 2)
#' @details
#' The function will detect the observations for which `x > Ct`, `x < 0` or `is.na(x)`. Then, if
#' \describe{
#' \item{`x < 0`}{If a value is below 0 for a day it will be imputed to be equal to `min(x)` for that day}.
#' \item{`x > Ct`}{If a value is above the maximum clear sky Ct it will be imputed to be `Ct*(1-threshold)`}.
#' \item{`is.na(x)`}{If a value is NA it will be imputed to be the average `mean(x)` for that day.}.
#' }
#'
#' @rdname clearsky_outliers
#' @name clearsky_outliers
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_outliers <- function(x, Ct, date, threshold = 0.0001, quiet = FALSE){
  # Initialize a dataset
  data <- dplyr::tibble(Ct = Ct, x = x)
  # Eventually add a date for non-stationary data
  if (!missing(date)){
    data <- dplyr::mutate(data,
                          date = as.Date(date),
                          Month = lubridate::month(date),
                          Day = lubridate::day(date))
  }
  # Number of observations
  nobs <- nrow(data)
  # Detect problems and violations
  outliers_na <- which(is.na(data$x))
  outliers_lo <- which(data$x <= 0 & !is.na(data$x))
  outliers_hi <- which(data$x >= data$Ct & !is.na(data$x))
  # Complete outliers index
  idx_outliers <- unique(c(outliers_na, outliers_lo, outliers_hi))
  # Check presence of outliers and impute them
  if (purrr::is_empty(idx_outliers)) {
    if (!quiet) message("No outliers!")
    data_clean <- data
  } else {
    # Verbose message
    if (!quiet) message("Outliers: ", length(idx_outliers), " (", format(length(idx_outliers)/nobs*100, digits = 2), " %)")
    # Dataset without outliers
    data_no_outliers <- data[-c(idx_outliers),]
    # Imputed dataset
    data_clean <- data
    # Impute outliers
    for (i in idx_outliers) {
      df_n <- data[i,]
      if (!missing(date)) {
        # Data for the same day and month (without outliers)
        df_day <- dplyr::filter(data_no_outliers, Month == df_n$Month & Day == df_n$Day & date != df_n$date)
      } else {
        df_day <- data_no_outliers
      }
      # Imputed data depending on outliers type
      if (i %in% outliers_na) {
        # Outlier is an NA
        data_clean[i,]$x <- mean(df_day$x, na.rm = TRUE)
      } else if (i %in% outliers_lo) {
        # Outlier is under minimum value for the day
        data_clean[i,]$x <- min(df_day$x, na.rm = TRUE)
      } else if (i %in% outliers_hi) {
        # Outlier is above maximum value for the day
        data_clean[i,]$x <- data_clean$Ct[i] * (1-threshold)
      }
    }
  }

  if (missing(date)) {
    date <- NA
  } else {
    date <- data$date[idx_outliers]
  }

  # Structure output data
  structure(
    list(
      # Complete imputed series
      x = data_clean$x,
      # Values before imputation
      original = data$x[idx_outliers],
      # Values after imputation
      imputed = data_clean$x[idx_outliers],
      # Index of imputed observations
      index = idx_outliers,
      # Index of type of outlier
      index_type = list(na = outliers_na, lo = outliers_lo, hi = outliers_hi),
      # Dates of imputed observations
      date = date,
      # Number of outliers
      n = length(idx_outliers),
      # Error of imputed series with respect to original data
      MAPE = mean(abs((data$x - data_clean$x)/data$x))*100,
      MSE = sd(data$x - data_clean$x),
      # Threshold for imputation
      threshold = threshold
    )
  )
}

