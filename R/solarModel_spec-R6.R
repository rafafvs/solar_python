#' Control function for a `solarModel` object
#'
#' Control function for a `solarModel` object that contains all the setups used for the estimation.
#'
#' @examples
#' control <- solarModel_spec$new()
#' @rdname solarModel_spec
#' @name solarModel_spec
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_spec <- R6::R6Class("solarModel_spec",
                               public = list(
                                 #' @description
                                 #' Initialize a `solarModel_spec` object.
                                 initialize = function(){
                                   self$set_params()
                                   self$set_transform()
                                   self$set_clearsky()
                                   self$set_seasonal.mean()
                                   self$set_mean.model()
                                   self$set_seasonal.variance()
                                   self$set_variance.model()
                                   self$set_mixture.model()
                                 },
                                 #' @description
                                 #' Specification function for a `solarModel`
                                 #' @param place Character, name of an element in the `CAMS_data` list.
                                 #' @param target Character, target variable to model. Can be `GHI` or `clearsky`.
                                 #' @param min_date Character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
                                 #' @param max_date Character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
                                 #' @param from Character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
                                 #' If `missing` will be used the minimum data available after filtering for `min_date`.
                                 #' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
                                 #' If `missing` will be used the maximum data available after filtering for `max_date`.
                                 #' @param data data for the selected location.
                                 specification = function(place, target = "GHI", min_date, max_date, from, to, data){
                                   # Match the target variable to model
                                   target <- match.arg(target, choices = c("GHI", "clearsky"))
                                   # Extract CAMS data for the selected location
                                   if (missing(data)){
                                     # Match a location in the dataset
                                     place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
                                     data <- solarr::CAMS_data[[place]]
                                   }
                                   # Minimum date for the complete data
                                   if (missing(min_date) || is.null(min_date) || is.na(min_date)) {
                                     min_date <- min(data$date, na.rm = TRUE)
                                   } else {
                                     min_date <- as.Date(min_date)
                                   }
                                   # Maximum date for the complete data
                                   if (missing(max_date) || is.null(max_date) || is.na(max_date)) {
                                     max_date <- max(data$date, na.rm = TRUE)
                                   } else {
                                     max_date <- as.Date(max_date)
                                   }
                                   # Minimum date for train data
                                   if (missing(from) || is.null(from) || is.na(from)) {
                                     from <- min(data$date, na.rm = TRUE)
                                   } else {
                                     from <- as.Date(from)
                                   }
                                   # Maximum date for train data
                                   if (missing(to) || is.null(to) || is.na(to)) {
                                     to <- max(data$date, na.rm = TRUE)
                                   } else {
                                     to <- as.Date(to)
                                   }
                                   # Filter for min and max dates the complete dataset
                                   data <- dplyr::filter(data, date >= min_date & date <= max_date)
                                   # Increase clearsky to avoid NaN
                                   data$clearsky <- data$clearsky * self$clearsky_threshold
                                   # It may happen in CAMS data that clear sky value is just a little bit greater than GHI (~10-3)
                                   # Therefore, before using such time series it is convenient to impute GHI value such that they
                                   # became equal to the given CAMS clear sky
                                   # Detect and impute outliers
                                   outliers <- clearsky_outliers(data$GHI, data$clearsky, date = data$date, threshold = 0, quiet = self$quiet)
                                   # Update the time series of GHI with adjusted values
                                   data$GHI <- outliers$x
                                   # Label for data used for estimation
                                   data <- dplyr::mutate(data,
                                                         isTrain = ifelse(date >= from & date <= to, TRUE, FALSE),
                                                         weights = ifelse(isTrain, 1, 0))
                                   # Add the normalized weights
                                   data$weights <-  data$weights / sum(data$weights)
                                   # Train observations and percentage
                                   nobs_train <- length(data$isTrain[data$isTrain])
                                   # Compute percentage of train obs. on total obs.
                                   perc_train <- nobs_train / nrow(data)
                                   # Train observations and percentage
                                   nobs_test <- length(data$isTrain[!data$isTrain])
                                   # Compute percentage of test obs. on total obs.
                                   perc_test <- nobs_test / nrow(data)
                                   # Model dates
                                   model_dates = list(data = list(from = min_date, to = max_date, nobs = nrow(data), perc = 1),
                                                      train = list(from = from, to = to, nobs = nobs_train, perc = perc_train),
                                                      test = list(from = to, to = max_date, nobs = nobs_test, perc = perc_test))
                                   # Store the data
                                   private$..place = attr(data, "place")
                                   private$..coords = attr(data, "coords")
                                   private$..data = data
                                   private$..dates = model_dates
                                   private$..target = target
                                 },
                                 #' @description
                                 #' Generic controls
                                 #' @param stochastic_clearsky Logical, when `TRUE` the clear sky will be considered stochastic.
                                 #' @param clearsky_threshold Numeric, parameter > 1, used to scale up CAMS clearsky to avoid that clear sky radiaion and global horizontal radiation are equal.
                                 #' @param quiet Logical. When `TRUE` the function will not display any message. The dafault if `TRUE`.
                                 set_params = function(stochastic_clearsky = FALSE, clearsky_threshold = 1.01, quiet = FALSE){
                                   private$..stochastic_clearsky = stochastic_clearsky
                                   private$..clearsky_threshold = clearsky_threshold
                                   private$..quiet = quiet
                                 },
                                 #' @description
                                 #' Control parameters for the `solarTransform`. See \code{\link{solarTransform}} for more details.
                                 #' @param min_pos Integer, position of the minimum. For example when `2` the minimum is the second lowest value.
                                 #' @param max_pos Integer, position of the maximum. For example when `3` the maximum is the third greatest value.
                                 #' @param link Character, link function.
                                 #' @param delta transform params
                                 #' @param threshold Numeric. Threshold used to estimate the transformation parameters \deqn{\alpha} and \deqn{\beta}.
                                 #' The default is `0.01`. See \code{\link{solarTransform}} for more details.
                                 set_transform = function(min_pos = 1, max_pos = 1, link = "invgumbel", delta = 0.05, threshold = 0.01){
                                   private$..transform$min_pos = min_pos
                                   private$..transform$max_pos = max_pos
                                   private$..transform$delta = delta
                                   private$..transform$threshold = threshold
                                   private$..transform$link = match.arg(link, choices = c("invgumbel", "gumbel", "logis", "norm"))
                                 },
                                 #' @description
                                 #' List with specification's parameters of the clear sky model.
                                 #' @param control Named list with control parameters. See \code{\link{control_seasonalClearsky}} for more details.
                                 set_clearsky = function(control = control_seasonalClearsky()){
                                   private$..clearsky <- control
                                 },
                                 #' @description
                                 #' List with specification's parameters of the seasonal mean \eqn{\bar{Y}_t} for \eqn{Y_t}.
                                 #' @param order Integer. Specify the order of the seasonal mean \deqn{\bar{Y}_t}. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{a_0} will be included in the seasonal model, otherwise will be excluded. The default is `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly means will be computed on the deseasonalized series \deqn{\tilde{Y}_t = Y_t - \bar{Y}_t}
                                 #'  and it is subtracted to ensure that the time series is centered around zero for all the months. The dafault if `TRUE`.
                                 set_seasonal.mean = function(order = 1, period = 365, include.trend = FALSE, include.intercept = TRUE, monthly.mean = TRUE){
                                   private$..seasonal.mean = list(order = order, period = period,
                                                                  include.trend = include.trend, include.intercept = include.intercept,
                                                                  monthly.mean = monthly.mean)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the ARMA model for deseasonalized series \eqn{\tilde{Y}_t = Y_t - \bar{Y}_t}.
                                 #' @param arOrder Integer. An integer specifying the order of the AR component. The default is `1`.
                                 #' @param maOrder Integer. An integer specifying the order of the MA component. The default is `0`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{\phi_0} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 set_mean.model = function(arOrder = 1, maOrder = 0, include.intercept = FALSE){
                                   private$..mean.model = list(arOrder = arOrder, maOrder = maOrder, include.intercept = include.intercept)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the seasonal variance \eqn{\bar{\sigma}_t} for ARMA's residuals \eqn{e_t}
                                 #' @param order Integer. Specify the order of the seasonality of the seasonal variance. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param correction Logical. When `TRUE` the parameters of seasonal variance are corrected to ensure
                                 #'  that the standardize the residuals have exactly a unitary variance. The dafault if `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly std. deviations will be computed
                                 #'  on the standardized residuals  \deqn{\tilde{\varepsilon}_t} and used to standardize the time series
                                 #'  such that it has unitary variance for all the months. The default if `TRUE`.
                                 set_seasonal.variance = function(order = 1, period = 365, include.trend = FALSE, correction = TRUE, monthly.mean = TRUE){
                                   private$..seasonal.variance = list(order = order, period = period, correction = correction,
                                                                      include.trend = include.trend, monthly.mean = monthly.mean)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the GARCH variance \eqn{\sigma_t} for deseasonalized residuals \eqn{\tilde{e}_t = e_t/\bar{\sigma}_t}.
                                 #' @param archOrder Integer. An integer specifying the order of the ARCH component. The default is `1`.
                                 #' @param garchOrder Integer. An integer specifying the order of the GARCH component. The default is `1`.
                                 #' @param garch_variance Logical. When `TRUE` the GARCH model will be used to standardize the residuals otherwise will be excluded. The dafault if `TRUE`.
                                 set_variance.model = function(archOrder = 1, garchOrder = 1, garch_variance = TRUE){
                                   private$..garch_variance = garch_variance
                                   # Modify only if variance is not false
                                   if (archOrder == 0 & garchOrder == 0){
                                     private$..garch_variance <- FALSE
                                   } else {
                                     private$..garch_variance <- TRUE
                                   }
                                   private$..variance.model = list(archOrder = archOrder, garchOrder = garchOrder)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the Gaussian mixture model for GARCH residuals \eqn{u_t = \tilde{e}_t/\sigma_t}.
                                 #' @param abstol Numeric. Absolute level for convergence of the EM-algorithm. The default is `1e-20`.
                                 #' @param match.expectation Logical, when `TRUE` the mixture parameters ensures that the expected value is matched.
                                 #' @param match.variance Logical, when `TRUE` the mixture parameters ensures that the variance is matched.
                                 #' @param match.empiric Logical, when `TRUE` and `match.expectation = TRUE` or  `match.variance = TRUE` the mixture parameters
                                 #' will be estimated ensuring that mean and variance matches the empirical parameters. Otherwise if `FALSE` and
                                 #'  `match.expectation = TRUE` or `match.variance = TRUE` the target expectation will be zero and the target variance 1.
                                 #' @param method Character, package used to fit the parameters. Can be `mclust` or `mixtools`.
                                 #' @param maxit Integer. Maximum number of iterations for EM-algorithm. The default is `5000`.
                                 #' @param maxrestarts Integer. Maximum number of restarts when EM-algorithm does not converge. The default is `500`.
                                 set_mixture.model = function(abstol = 1e-20, match.expectation = TRUE, match.variance = FALSE,
                                                              match.empiric = FALSE, method = "mclust", maxit = 5000, maxrestarts = 500){
                                   method <- match.arg(method, choices = c("mclust", "mixtools"))
                                   private$..mixture.model = list(abstol = abstol,
                                                                  method = method,
                                                                  match.expectation = match.expectation,
                                                                  match.variance = match.variance,
                                                                  match.empiric = match.empiric,
                                                                  maxit = maxit,
                                                                  maxrestarts = maxrestarts)
                                 },
                                 #' @description
                                 #' Print method for `solarModel_spec` class.
                                 print = function(){
                                   # Seasonal mean order
                                   sm_order <- self$seasonal.mean$order
                                   # ARMA order
                                   arOrder <- self$mean.model$arOrder
                                   maOrder <- self$mean.model$maOrder
                                   # Seasonal variance order
                                   sv_order <- self$seasonal.variance$order
                                   # GARCH order
                                   archOrder <- self$variance.model$archOrder
                                   garchOrder <- self$variance.model$garchOrder
                                   green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                                   red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                                   msg_col <- function(x) ifelse(x, green(x), red(x))

                                   if (!is.na(self$place)) {
                                     # Complete data specifications
                                     data <- self$dates$data
                                     # Train data specifications
                                     train <- self$dates$train
                                     train$perc <- format(train$perc*100, digits = 4)
                                     # Test data specifications
                                     test <- self$dates$test
                                     test$perc <- format(test$perc*100, digits = 4)
                                     msg_0 <- paste0("--------------------- ", "solarModel", " (\033[4;35m", self$place, "\033[0m) ", "--------------------- \n")
                                     msg_1 <- paste0("\033[1;35m Target \033[0m: ", self$target, " \n\033[1;35m Coordinates\033[0m: ",
                                                     "(\033[1;35mLat\033[0m: ", self$coords$lat,
                                                     ", \033[1;35mLon\033[0m: ", self$coords$lon,
                                                     ", \033[1;35mAlt\033[0m: ", self$coords$alt, ") \n")
                                     msg_2 <- paste0("\033[1;35m Observations\033[0m: ", data$nobs, "\n")
                                     msg_3 <- paste0("---------------------------------------------------------------\n")
                                     msg_4 <- paste0("\033[1;34m Dates\033[0m: ", data$from, " - ", data$to, "\n")
                                     msg_5 <- paste0("  - \033[1;34mTrain\033[0m: ", train$from, " - ", train$to, " (", train$nobs, " points ~ ", train$perc, "%)", "\n")
                                     msg_6 <- paste0("  - \033[1;34mTest\033[0m: ", test$from, " - ", test$to, " (", test$nobs, " points ~ ", test$perc, "%)", "\n")
                                     cat(paste0(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5, msg_6))
                                   }
                                   msg_0 <- "------------------------\033[4;35m Specification \033[0m------------------------ \n"
                                   model_name <- paste0("S(", sm_order, ", ", sv_order, ")-ARMA(", arOrder, ", ", maOrder, ")")
                                   msg_1 <- paste0(" - Stochastic clearsky: ", msg_col(self$stochastic_clearsky), "\n")
                                   if (self$garch_variance) {
                                     model_name <- paste0(model_name,"-GARCH(", archOrder, ", ", archOrder, ")")
                                   }
                                   cat(c(msg_0, paste0(model_name, "\n"), msg_1))
                                   msg_0 <- "-------------------------\033[4;35m Mean Models \033[0m------------------------- \n"
                                   msg_1 <- paste0(" - Trend: ", msg_col(self$seasonal.mean$include.trend), "\n")
                                   msg_2 <- paste0(" - Intercept: ", msg_col(self$seasonal.mean$include.intercept), "\n")
                                   msg_3 <- paste0(" - Monthly correction: ", msg_col(self$seasonal.mean$monthly.mean), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3))
                                   msg_0 <- "-----------------------\033[4;35m Variance Models \033[0m----------------------- \n"
                                   msg_1 <- paste0(" - Trend: ", msg_col(self$seasonal.variance$include.trend), "\n")
                                   msg_2 <- paste0(" - Correction: ", msg_col(self$seasonal.variance$correction), "\n")
                                   msg_3 <- paste0(" - Monthly correction: ", msg_col(self$seasonal.variance$monthly.mean), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3))
                                   msg_0 <- "------------------------\033[4;35m Mixture Model \033[0m------------------------ \n"
                                   msg_1 <- paste0(" - Method: ", self$mixture.model$method, "\n")
                                   msg_2 <- paste0(" - Match Expectation: ", msg_col(self$mixture.model$match.expectation), "\n")
                                   msg_3 <- paste0(" - Match Variance: ", msg_col(self$mixture.model$match.variance), "\n")
                                   msg_4 <- paste0(" - Match Empiric: ", msg_col(self$mixture.model$match.empiric), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3, msg_4))
                                 }
                               ),
                               private = list(
                                 ..place = NA,
                                 ..coords = NA,
                                 ..data = NA,
                                 ..dates = NA,
                                 ..target = NA,
                                 ..transform = list(),
                                 ..clearsky = list(),
                                 ..seasonal.mean = list(),
                                 ..mean.model = list(),
                                 ..seasonal.variance = list(),
                                 ..variance.model = list(),
                                 ..mixture.model = list(),
                                 ..garch_variance = FALSE,
                                 ..stochastic_clearsky = FALSE,
                                 ..clearsky_threshold = 1.01,
                                 ..quiet = FALSE
                               ),
                               active = list(
                                 #' @field place Character, optional name of the location considered.
                                 place = function(){
                                   private$..place
                                 },
                                 #' @field target Character, name of the target variable to model. Can be `"GHI"` or `"clearsky"`.
                                 target = function(){
                                   private$..target
                                 },
                                 #' @field coords A named list with the coordinates of the location considered. Contains:
                                 #' \describe{
                                 #'  \item{lat}{Numeric, reference latitude in degrees.}
                                 #'  \item{lon}{Numeric, reference longitude in degrees.}
                                 #'  \item{alt}{Numeric, reference altitude in metres.}
                                 #'}
                                 coords = function(){
                                   private$..coords
                                 },
                                 #' @field dates A named list, with three sub-lists: `data` containing the information on the complete dataset,
                                 #' `train` containing the information on the train dataset and `test` containing the information on the test dataset.
                                 #' Each sub-list is structured as follows:
                                 #' \describe{
                                 #'  \item{from}{Character date, minmum date in the dataset.}
                                 #'  \item{to}{Character date, maximum date in the dataset.}
                                 #'  \item{nobs}{Integer scalar, number of observations contained in the dataset between `from` and `to`.}
                                 #'  \item{perc}{Numeric scalar, percentage of data in the dataset with respect to the complete data.}
                                 #'}
                                 dates = function(){
                                   private$..dates
                                 },
                                 #' @field data Tibble, dataset with CAMS solar radiation data.
                                 data = function(){
                                   private$..data
                                 },
                                 #' @field transform Named list, specification of the solar transform.
                                 transform = function(){
                                   private$..transform
                                 },
                                 #' @field clearsky Named list, specification of the clear sky model.
                                 clearsky = function(){
                                   private$..clearsky
                                 },
                                 #' @field seasonal.mean Named list, specification of the seasonal model.
                                 seasonal.mean = function(){
                                   private$..seasonal.mean
                                 },
                                 #' @field mean.model Named list, specification of the ARMA model.
                                 mean.model = function(){
                                   private$..mean.model
                                 },
                                 #' @field seasonal.variance Named list, specification of the seasonal variance model.
                                 seasonal.variance = function(){
                                   private$..seasonal.variance
                                 },
                                 #' @field variance.model Named list, specification of the GARCH model for deseasonalized residuals \eqn{\tilde{e}_t}.
                                 variance.model = function(){
                                   private$..variance.model
                                 },
                                 #' @field mixture.model Named list, specification of the Mixture model for GARCH residuals \eqn{u_t}.
                                 mixture.model = function(){
                                   private$..mixture.model
                                 },
                                 #' @field garch_variance Logical, when `TRUE` the GARCH model will be used otherwise no.
                                 garch_variance = function(){
                                   private$..garch_variance
                                 },
                                 #' @field clearsky_threshold Numeric, parameter > 1, used to scale up CAMS clearsky.
                                 clearsky_threshold = function(){
                                   private$..clearsky_threshold
                                 },
                                 #' @field stochastic_clearsky Logical, when `TRUE` the clear sky is considered stochastic.
                                 stochastic_clearsky = function(){
                                   private$..stochastic_clearsky
                                 },
                                 #' @field quiet Logical. When `TRUE` the function will not display any message. The dafault if `TRUE`.
                                 quiet = function(){
                                   private$..quiet
                                 }
                               ))

