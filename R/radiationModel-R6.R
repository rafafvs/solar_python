#' Radiation model
#'
#' @rdname radiationModel
#' @name radiationModel
#' @note Version 1.0.0
#' @export
radiationModel <- R6::R6Class("radiationModel",
                              # ====================================================================================================== #
                              #                                             Public slots
                              # ====================================================================================================== #
                              public = list(
                                #' @field theta Numeric, mean reversion parameter.
                                theta = NA,
                                #' @field lambda_S Numeric, market risk premium Q-measure.
                                lambda_S = 0,
                                #' @description
                                #' Initialize a `radiationModel` object
                                #' @param model A `solarModel` object. See \code{\link{solarModel}}.
                                #' @param correction Logical. When `TRUE` the mixture means will be corrected for the
                                #' discrepancy between the integral seasonal std. deviation and variance.
                                initialize = function(model, correction = FALSE){
                                  # Store the discrete model
                                  private$..model <- model$clone(TRUE)
                                  # 1) Estimate mean reversion parameter
                                  df_fit <- filter(private$..model$data, isTrain & weights != 0)
                                  self$theta <- martingale_method_seasonal(df_fit$Yt, df_fit$Yt_bar)
                                  # Convert AR parameter into mean reversion parameter (alternative)
                                  # self$theta <- -log(self$model$ARMA$coefficients[2])
                                  private$..model$update(c(phi_1 = exp(-self$theta)))
                                  private$..model$filter()
                                  private$..model$fit_seasonal_variance()
                                  private$..model$filter()
                                  private$..model$fit_NM_model()
                                  private$..model$update_NM_classification()
                                  private$..model$update_moments()
                                  private$..model$update_logLik()
                                  # 2) Reparametrize seasonal function with continuous time parameters
                                  self$parametrize_seasonal_variance()
                                  # 3) Update mixture parameters
                                  if (correction){
                                    self$correct_NM_coefficients()
                                  }
                                },
                                #' @description
                                #' Compute the parameters of the seasonal variance given OLS estimates.
                                parametrize_seasonal_variance = function(){
                                  # Reparametrize seasonal function with continuous time parameters
                                  reparam <- reparam_seasonal_function(self$model$seasonal_variance$coefficients, self$theta, omega = 2*base::pi/365)
                                  # Clone seasonal variance to store c parameters
                                  private$seasonal_variance <- self$model$seasonal_variance$clone(TRUE)
                                  private$seasonal_variance$extra_params$reparam <- reparam
                                  names(reparam$c_) <- names(self$model$seasonal_variance$coefficients)
                                  private$seasonal_variance$update(reparam$c_)
                                  private$..integral_variance <- integral_sigma2_formula(self$theta, reparam$gamma, omega = 2*base::pi/365)
                                  private$..integral_expectation <- integral_sigma_numeric(self$theta, reparam$c_, omega = 2*base::pi/365)
                                },
                                #' @description
                                #' Correct for the discrepancy between the integral seasonal std. deviation and variance.
                                correct_NM_coefficients = function(){
                                  # Adjustment for the mean that multiplies I
                                  t <- 1:365
                                  sigma_bar_J <- sqrt(private$..integral_variance(t-1, t, t))
                                  sigma_bar_I <- private$..integral_expectation(t-1, t, t)
                                  private$..k1 <- mean(sigma_bar_I/sigma_bar_J)
                                  # Original variance
                                  v <- self$model$NM_model$moments$variance
                                  # Original parameters
                                  sigma <- self$model$NM_model$sd
                                  probs <- self$model$NM_model$p
                                  # Adjusted means
                                  means <- self$model$NM_model$means * private$..k1
                                  # New expectation
                                  m1 <- means[,1] * probs[,1] + means[,2] * probs[,2]
                                  k2 <- (v + m1 - (means[,1] - means[,2])^2 * probs[,1] * probs[,2]) / (sigma[,1]^2 * probs[,1] + sigma[,2]^2 * probs[,2])
                                  # Adjusted std. dev
                                  private$..k2 <- sqrt(k2)
                                  sigma <- sigma * private$..k2
                                  # Update mixture parameters
                                  private$..model$NM_model$update(means = means)#, sd = sigma)
                                },
                                #' @description
                                #' Change the reference probability measure
                                #' @param measure Character, probability measure. Can be `P` or `Q`.
                                change_measure = function(measure){
                                  measure <- match.arg(measure, choices = c("P", "Q"))
                                  private$..measure <- measure
                                  if (measure == "Q"){
                                    private$..lambda <- self$lambda_S
                                  } else {
                                    private$..lambda <- 0
                                  }
                                },
                                #' @description
                                #' Clear sky radiation for a day of the year.
                                #' @param t_now Character, today date.
                                #' @return Clear sky radiation at time t_now.
                                Ct = function(t_now){
                                  self$model$seasonal_model_Ct$predict(number_of_day(t_now))
                                },
                                #' @description
                                #' Seasonal mean for the transformed variable \eqn{Y_t} for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for \eqn{Y_t} at time t_now.
                                Yt_bar = function(t_now){
                                  self$model$seasonal_model_Yt$predict(number_of_day(t_now))
                                },
                                #' @description
                                #' Seasonal mean for the solar radiation for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for Rt.
                                Rt_bar = function(t_now){
                                  self$model$Y_to_R(self$Yt_bar(t_now), t_now)
                                },
                                #' @description
                                #' Transformed variable instantaneous seasonal std. deviation \eqn{\bar{\sigma_{t}}}.
                                #' @param t_now Character, today date.
                                #' @return Seasonal std. deviation for Yt on date t_now.
                                sigma_bar = function(t_now){
                                  sqrt(private$seasonal_variance$predict(number_of_day(t_now)))
                                },
                                #' @description
                                #' Transformed variable mixture mean drift \eqn{\mu_(B)}.
                                #' @param t_now Character, today date.
                                #' @param B Integer. If `B = 1` it is returned the mean of the first component,
                                #' otherwise if `B = 1` the second.
                                #' @return Mixture seasonal drift for \eqn{Y_t} at time t_now.
                                mu_B = function(t_now, B = 1){
                                  result <- ifelse(B == 1, self$model$NM_model$mu1$predict(t_now), self$model$NM_model$mu2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Transformed variable mixture diffusion drift \eqn{\sigma_(B)}
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Mixture seasonal diffusion for \eqn{Y_t}.
                                sigma_B = function(t_now, B){
                                  result <- ifelse(B == 1, self$model$NM_model$sd1$predict(t_now), self$model$NM_model$sd2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Transformed variable drift \eqn{\mu_(Y)}.
                                #' @param Yt Numeric, transformed solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @param dt Numeric, time step.
                                #' @return Mixture drift for \eqn{Y_t}.
                                mu_Y = function(Yt, t_now, B = 1, dt = 1){
                                  # Number of the day
                                  n <- number_of_day(t_now)
                                  # Differential for seasonal function
                                  dYt_bar <- self$model$seasonal_model_Yt$predict(n+dt) - self$model$seasonal_model_Yt$predict(n)
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(n)
                                  # Drift for Yt
                                  drift_Y <- dYt_bar - self$theta * (Yt - self$Yt_bar(n)) + sigma_bar * self$mu_B(n, B)
                                  # !!! Change not already implemented for QT measure
                                  drift_Y + sigma_bar * self$sigma_B(n, B) * self$lambda * ifelse(self$measure == "Q", 1, 0)
                                },
                                #' @description
                                #' Transformed variable diffusion \eqn{\sigma_(Y)}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for \eqn{Y_t}.
                                sigma_Y = function(t_now, B = 1){
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(t_now)
                                  # Diffusion for Rt process
                                  sigma_bar * self$sigma_B(t_now, B)
                                },
                                #' @description
                                #' Solar radiation drift \eqn{\mu_(R)}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @param dt Numeric, time step.
                                #' @return Drift for \eqn{R_t}.
                                mu_R = function(Rt, t_now, B = 1, dt = 1){
                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Approximate the drift for Ct
                                  n <- number_of_day(t_now)
                                  Ct <- self$Ct(n)
                                  # Differential
                                  dCt <- self$Ct(n + dt) - Ct
                                  # Clearness index
                                  Kt <- 1 - self$model$transform$alpha - self$model$transform$beta * exp(-exp(Yt))
                                  Kt * dCt / dt + Ct * self$model$transform$beta * exp(Yt - exp(Yt)) * (self$mu_Y(Yt, t_now, B, dt) + 0.5 * (1 - exp(Yt)) * self$sigma_Y(t_now, B)^2)
                                },
                                #' @description
                                #' Solar radiation diffusion \eqn{\sigma_(R)}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for \eqn{R_t}.
                                sigma_R = function(Rt, t_now, B = 1){
                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Diffusion for Rt process
                                  self$Ct(t_now) * self$model$transform$beta * exp(Yt - exp(Yt)) * self$sigma_Y(t_now, B)
                                },
                                #' @description
                                #' Compute the integral for expectation \eqn{\mu_(t,T)}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_expectation = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma <- private$..integral_expectation(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Compute the integral for variance \eqn{\sigma^2_(t,T)}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_variance = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma2 <- private$..integral_variance(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Integral mixture drift of both component of \eqn{Y_t}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Mixture expected value for both component of \eqn{Y_t}.
                                e_mix_drift = function(t_now, t_hor, df_date){
                                  if (missing(df_date)){
                                    # Create a sequence of dates
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Adjust parameters for Q-measure
                                  df_params <- self$model$NM_model$coefficients
                                  # Expected drift
                                  df_params$e_mu_B <- df_params$mu1 * df_params$p1 + df_params$mu2 * df_params$p2
                                  # Expected diffusion drift
                                  df_params$e_sigma_B <- df_params$sd1 * df_params$p1 + df_params$sd2 * df_params$p2
                                  # Combine the datasets
                                  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
                                  # !! Maybe can be removed !! tocheck
                                  df <- df[order(df$n),]
                                  # Compute the integral
                                  df$int <- private$..integral_expectation(df$n, df$N, df$tau)
                                  # Index for last day
                                  nrows <- nrow(df)
                                  # Total drift from time t up to T-1 for expectation under P
                                  mix_drift_mu <- df[-nrows,]$int * df[-nrows,]$e_mu_B
                                  # Conditional drift from time T-1 up to T for expectation under P
                                  mix_drift_mu_1 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu1
                                  mix_drift_mu_2 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu2
                                  # Total drift from time t up to T-1 for expectation under Q
                                  mix_drift_sigma_B <- df[-nrows,]$int * df[-nrows,]$e_sigma_B
                                  # Conditional drift from time T-1 up to T for expectation under Q
                                  mix_drift_sd_1 <- sum(mix_drift_sigma_B) + df[nrows,]$int * df[nrows,]$sd1
                                  mix_drift_sd_2 <- sum(mix_drift_sigma_B) + df[nrows,]$int * df[nrows,]$sd2
                                  # Drift depending on lambda
                                  mix_drift_1 <- mix_drift_mu_1 + mix_drift_sd_1 * self$lambda
                                  mix_drift_2 <- mix_drift_mu_2 + mix_drift_sd_2 * self$lambda

                                  list(e_drift_1 = mix_drift_1, e_drift_2 = mix_drift_2,
                                       e_mu1 = mix_drift_mu_1, e_mu2 = mix_drift_mu_2,
                                       e_sd1 = mix_drift_sd_1, e_sd2 = mix_drift_sd_2)
                                },
                                #' @description
                                #' Integral mixture diffusion of both component of \eqn{Y_t}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Mixture expected value for both component of \eqn{Y_t}.
                                e_mix_diffusion = function(t_now, t_hor, df_date){
                                  if (missing(df_date)){
                                    # Create a sequence of dates
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Adjust parameters for Q-measure
                                  df_params <- self$model$NM_model$coefficients %>%
                                    dplyr::mutate(mu_diff = mu1 - mu2,
                                                  sigma_diff = sd1 - sd2,
                                                  sigma2_1 = sd1^2,
                                                  sigma2_2 = sd2^2)
                                  # Combine the datasets
                                  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
                                  df <- df[order(df$n),]
                                  # Compute the integral
                                  df$int_sigma2 <- private$..integral_variance(df$n, df$N, df$tau)
                                  # Index for last day
                                  nrows <- nrow(df)
                                  df_t <- df[-nrows,]
                                  df_T <- df[nrows,]

                                  # Drift variance P in t, T-1
                                  variance_drift_P <- sum(df_t$int_sigma2 * df_t$mu_diff^2 * df_t$p1 * df_t$p2)
                                  # Diffusion variance Q in t, T-1
                                  variance_drift_Q <- sum(df_t$int_sigma2 * df_t$sigma_diff^2 * df_t$p1 * df_t$p2)
                                  # Diffusion second moment in t, T-1
                                  variance_diffusion <- sum(df_t$int_sigma2 * (df_t$sigma2_1 *df_t$p1 + df_t$sigma2_2 * df_t$p2))
                                  # Total variance in t, T-1
                                  common_variance <- variance_drift_P + variance_diffusion + variance_drift_Q * self$lambda^2
                                  # Realized variance for each component between T-1 and T
                                  variance_1 <- common_variance + df_T$int_sigma2 * df_T$sigma2_1
                                  variance_2 <- common_variance + df_T$int_sigma2 * df_T$sigma2_2

                                  list(variance_1 = variance_1, variance_2 = variance_2,
                                       variance_drift_P = variance_drift_P, variance_drift_Q = variance_drift_Q,
                                       common_variance = common_variance + df_T$int_sigma2, variance_diffusion = variance_diffusion,
                                       last_1 = df_T$int_sigma2 * df_T$sigma2_1, last_2 = df_T$int_sigma2 * df_T$sigma2_2)
                                },
                                #' @description
                                #' Conditional expectation of \eqn{Y_T} given \eqn{Y_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Conditional mean for \eqn{Y_t}
                                M_Y = function(Rt, t_now, t_hor, df_date){
                                  # Create once the sequence of dates
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Time to maturity in days
                                  tau <- max(df_date$tau) - min(df_date$n)
                                  # Compute the drifts
                                  mix_drift <- self$e_mix_drift(df_date = df_date)
                                  # Forecast for Yt
                                  Y_forecast <- self$Yt_bar(t_hor) + (Yt - self$Yt_bar(t_now)) * exp(-self$theta * tau)
                                  # Expected value for each component
                                  c(up = Y_forecast + mix_drift[[1]], dw = Y_forecast + mix_drift[[2]], e_Yt = Y_forecast)
                                },
                                #' @description
                                #' Conditional variance of \eqn{Y_T}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Conditional variance for \eqn{Y_t}
                                S_Y = function(t_now, t_hor, df_date){
                                  # Create once the sequence of dates
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Compute the diffusions
                                  mix_diffusion <- self$e_mix_diffusion(df_date = df_date)
                                  # Compute the variances
                                  c(up = mix_diffusion[["variance_1"]], dw = mix_diffusion[["variance_2"]])
                                },
                                #' @description
                                #' Compute the conditional moments
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param quiet Logical, when `TRUE` there won't be displayed any messages.
                                Moments = function(t_now, t_hor, quiet = FALSE){
                                  # Convert in dates
                                  t_now <- as.Date(t_now)
                                  t_hor <- as.Date(t_hor)
                                  # Time to maturity
                                  tau <- as.numeric((t_hor - t_now)[1])
                                  # Dataset
                                  data <- self$model$data
                                  # Add strike price in terms of Y
                                  data$K_Y <- self$model$R_to_Y(data$GHI_bar, data$date)
                                  # Add realized payoff
                                  data$payoff <- (data$GHI_bar - data$GHI) * ifelse(data$GHI_bar > data$GHI, 1, 0)
                                  # Lagged values
                                  data$L1_Yt_bar <- lag(data$Yt_bar, tau)
                                  data$L1_Yt <- lag(data$Yt, tau)
                                  # Filter dataset
                                  data <- dplyr::filter(data, date >= min(t_hor) & date <= max(t_hor))
                                  nstep <- 1
                                  # 1) Compute the fixed parameters
                                  start_time <- Sys.time()
                                  if(!quiet) message("Fitting Y forecast (", nstep, "/3) ...", appendLF = FALSE)
                                  # Forecasted Y
                                  data$Y_forecast <- data$Yt_bar + (data$L1_Yt - data$L1_Yt_bar) * exp(-self$theta * tau)
                                  end_time <- Sys.time()
                                  if(!quiet) message("Done in ", difftime(end_time, start_time, units = "secs"), " secs \r")
                                  nstep <- nstep + 1
                                  # 2) Drifts
                                  start_time <- Sys.time()
                                  if(!quiet) message("Fitting Mixture drift (", nstep, "/3) ...", appendLF = FALSE)
                                  e_mix_drift <- purrr::map2(t_now, t_hor, ~unlist(self$e_mix_drift(.x, .y)))
                                  data$e_mu1 <- purrr::map_dbl(e_mix_drift, ~.x["e_mu1"])
                                  data$e_mu2 <- purrr::map_dbl(e_mix_drift, ~.x["e_mu2"])
                                  data$e_sd1 <- purrr::map_dbl(e_mix_drift, ~.x["e_sd1"])
                                  data$e_sd2 <- purrr::map_dbl(e_mix_drift, ~.x["e_sd2"])
                                  end_time <- Sys.time()
                                  if(!quiet) message("Done in ", difftime(end_time, start_time, units = "secs"), " secs. \r")
                                  nstep <- nstep + 1
                                  # 3) Diffusion
                                  start_time <- Sys.time()
                                  if(!quiet) message("Fitting Mixture diffusions (", nstep, "/3) ...", appendLF = FALSE)
                                  e_mix_diffusion <- purrr::map2(t_now, t_hor, ~unlist(self$e_mix_diffusion(.x, .y)))
                                  data$variance_drift_P <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_drift_P"])
                                  data$variance_drift_Q <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_drift_Q"])
                                  data$variance_diffusion <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_diffusion"])
                                  data$last_1 <- purrr::map_dbl(e_mix_diffusion, ~.x["last_1"])
                                  data$last_2 <- purrr::map_dbl(e_mix_diffusion, ~.x["last_2"])
                                  end_time <- Sys.time()
                                  if(!quiet) message("Done in ", difftime(end_time, start_time, units = "secs"), " secs. \r")
                                  # Add bounds parameters
                                  data$alpha <- self$model$transform$alpha
                                  data$beta <- self$model$transform$beta
                                  data$date_from <- data$date - tau
                                  # Select only relevant variables
                                  dplyr::select(data, t_now = "date_from", t_hor = "date",
                                                Y_forecast, e_mu1, e_mu2, e_sd1, e_sd2,
                                                variance_drift_P, variance_drift_Q, variance_diffusion,
                                                last_1, last_2, p1, alpha, beta, payoff, Ct, Yt_bar, GHI, GHI_bar, K_Y)
                                },
                                #' @description
                                #' Conditional density of \eqn{Y_T} given \eqn{Y_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional density for \eqn{Y_T}.
                                pdf_Y = function(Rt, t_now, t_hor, B){
                                  # Create once the sequence of dates
                                  df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  # Expected values Y at horizon
                                  M_Y <- unlist(self$M_Y(Rt, t_now, t_hor, df_date = df_date))
                                  # Variance of Y at horizon
                                  S_Y <- unlist(self$S_Y(t_now, t_hor, df_date = df_date))
                                  # Mixture probability at horizon time
                                  p <- self$model$NM_model$prob$predict(t_hor)
                                  # Densities
                                  if (missing(B)) {
                                    function(x){p * dnorm(x, M_Y[1], sqrt(S_Y[1])) + (1 - p) * dnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  } else if (B == 1) {
                                    function(x) {dnorm(x, M_Y[1], sqrt(S_Y[1]))}
                                  } else {
                                    function(x) {dnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  }
                                },
                                #' @description
                                #' Conditional distribution of \eqn{Y_T} given \eqn{Y_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture distribution, otherwise the component distribution non weighted.
                                #' @return Conditional distribution for \eqn{Y_T}.
                                cdf_Y = function(Rt, t_now, t_hor, B){
                                  # Create once the sequence of dates
                                  df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  # Expected values Y at horizon
                                  M_Y <- unlist(self$M_Y(Rt, t_now, t_hor, df_date = df_date))
                                  # Variance of Y at horizon
                                  S_Y <- unlist(self$S_Y(t_now, t_hor, df_date = df_date))
                                  # Mixture probability at horizon time
                                  p <- self$model$NM_model$prob$predict(t_hor)
                                  if (missing(B)) {
                                    function(x){p * pnorm(x, M_Y[1], sqrt(S_Y[1])) + (1 - p) * pnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  } else if (B == 1) {
                                    function(x){pnorm(x, M_Y[1], sqrt(S_Y[1]))}
                                  } else {
                                    function(x){pnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  }
                                },
                                #' @description
                                #' Conditional density of \eqn{R_T} given \eqn{R_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional density for \eqn{R_T}
                                pdf_R = function(Rt, t_now, t_hor, B){
                                  C_T <- self$Ct(t_hor)
                                  pdf_Y <- self$pdf_Y(Rt, t_now, t_hor, B)
                                  function(x){
                                    z_x <- (1 - x/C_T - self$model$transform$alpha) / self$model$transform$beta
                                    u_x <- suppressWarnings(log(-log(z_x)))
                                    num <- pdf_Y(u_x)
                                    den <- C_T * self$model$transform$beta * log(z_x^z_x)
                                    probs <- -num/den
                                    probs[is.nan(probs)] <- 0
                                    return(probs)
                                  }
                                },
                                #' @description
                                #' Conditional distribution of \eqn{R_T} given \eqn{R_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture distribution, otherwise the component distribution non weighted.
                                #' @return Conditional distribution for \eqn{R_T}
                                cdf_R = function(Rt, t_now, t_hor, B){
                                  C_T <- self$Ct(t_hor)
                                  cdf_Y <- self$cdf_Y(Rt, t_now, t_hor, B)
                                  function(x){
                                    z_x <- (1 - x/C_T - self$model$transform$alpha) / self$model$transform$beta
                                    u_x <- suppressWarnings(log(-log(z_x)))
                                    probs <- cdf_Y(u_x)
                                    probs[is.nan(probs)] <- 0
                                    return(probs)
                                  }
                                },
                                #' @description
                                #' Conditional expected value of \eqn{R_T} given \eqn{R_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @param moment Integer, scalar. Moment order. The default is 1, i.e. the expectation.
                                #' @return Conditional moment for solar radiation
                                e_GHI = function(Rt, t_now, t_hor, B, moment = 1){
                                  pdf_Y <- self$pdf_Y(Rt, t_now, t_hor, B)
                                  C_T <- self$Ct(t_hor)
                                  Rt <- function(y) (C_T * (1 - self$model$transform$alpha - self$model$transform$beta * exp(-exp(y))))^moment * pdf_Y(y)
                                  moment_Rt <- integrate(Rt, lower = -Inf, upper = Inf)$value
                                  return(moment_Rt)
                                },
                                #' @description
                                #' Conditional variance value of \eqn{R_T} given \eqn{R_t}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional variance for \eqn{R_T}
                                v_GHI = function(Rt, t_now, t_hor, B){
                                  self$e_GHI(Rt, t_now, t_hor, B, 2) - self$e_GHI(Rt, t_now, t_hor, B, 1)^2
                                },
                                #' @description
                                #' Method print for `radiationModel` object.
                                print = function(){
                                  self$model$print()
                                  cat("MRP: ", self$lambda, "\n")
                                  cat("Measure: ", self$measure, "\n")
                                  cat("Version: ", private$version, "\n")
                                }
                              ),
                              # ====================================================================================================== #
                              #                                             Private slots
                              # ====================================================================================================== #
                              private = list(
                                version = "1.0.0",
                                ..model = NA,
                                ..measure = "P",
                                ..lambda = 0,
                                ..k1 = 1,
                                ..k2 = 1,
                                ..integral_variance = NA,
                                ..integral_expectation = NA,
                                seasonal_variance = NA
                              ),
                              # ====================================================================================================== #
                              #                                             Active slots
                              # ====================================================================================================== #
                              active = list(
                                #' @field model An object of the class `solarModel`.
                                model = function(){
                                  private$..model
                                },
                                #' @field measure Character, reference probability measure used.
                                measure = function(){
                                  private$..measure
                                },
                                #' @field lambda Numeric, market risk premium used.
                                lambda = function(){
                                  private$..lambda
                                }
                              )
)
