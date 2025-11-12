#' Compute the conditional moments
#'
#' @param data Slot `data` of a `solarModel` object. See the function \code{\link{solarModel}} for details.
#' @param theta Numeric, shift parameter for the Gaussian mixture residuals.
#' @param control_model An object with the class `solarModel_spec`. See the function \code{\link{solarModel_spec}} for details.
#'
#' @keywords solarMoments
#' @details Version 1.0.0.
#' @rdname solarMoments_conditional
#' @name solarMoments_conditional
#' @export
solarMoments_conditional = function(data, theta = 0, control_model){
  if (!control_model$garch_variance) {
    data$sigma <- 1
  }
  # Compute conditional moments
  data <- dplyr::mutate(data,
                        # Change of measure
                        theta = theta,
                        # Conditional expectation Yt
                        e_Yt = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond,
                        # Conditional std. deviation Yt
                        sd_Yt = sigma * sigma_bar * sigma_uncond,
                        # Conditional moments Yt (state up)
                        M_Y1 = e_Yt + sd_Yt * (mu1 + sd1 * theta),
                        S_Y1 = sd_Yt * sd1,
                        # Conditional moments Yt (state dw)
                        M_Y0 = e_Yt + sd_Yt * (mu2 + sd2 * theta),
                        S_Y0 = sd_Yt * sd2)
  # Store only relevant variables
  data <- dplyr::select(data, date, Year, Month, Day, e_Yt, sd_Yt, M_Y1, S_Y1, M_Y0, S_Y0, theta)
  return(data)
}

#' Compute the unconditional moments
#'
#' @param data Slot `data` of a `solarModel` object. See the function \code{\link{solarModel}} for details.
#' @param ARMA Slot `ARMA` of a `solarModel` object.
#' @param GARCH Slot `GARCH` of a `solarModel` object.
#' @param theta Numeric, shift parameter for the Gaussian mixture residuals.
#'
#' @keywords solarMoments
#' @details Version 1.0.0.
#' @rdname solarMoments_unconditional
#' @name solarMoments_unconditional
#' @export
solarMoments_unconditional = function(data, ARMA, GARCH, theta = 0){
  # ARMA variance
  arma_variance <- sqrt(ARMA$variance(1000)[1000])
  # GARCH long term std. deviation
  GARCH_vol <- GARCH$omega / (1 - sum(GARCH$coefficients[-1]))

  # Compute unconditional moments
  data <- dplyr::mutate(data,
                        # Change of measure
                        theta = theta,
                        # Unconditional expectation Yt
                        e_Yt = ARMA$intercept + Yt_bar + Yt_tilde_uncond,
                        # Unconditional std. deviation Yt
                        sd_Yt = sigma_bar * sigma_uncond * GARCH_vol * arma_variance,
                        # Unconditional moments Yt (state up)
                        M_Y1 = e_Yt + sd_Yt * (mu1 + sd1 * theta),
                        S_Y1 = sd_Yt * sd1,
                        # Unconditional moments Yt (state dw)
                        M_Y0 = e_Yt + sd_Yt * (mu2 + sd2 * theta),
                        S_Y0 = sd_Yt * sd2)
  # Store only relevant variables
  data <- dplyr::select(data, date, Year, Month, Day, e_Yt, sd_Yt, M_Y1, S_Y1, M_Y0, S_Y0, theta)
  return(data)
}

#' Compute the generic conditional moments of a solarModel object
#'
#' @param t_now Date for today.
#' @param t_hor Horizon date.
#' @param data Slot `data` of a `solarModel` object. See the function \code{\link{solarModel}} for details.
#' @param ARMA Slot `ARMA` of a `solarModel` object.
#' @param GARCH Slot `GARCH` of a `solarModel` object.
#' @param NM_model Slot `NM_model` of a `solarModel` object.
#' @param transform Slot `transform` of a `solarModel` object.
#' @param theta Numeric, vector of shift parameters for the Gaussian mixture residuals.
#'
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' t_now = "2019-07-11"
#' t_hor = "2019-10-19"
#' data = model$data
#' ARMA = model$ARMA
#' GARCH  = model$GARCH
#' NM_model = model$NM_model
#' transform = model$transform
#' theta = 0
#' solarMoments(t_now, t_hor, data, ARMA, GARCH, NM_model, transform, theta = 0, quiet = FALSE)
#' filter(model$moments$conditional, date == t_hor)
#' filter(model$moments$unconditional, date == t_hor)
#'
#' t_seq <- seq.Date(as.Date("2013-01-01"), as.Date("2013-12-31"), 1)
#' mom <- purrr::map_df(t_seq, ~solarr::solarMoments(.x-150, .x, data, ARMA, GARCH, NM_model, transform, theta = 0, quiet = FALSE))
#' solarOption_model(model, mom, control_options = control_solarOption(nyears = c(2005, 2024)))
#' solarOption_historical(model, control_options = control_solarOption(nyears = c(2005, 2024)))
#'
#' @keywords solarMoments
#' @details Version 1.0.0.
#' @rdname solarMoments
#' @name solarMoments
#' @export
solarMoments <- function(t_now, t_hor, data, ARMA, GARCH, NM_model, transform, theta = 0, quiet = FALSE){
  if(!quiet) cli::cli_alert_info(paste0("Forecast: ", t_hor, " given ", t_now))
  # Conditioning date
  t_now <- as.Date(t_now)
  # Horizon date
  t_hor <- as.Date(t_hor)
  # Maximum order of AR / GARCH
  lag_max <- max(c(ARMA$order, GARCH$order))
  # Filter data between (t_now - lag_max + 1) and t_hor
  data <- dplyr::filter(data, date >= (t_now - lag_max + 1) & date <= t_hor)
  # Unknown data till maturity
  df_tT <- data[-c(1:lag_max),]
  # Known data at time t_now (used as vector only in state-space forecast)
  df_t <- data[c(1:lag_max),]
  # Must be reversed ordered, i.e. most recent (t, t-1, t-2, ...)
  df_t <- df_t[order(df_t$date, decreasing = TRUE),]
  # Forecasting horizon
  h <- nrow(df_tT)
  # Store alpha and beta
  df_tT$alpha <- transform$alpha
  df_tT$beta <- transform$beta
  # *******************************************************************************
  #  0)  Forecast mean and variance of Yt_tilde
  # *******************************************************************************
  # Companion matrix
  A <- ARMA$A
  # Residuals vector for mean
  b <- matrix(ARMA$b, ncol = 1)
  # Residuals matrix for variance
  B <- b %*% t(b)
  # Extract first component
  e1 <- matrix(rep(1, length(b)), ncol = 1)
  e1[-1] <- 0
  # Conditioning values
  Y0 <- df_t$Yt_tilde[1:ARMA$order[1]]
  eps0 <- c()
  if (ARMA$order[2] > 0){
    eps0 <- df_t$eps[1:ARMA$order[2]]
  }
  # State vector
  Xt <- c(Y0, eps0)
  # Initialize the variables
  df_tT$psi_j <- df_tT$psi2_j <- NA
  df_tT$intercept <- 0
  j <- 1
  for(j in 1:h){
    # Pre compute the powers
    A_pow_j <- pow_matrix(A, j)
    A_pow_hj <- pow_matrix(A, h - j)
    # Compute weights for expectations
    df_tT$psi_j[j] <- t(e1) %*% (A_pow_hj %*% b)
    # Intercept contribution
    df_tT$intercept[j] <- t(e1) %*% (A_pow_j %*% (e1 * (ARMA$intercept)))
    # Variance
    df_tT$psi2_j[j] <- t(e1) %*% (A_pow_hj %*% B %*% t(A_pow_hj)) %*% e1
    # Forecasted value
    df_tT$Yt_tilde_hat[j] <- ARMA$intercept + sum(df_tT$intercept[1:j]) + t(e1) %*% A_pow_j %*% Xt
  }
  # *******************************************************************************
  #  2)  Forecast seasonal variance
  # *******************************************************************************
  df_tT$sigma_bar <- df_tT$sigma_bar * df_tT$sigma_uncond
  solarMoments_path(df_tT, GARCH, NM_model, theta = theta)
}

#' Condition the moments for a specific Bernoulli realization at a time t_cond
#'
#' @inheritParams solarMoments
#' @param B conditioning value for the Bernoulli state at time thor
#' @param t_cond conditioning date
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' GARCH  = model$GARCH
#' NM_model = model$NM_model
#' # Compute the moments
#' moments <- model$Moments("2012-01-01", "2012-03-05")
#'
#' # Condition the moments on a realization of B
#' moments
#' solarOption_model(model, moments)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = 1, "2012-02-03")
#' solarOption_model(model, mom_cond)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = 1, "2012-03-01")
#' solarOption_model(model, mom_cond)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = 1, "2012-03-04")
#' solarOption_model(model, mom_cond)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = 0, "2012-03-04")
#' solarOption_model(model, mom_cond)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = c(0, 0), c("2012-03-03", "2012-03-04"))
#' solarOption_model(model, mom_cond)
#'
#' mom_cond <- solarMoments_path(moments, GARCH, NM_model, theta = 0, B = c(1, 1), c("2012-03-03", "2012-03-04"))
#' solarOption_model(model, mom_cond)
#'
#' @keywords solarMoments
#' @details Version 1.0.0.
#' @rdname solarMoments_path
#' @name solarMoments_path
#' @export
solarMoments_path <- function(moments, GARCH, NM_model, theta = 0, B = 1, t_cond){
  if (!is.list(moments$psi_j)){
    df_tT <- moments
  } else {
    df_tT <- moments$psi_j[[1]]
  }
  # Forecasting horizon
  h <- nrow(df_tT)
  # Add mixture parameters
  df_tT <- dplyr::select(df_tT, -mu1, -mu2, -sd1, -sd2, -p1, -p2)
  df_tT <- dplyr::left_join(df_tT, NM_model$coefficients, by = "Month")
  if (any(theta != 0)){
    df_tT$mu1 <- df_tT$mu1 + df_tT$sd1 * theta
    df_tT$mu2 <- df_tT$mu2 + df_tT$sd2 * theta
  }
  # Compute mixture moments
  df_tT$mean <- (df_tT$mu1 * df_tT$p1) + (df_tT$mu2 * df_tT$p2)
  # Second moment
  df_tT$m2 <- (df_tT$mu1^2 + df_tT$sd1^2) * df_tT$p1 + (df_tT$mu2^2 + df_tT$sd2^2) * df_tT$p2
  # Fourth moment
  df_tT$m4 <- (3 * df_tT$sd1^4 + 6 * df_tT$mu1^2 * df_tT$sd1^2 + df_tT$mu1^4) * df_tT$p1 + (3 * df_tT$sd2^4 + 6 * df_tT$mu2^2 * df_tT$sd2^2 + df_tT$mu2^4) * df_tT$p2
  # Variance
  df_tT$variance <- df_tT$m2 - df_tT$mean^2

  # *******************************************************************************
  #  1) Conditioning variable
  # *******************************************************************************
  if (!missing(t_cond)) {
    t_cond <- as.Date(t_cond)
    for(t in 1:length(t_cond)){
      # Second moment
      df_tT$m2 <- ifelse(df_tT$date == t_cond[t], ifelse(B[t] == 1, df_tT$mu1^2 + df_tT$sd1^2, df_tT$mu2^2 + df_tT$sd2^2), df_tT$m2)
      # Fourth moment
      df_tT$m4 <- ifelse(df_tT$date == t_cond[t], ifelse(B[t] == 1, 3 * df_tT$sd1^4 + 6 * df_tT$mu1^2 * df_tT$sd1^2 + df_tT$mu1^4, 3 * df_tT$sd2^4 + 6 * df_tT$mu2^2 * df_tT$sd2^2 + df_tT$mu2^4), df_tT$m4)
      # Expectation
      df_tT$mean <- ifelse(df_tT$date == t_cond[t], ifelse(B[t] == 1, df_tT$mu1, df_tT$mu2), df_tT$mean)
      # Variance
      df_tT$variance <- ifelse(df_tT$date == t_cond[t], ifelse(B[t] == 1, df_tT$sd1^2, df_tT$sd2^2), df_tT$variance)
    }
  }
  # *******************************************************************************
  #  1) Forecast GARCH moments
  # *******************************************************************************
  if (h == 1) {
    df_tT$e_sigma2 <- df_tT$sigma^2
    df_tT$e_sigma4 <- df_tT$sigma^4
  } else {
    # Second moment GARCH variance (exact)
    df_tT$e_sigma2 <- e_sigma2_h_mix(h - 1, GARCH$omega, GARCH$alpha, GARCH$beta, e_x2 = df_tT$m2, df_tT$sigma[1]^2)
    # Second moment of GARCH std. dev (exact)
    df_tT$e_sigma4 <- e_sigma4_h_mix(h - 1, GARCH$omega, GARCH$alpha, GARCH$beta, e_x2 = df_tT$m2, e_x4 = df_tT$m4, df_tT$sigma[1]^4)
  }
  # Variance of GARCH (exact)
  df_tT$v_sigma2 <- df_tT$e_sigma4 - df_tT$e_sigma2^2
  #df_tT$v_sigma2
  # First moment of GARCH std.dev (approximated)
  df_tT$e_sigma1 <-  df_tT$e_sigma2^(1/2) - (1/8) * df_tT$v_sigma2 / sqrt(df_tT$e_sigma2)^3
  # Variance of GARCH (exact)
  df_tT$v_sigma <- df_tT$v_sigma2 / (4*df_tT$e_sigma2)

  df_tT$cv_psi_ij <- 0
  if (h > 1) {
    for(i in 2:nrow(df_tT)){
      cv_psi_ij <- cov_sigma2_h_mix(i, 1, GARCH$omega, GARCH$alpha, GARCH$beta, e_x2 = df_tT$m2, e_x4 = df_tT$m4, df_tT$sigma[1]^2)
      cv_psi_ij <- cv_psi_ij / (4 * sqrt(df_tT$e_sigma2[i] * df_tT$e_sigma2[2:i]))
      cv_psi_ij <- cv_psi_ij * (df_tT$psi_j[i] * df_tT$psi_j[2:i])
      df_tT$cv_psi_ij[i] <- sum(cv_psi_ij)
    }
  }
  # *******************************************************************************
  #  3)  Compute the series of psi
  # *******************************************************************************
  # Compute weights for expectations
  df_tT$psi_y <- df_tT$psi_j * df_tT$sigma_bar * df_tT$e_sigma1 * df_tT$mean
  # Variance
  df_tT$psi2_y <- df_tT$psi2_j * df_tT$sigma_bar^2 * (df_tT$e_sigma2 * df_tT$m2 - df_tT$e_sigma1^2 * df_tT$mean^2)
  # Forecasted mean value
  df_tT$e_Yt_tilde <- df_tT$Yt_tilde_hat + cumsum(df_tT$psi_y)
  # Mean value
  df_tT$e_Yt <- df_tT$e_Yt_tilde + df_tT$Yt_tilde_uncond + df_tT$Yt_bar
  # Forecasted total variance
  df_tT$Sigma2 <- cumsum(df_tT$psi2_y)
  # Store data to compute covariances
  data_psi <-  df_tT %>%
    dplyr::select(date, Month, e_sigma1, e_sigma2, e_sigma4, v_sigma, v_sigma2, Yt_tilde_hat, Yt_tilde_uncond, Yt_bar,
                  psi_j, psi_y, psi2_j, psi2_y,  mu1, mu2, sd1, sd2, p1, p2, mean, m2, m4, variance, sigma_bar, sigma)
  # Last value
  df_T <- tail(df_tT, 1)
  # *******************************************************************************
  #  4)  Approximate the multinomial mixture with a 2-component GM
  # *******************************************************************************
  # Next step parameters
  # Approximated mixture means
  M_Y <- df_T$Yt_tilde_hat + df_T$Yt_tilde_uncond + df_T$Yt_bar + df_T$sigma_bar * df_T$e_sigma1 * c(df_T$mu1, df_T$mu2)
  # Approximated mixture variances
  S2_Y <- c(df_T$sd1^2, df_T$sd2^2) * df_T$sigma_bar^2 * df_T$e_sigma2 + df_T$v_sigma * df_T$sigma_bar^2 * c(df_T$mu1^2, df_T$mu2^2)
  if (h > 1){
    # Add conditional covariances
    rho2_U <- sum((df_tT$cv_psi_ij * df_tT$mean)[-h]) + df_T$cv_psi_ij * c(df_T$mu1, df_T$mu2)
    # Approximated mixture means
    M_Y <- M_Y + sum(df_tT$psi_y[-h])
    # Approximated mixture variances
    S2_Y <- S2_Y + sum(df_tT$psi2_y[-h]) + 2*rho2_U
    df_tT$Sigma2 <- df_tT$Sigma2 + cumsum((df_tT$cv_psi_ij * df_tT$mean))
  }
  S_Y <- sqrt(S2_Y)

  # GM_moments(M_Y, S_Y, c(df_T$p1, 1 - df_T$p1))
  # Target second moment constraint
  # m2_target <- df_T$Sigma2
  # Target third moment constraint
  # m3_target <- df_T$Omega
  # Calibrate the variance
  #p <- c(df_T$p1, 1 - df_T$p1)
  # Central moments
  #delta_k <- M_Y - sum(M_Y*p)
  # Solve the variances to match the second and third moment
  #D <- 3 * p[1] * p[2] * (delta_k[2] - delta_k[1])
  #D1 <- (m2_target - p[1] * delta_k[1]^2 - p[2] * delta_k[2]^2) * 3 * delta_k[2] * p[2] - p[2] * (m3_target - p[1] * delta_k[1]^3 - p[2] * delta_k[2]^3)
  #D0 <- p[1] * (m3_target - p[1] * delta_k[1]^3 - p[2] * delta_k[2]^3)  - (m2_target - p[1] * delta_k[1]^2 - p[2] * delta_k[2]^2) * 3 * delta_k[1] * p[1]
  #S2_Y_star <- c(D1/D, D0/D)
  t_hor <- max(df_T$date)
  moments <- dplyr::filter(moments, date == t_hor)
  dplyr::tibble(
    date = t_hor,
    h = h,
    Year = lubridate::year(t_hor),
    Month = lubridate::month(t_hor),
    Day = lubridate::day(t_hor),
    e_Yt = df_T$e_Yt,
    sd_Yt = ifelse(h == 1, df_T$sigma_bar * df_T$sigma, sqrt(df_T$Sigma2)),
    M_Y1 = M_Y[1],
    S_Y1 = S_Y[1],
    M_Y0 = M_Y[2],
    S_Y0 = S_Y[2],
    psi_j = list(data_psi),
    p1 = moments$p1[1],
    Ct = moments$Ct[1],
    theta = tail(theta, n = 1),
    GHI_bar = moments$GHI_bar[1],
    alpha = moments$alpha[1],
    beta = moments$beta[1]
  )
}

#' Calibrate theta to match a certain level of solar radiation
#'
#' @param e_RT_target Numeric, vector of target expectation to match.
#' @inheritParams solarMoments
#'
#' @keywords solarModel
#' @details Version 1.0.0.
#' @rdname solarModel_calibrate_theta
#' @name solarModel_calibrate_theta
#' @export
solarModel_calibrate_theta <- function(model, moments, e_RT_target, quiet = FALSE){
  # Time to maturity
  h <- nrow(moments)
  if (length(e_RT_target) != h) stop("The length of `e_RT_target` is not equal to the number of rows of `moments`.")
  # Mixture model
  NM_model <- model$NM_model
  # GARCH model
  GARCH = model$GARCH
  # Loss function for theta_t
  loss_function <- function(theta_t, mom, theta, e_RT_target){
    # Moments under P(theta)
    mom_theta <- solarMoments_path(mom, NM_model = NM_model, GARCH = GARCH, theta = c(theta, theta_t))
    # Expected value of Rt
    e_RT <- solarModel_forecast(model, mom_Q, ci = 0.1)$e_Rt
    # Loss
    (e_RT - e_RT_target)^2
  }
  # Moments under P-theta
  mom_theta <- list()
  # Sequential calibration of theta
  theta <- c()
  for(t in 1:h){
    # Moments for day t
    mom_t <- moments[t,]
    # Optimal theta
    theta_t <- optim(0, loss_function, lower = -1, upper = 1,
                     method = "Brent", mom = mom_t, theta = theta, e_RT_target = e_RT_target[t])$par
    # Add the calibrated theta
    if(!quiet) cli::cli_alert_success(paste0("(", t, "/", h, ")", " Optimal theta_t: ", theta_t))
    theta <- c(theta, theta_t)
    # Moments under P-theta
    mom_theta[[t]] <- solarMoments_path(mom_t, NM_model = NM_model, GARCH = GARCH, theta = theta)
    mom_theta[[t]]$theta_t <- list(theta = theta)
  }
  # Create a unique dataset for the moments
  mom_theta <- dplyr::bind_rows(mom_theta)
  structure(
    list(
      moments = mom_theta,
      theta = theta
    )
  )
}

#' Calibrate the time series of theta to match a certain level of Option price
#'
#' @param P_target Numeric, vector of target prices to match.
#' @inheritParams solarMoments
#' @inheritParams solarOption_historical
#'
#' @keywords solarOption
#' @details Version 1.0.0.
#' @rdname solarOption_calibrate_theta
#' @name solarOption_calibrate_theta
#' @export
solarOption_calibrate_theta <- function(model, moments, P_target, put = TRUE, quiet = FALSE, control_options = control_solarOption()){
  # Time to maturity
  h <- nrow(moments)
  if (length(P_target) != h) stop("The length of `P_target` is not equal to the number of rows of `moments`.")
  # Mixture model
  NM_model <- model$NM_model
  # GARCH model
  GARCH = model$GARCH
  # Loss function for theta_t
  loss_function <- function(theta_t, mom, theta, P_target){
    # Moments under P(theta)
    mom_theta <- solarMoments_path(mom, NM_model = NM_model, GARCH = GARCH, theta = c(theta, theta_t))
    # Price under P(theta)
    price_P_theta <- solarOption_model(model, mom_theta, put = put, control_options = control_options)$payoff_year$premium
    # Loss
    (price_P_theta - P_target)^2
  }
  # Moments under P-theta
  mom_theta <- list()
  # Sequential calibration of theta
  theta <- c()
  for(t in 1:h){
    # Moments for day t
    mom_t <- moments[t,]
    # Optimal theta
    theta_t <- optim(0, loss_function, lower = -1, upper = 1,
                     method = "Brent", mom = mom_t, theta = theta, P_target = P_target[t])$par
    # Add the calibrated theta
    if(!quiet) cli::cli_alert_success(paste0("(", t, "/", h, ")", " Optimal theta_t: ", theta_t))
    theta <- c(theta, theta_t)
    # Moments under P-theta
    mom_theta[[t]] <- solarMoments_path(mom_t, NM_model = NM_model, GARCH = GARCH, theta = theta)
    mom_theta[[t]]$theta_t <- list(theta = theta)
  }
  # Create a unique dataset for the moments
  mom_theta <- dplyr::bind_rows(mom_theta)
  structure(
    list(
      moments = mom_theta,
      theta = theta
    )
  )
}


