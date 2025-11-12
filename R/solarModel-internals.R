#' Compute the covariance
#'
#' @rdname solarModel_covariance
#' @name solarModel_covariance
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_covariance <- function(t_now, mom_t, mom_T, GARCH, NM_model, theta = 0, tol = 0.01){

  t_cond <- mom_t$date
  t_hor <- mom_T$date
  print(paste0("Covariance ", t_cond, " - ", t_hor))
  mom_T_1 <- solarMoments_path(mom_T, GARCH, NM_model, theta = theta, B = 1, t_cond)
  mom_T_0 <- solarMoments_path(mom_T, GARCH, NM_model, theta = theta, B = 0, t_cond)

  # Mean YT when Bt = 1, BT = 1
  mom_T$M_YT_11 <- mom_T_1$M_Y1
  # Mean Yt when Bt = 1, BT = 0
  mom_T$M_YT_10 <- mom_T_1$M_Y0
  # Mean YT when Bt = 0, BT = 1
  mom_T$M_YT_01 <- mom_T_0$M_Y1
  # Mean Yt when Bt = 0, BT = 0
  mom_T$M_YT_00 <- mom_T_0$M_Y0

  # Variance YT when Bt = 1, BT = 1
  mom_T$S_YT_11 <- mom_T_1$S_Y1
  # Variance YT when Bt = 1, BT = 0
  mom_T$S_YT_10 <- mom_T_1$S_Y0
  # Variance YT when Bt = 0, BT = 1
  mom_T$S_YT_01 <- mom_T_0$S_Y1
  # Variance YT when Bt = 0, BT = 0
  mom_T$S_YT_00 <- mom_T_0$S_Y0

  # Mean Yt when Bt = 1, BT = 1 equal to mean Yt when Bt = 1, BT = 0
  mom_T$M_Yt_11 <- mom_T$M_Yt_10 <- mom_t$M_Y1
  # Mean Yt when Bt = 1, BT = 1 equal to mean Yt when Bt = 1, BT = 0
  mom_T$M_Yt_01 <- mom_T$M_Yt_00 <- mom_t$M_Y0
  # Variance Yt when Bt = 1, BT = 1 equal to Variance Yt when Bt = 1, BT = 0
  mom_T$S_Yt_11 <- mom_T$S_Yt_10 <- mom_t$S_Y1
  # Variance Yt when Bt = 0, BT = 1 equal to Variance Yt when Bt = 0, BT = 0
  mom_T$S_Yt_01 <- mom_T$S_Yt_00 <- mom_t$S_Y0

  # Covariance
  data_psi_t <- mom_t$psi_j[[1]]
  data_psi_T_1 <- dplyr::filter(mom_T_1$psi_j[[1]], date <= mom_t$date)
  data_psi_T_0 <- dplyr::filter(mom_T_0$psi_j[[1]], date <= mom_t$date)
  t <- nrow(data_psi_t)
  data_psi_t$psi_t_1 <- data_psi_t$psi_j * data_psi_t$sigma_bar^2 * c(data_psi_t$variance[-t], data_psi_t$sd1[t]^2)
  data_psi_t$psi_t_0 <- data_psi_t$psi_j * data_psi_t$sigma_bar^2 * c(data_psi_t$variance[-t], data_psi_t$sd2[t]^2)

  # Covariance Bt = 0, BT = 1 equivalent to Bt = 0, BT = 1
  cv_given_1 <- sum(data_psi_t$psi_t_1 * data_psi_T_1$psi_j)
  # sum(data_psi_t$psi_t_1 * data_psi_T_0$psi_j * data_psi_T_1$mean)
  cv_given_0 <- sum(data_psi_t$psi_t_0 * data_psi_T_0$psi_j)
  # sum(data_psi_t$psi_t_0 * data_psi_T_0$psi_j * data_psi_T_1$mean)

  # Conditional pdf Bt = 1, BT = 1
  mu_11 <- c(mom_T$M_YT_11, mom_T$M_Yt_11)
  sigma11 <- matrix(c(mom_T$S_YT_11^2, cv_given_1, cv_given_1, mom_T$S_Yt_11^2), 2, 2, byrow = TRUE)
  pdf_11 <- function(x) mvtnorm::dmvnorm(x, mean = mu_11, sigma = sigma11)
  # Conditional pdf Bt = 0, BT = 1
  mu_01 <- c(mom_T$M_YT_01, mom_T$M_Yt_01)
  sigma01 <- matrix(c(mom_T$S_YT_01^2, cv_given_0, cv_given_0, mom_T$S_Yt_01^2), 2, 2, byrow = TRUE)
  pdf_01 <- function(x) mvtnorm::dmvnorm(x, mean = mu_01, sigma = sigma01)
  # Conditional pdf Bt = 1, BT = 0
  mu_10 <- c(mom_T$M_YT_10, mom_T$M_Yt_10)
  sigma10 <- matrix(c(mom_T$S_YT_10^2, cv_given_1, cv_given_1, mom_T$S_Yt_10^2), 2, 2, byrow = TRUE)
  pdf_10 <- function(x) mvtnorm::dmvnorm(x, mean = mu_10, sigma = sigma10)
  # Conditional pdf Bt = 0, BT = 0
  mu_00 <- c(mom_T$M_YT_00, mom_T$M_Yt_00)
  sigma00 <- matrix(c(mom_T$S_YT_00^2, cv_given_0, cv_given_0, mom_T$S_Yt_00^2), 2, 2, byrow = TRUE)
  pdf_00 <- function(x) mvtnorm::dmvnorm(x, mean = mu_00, sigma = sigma00)
  # Marginals
  pdf_Yt <- function(x) dmixnorm(x, c(mom_t$M_Y1, mom_t$M_Y0), c(mom_t$S_Y1, mom_t$S_Y0), c(mom_t$p1, 1-mom_t$p1))
  pdf_YT <- function(x) dmixnorm(x, c(mom_T$M_Y1, mom_T$M_Y0), c(mom_T$S_Y1, mom_T$S_Y0), c(mom_T$p1, 1-mom_T$p1))
  # Joint mixture pdf
  joint_pdf_YtT <- function(x) mom_T$p1 * mom_t$p1 * pdf_11(x) + mom_T$p1 * (1-mom_t$p1) * pdf_01(x) + (1-mom_T$p1) * mom_t$p1 * pdf_10(x) + (1-mom_T$p1) * (1-mom_t$p1) * pdf_00(x)
  # Store the pdfs
  mom_T$pdf <- list(pdf = list(Yt = pdf_Yt, YT = pdf_YT, joint = joint_pdf_YtT))
  # Compute covariance
  #mom_T$e_RT <- integrate(function(x) transform$iRY(x, mom_T$Ct)*pdf_YT(x), lower = -Inf, upper = Inf)$value
  #mom_T$e_Rt <- integrate(function(x) transform$iRY(x, mom_t$Ct)*pdf_Yt(x), lower = -Inf, upper = Inf)$value
  #mom_T$e_Rt_RT <- cubature::hcubature(function(x) transform$iRY(x[1], mom_T$Ct)*transform$iRY(x[2], mom_t$Ct)*joint_pdf_YtT(x),
  #                   tol = tol, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))$integral
  mom_T$t_cond <- mom_t$date
  # Covariance
  #mom_T$cv_Rt_RT <- mom_T$e_Rt_RT - mom_T$e_Rt * mom_T$e_RT
  # mom_T <- dplyr::select(mom_T, date, t_cond, Year:beta, pdf:cv_Rt_RT)
  return(mom_T)
}

#' Produce a forecast from a solarModel object
#'
#' @param lambda Numeric scalar, Sugeno parameter.
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' moments <- model$moments$conditional[14,]
#' object <- solarModel_predict(model, moments, ci = 0.01)
#' object
#'
#' @rdname solarModel_predict
#' @name solarModel_predict
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_predict <- function(model, moments, lambda = 0, ci = 0.01){
  # Moments
  df_n <- moments
  # Create datasets for the moments
  comb <- dplyr::tibble(mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), probs = c(df_n$p1, 1-df_n$p1))
  if (lambda == 0) {
    # Mixture pdf
    pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
    pdf_Yt_up <- function(x) dmixnorm(x, comb$mean[1], comb$sd[1], 1)
    pdf_Yt_dw <- function(x) dmixnorm(x, comb$mean[2], comb$sd[2], 1)
    # Mixture cdf
    cdf_Yt <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
    cdf_Yt_up <- function(x) pmixnorm(x, comb$mean[1], comb$sd[1], 1)
    cdf_Yt_dw <- function(x) pmixnorm(x, comb$mean[2], comb$sd[2], 1)
  } else {
    # Distribution and Density of Yt
    pdf_Yt_ <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
    pdf_Yt_up_ <- function(x) dmixnorm(x, comb$mean[1], comb$sd[1], 1)
    pdf_Yt_dw_ <- function(x) dmixnorm(x, comb$mean[2], comb$sd[2], 1)
    # Mixture cdf
    cdf_Yt_ <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
    cdf_Yt_up_ <- function(x) pmixnorm(x, comb$mean[1], comb$sd[1], 1)
    cdf_Yt_dw_ <- function(x) pmixnorm(x, comb$mean[2], comb$sd[2], 1)
    # Sugeno pdf
    pdf_Yt <- dsugeno(pdf_Yt_, cdf_Yt_, lambda)
    pdf_Yt_up <- dsugeno(pdf_Yt_up_, cdf_Yt_up_, lambda)
    pdf_Yt_dw <- dsugeno(pdf_Yt_dw_, cdf_Yt_dw_, lambda)
    # Sugeno cdf
    cdf_Yt <- psugeno(cdf_Yt_, lambda)
    cdf_Yt_up <- psugeno(cdf_Yt_up_, lambda)
    cdf_Yt_dw <- psugeno(cdf_Yt_dw_, lambda)
  }
  # Expected value of Rt^q
  e_Rt_q <- function(q = 1, pdf_Yt) integrate(function(x) model$transform$iRY(x, df_n$Ct)^q * pdf_Yt(x), lower = -Inf, upper = Inf)$value

  # Expected values
  # Mixture
  df_n$e_Rt <- e_Rt_q(1, pdf_Yt)
  # Sunny state
  df_n$e_Rt_up <- e_Rt_q(1, pdf_Yt_up)
  # Cloudy state
  df_n$e_Rt_dw <- e_Rt_q(1, pdf_Yt_dw)

  # Variances
  # Mixture
  df_n$v_Rt <- e_Rt_q(2, pdf_Yt) - df_n$e_Rt^2
  # Sunny state
  df_n$v_Rt_up <- e_Rt_q(2, pdf_Yt_up) - df_n$e_Rt_up^2
  # Cloudy state
  df_n$v_Rt_dw <- e_Rt_q(2, pdf_Yt_dw)- df_n$e_Rt_dw^2

  # Confidence interval
  # Mixture
  df_n$ci_mix_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt, link = model$spec$transform$link)
  df_n$ci_mix_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt, link = model$spec$transform$link)
  # Sunny state
  df_n$ci_up_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_up, link = model$spec$transform$link)
  df_n$ci_up_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_up, link = model$spec$transform$link)
  # Cloudy state
  df_n$ci_dw_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_dw, link = model$spec$transform$link)
  df_n$ci_dw_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_dw, link = model$spec$transform$link)

  # Number of points for the grid
  n_points <- 100
  # Compute bounds for GHI
  lower_Rt = df_n$Ct*model$transform$bounds("Kt")[1]
  upper_Rt = df_n$Ct*model$transform$bounds("Kt")[2]
  # Grid for PDF plot
  grid_x <- seq(lower_Rt, upper_Rt, length.out = n_points+2)[-c(1,n_points+2)]
  grid <- dplyr::tibble(x = grid_x)
  # Density GHI
  pdf_Rt <- function(x, pdf_Yt) dsolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt,
                                          link = model$spec$transform$link)
  # Density points (Mixture)
  grid$pdf_Rt_mix <- pdf_Rt(grid$x, pdf_Yt)
  # Density points (Mixture, up)
  grid$pdf_Rt_mix_up <- pdf_Rt(grid$x, pdf_Yt_up) * df_n$p1
  # Density points (Mixture, dw)
  grid$pdf_Rt_mix_dw <- pdf_Rt(grid$x, pdf_Yt_dw) * (1 - df_n$p1)

  # Add the value of the points for plotting
  # Expected value and variance
  # Mixture
  df_n$pdf_e_Rt <- pdf_Rt(df_n$e_Rt, pdf_Yt)
  # Confidence intervals
  # Mixture
  df_n$pdf_ci_mix_lo <- pdf_Rt(df_n$ci_mix_lo, pdf_Yt)
  df_n$pdf_ci_mix_hi <- pdf_Rt(df_n$ci_mix_hi, pdf_Yt)
  # Sunny state
  df_n$pdf_ci_up_lo <- pdf_Rt(df_n$ci_up_lo, pdf_Yt)
  df_n$pdf_ci_up_hi <- pdf_Rt(df_n$ci_up_hi, pdf_Yt)
  # Cloudy state
  df_n$pdf_ci_dw_lo <- pdf_Rt(df_n$ci_dw_lo, pdf_Yt)
  df_n$pdf_ci_dw_hi <- pdf_Rt(df_n$ci_dw_hi, pdf_Yt)

  # Realized GHI
  Rt <- dplyr::filter(model$data, date == df_n$date)$GHI
  df_n$Rt <- ifelse(purrr::is_empty(Rt), NA_integer_, Rt)

  structure(
    list(
      grid = grid,
      df_n = df_n,
      ci = ci
    ),
    class = c("solarModelForecast", "list")
  )
}

#' Iterate the forecast on multiple dates
#'
#' @rdname solarModel_forecast
#' @name solarModel_forecast
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_forecast <- function(model, moments, ci = 0.1, lambda = 0){
  fun <- function(df_n){
    safe_forecaster <- purrr::safely(solarModel_predict)
    smf <- safe_forecaster(model, moments = df_n, ci = ci, lambda = lambda)$result
    if (is.null(smf)){
      return(NULL)
    }
    return(smf[-1][[1]])
  }
  out <- purrr::map_df(1:nrow(moments), ~fun(moments[.x,]))
  return(out)
}

#' Match solarModel parameters in vector form
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' vec_params <- c(theta = 1, alpha1 = 10)
#' solarModel_match_params(vec_params, model$coefficients)
#'
#' @rdname solarModel_match_params
#' @name solarModel_match_params
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_match_params <- function(vec_params, params){
  # List of names of the original parameters
  names_params <- purrr::map(params, ~names(.x))
  # Names of the parameters in vector form
  names_vec_params <- names(vec_params)
  for(i in 1:length(vec_params)){
    # Parameter to use
    target_param <- vec_params[i]
    # Parameter name
    names_target_param <- names(target_param)
    # Check the parameter inside the list of parameters
    idx_list <- purrr::map_lgl(names_params, ~sum(.x %in% names_target_param) == 1)
    if (purrr::is_empty(which(idx_list))) {
      cli::cli_alert_warning(paste0('Parameter: "', names_target_param, '" not found in the model specification. Ignored!'))
      next
    }
    idx_list <- which(idx_list)
    # Position of the parameter inside the list
    idx_position <- names(params[[idx_list]]) == names_target_param
    # Update the parameter
    params[[idx_list]][idx_position] <- target_param[[1]]
  }
  return(params)
}


