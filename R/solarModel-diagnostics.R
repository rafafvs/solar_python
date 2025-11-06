#' Distribution test
#'
#' Evaluate a Kolmogorov-Smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution.
#'
#' @param model An object of the class `solarModel`
#' @param test Character, null hypothesis for the residuals distribution `"gm"` for Gaussian mixture and `"norm"` for normality.
#' @inheritParams ks_test
#' @param type Type of test.
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_distribution(model)
#'
#' @rdname solarModel_test_distribution
#' @name solarModel_test_distribution
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_distribution <- function(model, H0 = c("gm", "norm"), ci = 0.05, min_quantile = 0.025, max_quantile = 0.985,
                                         type = c("train", "test", "full")){
  # Null on the distribution
  test_H0 <- match.arg(H0, choices = c("gm", "norm"))
  # Type of data
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract data
  data <- model$data
  # Filter to exclude or not test data
  if (type == "train") {
    data <- data[data$isTrain & data$weights != 0,]
  } else if (type == "test") {
    data <- data[!data$isTrain,]
  }
  # Monthly tests
  tests <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(data, Month == nmonth)$u_tilde
    if (test_H0 == "norm"){
      # Normal CDF
      cdf_Yt <- function(x) pnorm(x, mean(x), sd(x))
    } else {
      # Gaussian mixture model
      gm <- model$NM_model$model[[nmonth]]
      # Mixture CDF
      cdf_Yt <- function(x) pmixnorm(x, gm$means, gm$sd, gm$p)
    }
    # Distribution test
    tests[[nmonth]] <- ks_test(x, cdf_Yt, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile)
  }
  tests <- dplyr::bind_cols(Month = 1:12, test_H0 = test_H0, data = type, dplyr::bind_rows(tests))
  tests$result <- ifelse(tests$H0 == "Non-Rejected", "Passed", "Not-passed")
  return(tests)
}

#' Autocorrelation test
#'
#' Evaluate the autocorrelation in the components of a `solarModel` object.
#'
#' @param model An object of the class `solarModel`
#' @param lag.max Numeric, scalar. Maximum lag to consider for the test.
#' @param ci Numeric, scalar. Minimum p-value to consider the test `"passed"`.
#' @param method Character, type of test. Can be `"bg"` for Breusch-Godfrey, `"bp"` for Box-pierce and `"lb"` for BLjung-Box.
#' @param train_data Logical, when `TRUE` only train data will be used to evaluate the test, otherwise all the available sample.
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_autocorr(model, method = "lb")
#'
#' @rdname solarModel_test_autocorr
#' @name solarModel_test_autocorr
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_autocorr <- function(model, lag.max = 3, ci = 0.05, method = c("bg", "bp", "lb"), type = c("train", "test", "full")){
  # Type of test
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract data
  data <- model$data
  if (type == "train") {
    data <- data[data$isTrain & data$weights != 0,]
  } else if (type == "test") {
    data <- data[!data$isTrain,]
  }
  # Standardized mixture residuals
  data <- dplyr::left_join(data, dplyr::select(model$NM_model$moments, Month, mean, variance), by = c("Month"))
  data$u <- data$u_tilde
  data$u_tilde <- (data$u - data$mean) / sqrt(data$variance)
  # Squared random variables
  data$Yt_tilde2 <- data$Yt_tilde^2
  data$eps2 <- data$eps^2
  data$eps_tilde2 <- data$eps_tilde^2
  data$u_2 <- data$u^2
  data$u_tilde2 <- data$u_tilde^2
  # Autocorrelation test Breusch-Godfrey
  bg_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    formula <- as.formula(paste0(target, "~ 1"))
    test <- lmtest::bgtest(formula, order = lag.max, data = data)
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Autocorrelation test Ljung-Box
  lb_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    test <- Box.test(data[[target]], lag = lag.max, type = "Ljung-Box")
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Autocorrelation test Box-pierce
  bp_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    test <- Box.test(data[[target]], lag = lag.max, type = "Box-Pierce")
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Residuals
  targets <- c("Yt_tilde", "eps", "eps_tilde", "u", "u_tilde")
  expected <- c("Rejected", "Not-rejected", "Not-rejected", "Not-rejected", "Not-rejected")
  # Squared residuals
  targets <- c(targets, "Yt_tilde2", "eps2", "eps_tilde2", "u_2", "u_tilde2")
  expected <- c(expected, "Rejected", "Rejected", "Rejected", "Not-rejected", "Not-rejected")
  names(expected) <- targets

  tests <- list()
  method <- match.arg(method, choices = c("bg", "bp", "lb"))
  for(target in targets){
    if (method == "bg"){
      tests[[target]] <- bg_test(target, data, lag.max, ci, expected = expected[target])
    } else if (method == "bp") {
      tests[[target]] <- bp_test(target, data, lag.max, ci, expected = expected[target])
    } else if (method == "lb") {
      tests[[target]] <- lb_test(target, data, lag.max, ci, expected = expected[target])
    }
  }
  tests <- dplyr::bind_rows(tests)
  return(tests)
}

#' Autocorrelation and distribution tests
#'
#' Evaluate a Kolmogorov-Smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution and a Breush-pagan or Box-pierce
#' test on the residuals.
#' @param lags Numeric vector. Lags on which perform the autocorrelation tests. Can be more than one.
#' @inheritParams solarModel_test_distribution
#' @inheritParams solarModel_test_autocorr
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_tests(model, train_data = TRUE)
#' @rdname solarModel_tests
#' @name solarModel_tests
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_tests <- function(model, lags = c(7), ci = 0.05, min_quantile = 0.025, max_quantile = 0.985,
                             method = "bg", type = c("train", "test", "full")){
  # Tests for absence of autocorrelation
  autocorr_tests <- list()
  for(i in 1:length(lags)){
    autocorr_tests[[i]] <- solarModel_test_autocorr(model, lag.max = lags[i], ci = ci, method, type)
  }
  # Standard names
  names(autocorr_tests) <- paste0("lag_", lags)

  # Tests for normal distribution
  normality_test <- solarModel_test_distribution(model, H0 = "norm", ci, min_quantile, max_quantile, type)
  # Tests for gaussian mixture distribution
  mixture_test <- solarModel_test_distribution(model, H0 = "gm", ci, min_quantile, max_quantile, type)
  structure(
    list(
      autocorr = autocorr_tests,
      normality = normality_test,
      mixture = mixture_test
    )
  )
}

#' Probability Integral Transform
#'
#' @rdname solarModel_test_PIT
#' @name solarModel_test_PIT
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_PIT <- function(model, ci = 0.05, type = c("train", "test", "full")){
  # Type of computation
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  Rt <- model$data$GHI
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
    Rt <- Rt[model$data$isTrain]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
    Rt <- Rt[!model$data$isTrain]
  }
  u <- c()
  link <- model$spec$transform$link
  for(i in 1:nrow(moments)){
    # Moments
    df_n <- moments[i,]
    # Distribution of Yt
    cdf_Y <- function(x) pmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
    # Grades of Rt
    u[i] <- psolarGHI(Rt[i], df_n$Ct, df_n$alpha, df_n$beta, cdf_Y, link = link)
  }
  moments$u <- u
  structure(
    list(
      data = dplyr::select(moments, date, Year, Month, Day, u),
      test = dplyr::bind_cols(link = model$spec$transform$link,
                              ks_test(u, function(x) punif(x), ci = ci))
    )
  )
}

#' Compute the Log-predictive density of a solarModel
#'
#' @rdname solarModel_test_LPD
#' @name solarModel_test_LPD
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_LPD <- function(model, type = c("train", "test", "full")){
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  LPD <- model$logLik(moments, target = "GHI")
  LPD <- LPD[!is.infinite(LPD)]

  dplyr::tibble(
    type = type,
    link = model$spec$transform$link,
    LPD =  mean(LPD, na.rm = TRUE)
  )
}

#' Compute metrics to test forecasts
#'
#' @rdname solarModel_test_forecast
#' @name solarModel_test_forecast
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_forecast <- function(model, ci = 0.1, type = c("train", "test", "full") ){
  # Type of dataset
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  # Compute forecasts
  forecast <- solarModel_forecast(model, moments, ci = ci)

  # Evaluate the goodness of fit with different metrics
  errors <- forecast$Rt - forecast$e_Rt
  # SSE: sum of the squared errors
  SSE <- sum(errors^2)
  # RMSE: root mean squared error
  RMSE <- sqrt(mean(errors^2))
  # MAE: mean absolute error
  MAE <- mean(abs(errors))
  # MAPE: mean absolute percentage error
  MAPE <- mean(abs(errors) / forecast$Rt * 100)
  # Violations of Upper VaR
  viol_ci_hi <- mean(forecast$Rt > forecast$ci_mix_hi)
  # Violations of Lower VaR
  viol_ci_lo <- mean(forecast$Rt < forecast$ci_mix_lo)

  dplyr::tibble(type = type,
                SSE = SSE,
                RMSE = RMSE,
                MAE = MAE,
                MAPE = MAPE,
                ci = ci*2,
                viol_ci_hi = viol_ci_hi,
                viol_ci_lo = viol_ci_lo)
}


#' Compute metrics to test option pricing
#'
#' @rdname solarModel_test_pricing
#' @name solarModel_test_pricing
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_pricing <- function(model, type = c("train", "test", "full"), control = control_solarOption()){
  # Type of dataset
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  # Option prices
  put_prices <- solarOption_model(model, moments, put = TRUE, control_options = control)
  call_prices <- solarOption_model(model, moments, put = FALSE, control_options = control)
  # PUT / Call Stats
  put_call_stats <- function(price, exercise, Gamma){
    errors <- price - Gamma
    dplyr::tibble(
      SSE = sum(errors^2),
      premium = sum(price),
      payoff = sum(Gamma),
      diff = premium - payoff,
      cPt = mean(price[exercise != 0]),
      cGamma = mean(Gamma[exercise != 0]),
      cdiff = cGamma - cPt
    )
  }
  # PUT Stats
  put_stats <- put_call_stats(put_prices$payoff$premium, put_prices$payoff$exercise, put_prices$payoff$payoff)
  # Call Stats
  call_stats <- put_call_stats(call_prices$payoff$premium, call_prices$payoff$exercise, call_prices$payoff$payoff)

  # SoRadIDX Put stats
  # Average payoff
  Pt_hist = solarOption_historical(model, put = TRUE, control_options = control)$payoff_year$premium
  # Average premium
  Pt <- put_prices$payoff_year$premium
  # Realized payoff
  Gamma_p <- put_prices$payoff_year$payoff
  # Stats
  put_idx_stats <- bind_cols(Pt_hist = Pt_hist, Pt = Pt, Gamma = Gamma_p) %>%
    mutate(diff_Gamma_Pt = Gamma - Pt, diff_Gamma_Pt_hist = Gamma - Pt_hist, diff_Pt_hist_Pt = Pt_hist-Pt)

  # SoRadIDX Call stats
  Ct_hist = solarOption_historical(model, put = FALSE, control_options = control)$payoff_year$premium
  # Average premium
  Ct <- call_prices$payoff_year$premium
  # Realized payoff
  Gamma_c <- call_prices$payoff_year$payoff
  # Stats
  call_idx_stats <- bind_cols(Ct_hist = Ct_hist, Ct = Ct, Gamma = Gamma_c)%>%
    mutate(diff_Gamma_Ct = Gamma - Ct, diff_Gamma_Ct_hist = Gamma - Ct_hist, diff_Ct_hist_Ct = Ct_hist-Ct)

  structure(
    list(
      put = bind_cols(type = type, K = control$K, n = nrow(moments), put_stats),
      call = bind_cols(type = type, K = control$K, n = nrow(moments), call_stats),
      put_idx = bind_cols(type = type, K = control$K, n = nrow(moments), put_idx_stats),
      call_idx = bind_cols(type = type, K = control$K, n = nrow(moments), call_idx_stats)
    )
  )
}

