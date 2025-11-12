#' Payoff of solar options on historical data
#'
#' @param model An object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param nmonths Numeric vector, months in which the payoff should be computed. Can vary from 1 (January) to 12 (December).
#' @param put Logical, when `TRUE`, the default, will be computed the price for a `put` contract, otherwise for a `call` contract.
#' @param control_options Named list with control parameters. See \code{\link{control_solarOption}} for more details.
#'
#' @return An object of the class `solarOptionPayoff`.
#' @examples
#' # Solar model
#' model <- solarModel$new(spec)
#' model$fit()
#' # Historical payoff (put)
#' solarOption_historical(model, put=TRUE)
#' # Historical payoff (call)
#' solarOption_historical(model, put=FALSE)
#'
#' @rdname solarOption_historical
#' @name solarOption_historical
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_historical <- function(model, nmonths = 1:12, put = TRUE, control_options = control_solarOption()){
  # Options control parameters
  K <- control_options$K
  # Target and seasonal mean
  target <- model$spec$target
  target_bar <- paste0(target, "_bar")
  # Complete data
  data <- dplyr::select(model$data, date, n, Year, Month, Day, tidyr::any_of(c(target, target_bar)), Ct)
  # Filter for control years
  data <- dplyr::filter(data, date >= control_options$from & date <= control_options$to)
  # Filter for selected months
  data <- dplyr::filter(data, Month %in% nmonths)
  # Rename target variable
  data$Rt <- data[[target]]
  # Compute strike price
  data$strike <- data[[target_bar]] * exp(K)
  # Select only relevant variables
  df_payoff <- dplyr::select(data, date, n, Year, Month, Day, Rt, strike, Ct)
  # Type of option
  df_payoff$side <- ifelse(put, "put", "call")
  # Daily option payoffs
  if (put) {
    # Exercise
    df_payoff$exercise <- ifelse(df_payoff$Rt < df_payoff$strike, 1, 0)
    # Realized payoff for a put
    df_payoff$payoff <- (df_payoff$strike - df_payoff$Rt) * df_payoff$exercise
    # Maximum payoff
    max_payoff <- df_payoff$strike - df_payoff$Ct * (1 - model$transform$alpha - model$transform$beta)
  } else {
    # Exercise
    df_payoff$exercise <- ifelse(df_payoff$Rt > df_payoff$strike, 1, 0)
    # Realized payoff for a call
    df_payoff$payoff <- (df_payoff$Rt - df_payoff$strike) * df_payoff$exercise
    # Maximum payoff
    max_payoff <- df_payoff$Ct * (1 - model$transform$alpha) - df_payoff$strike
  }
  # Aggregation of Historical Payoffs
  df_payoff$premium <- df_payoff$payoff
  # Compute maximum payoff for each day
  df_payoff$max_payoff <- max_payoff
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise, max_payoff)
  # Output structure
  solarOptionPayoff(df_payoff, control_options$leap_year)
}

#' Compute average premium on simulated Data
#'
#' @param scenario object with the class `solarModelScenario`. See the function \code{\link{solarModel_scenarios}} for details.
#' @param nsim number of simulation to use for computation.
#' @inheritParams solarOption_historical
#'
#' @param control_options control function, see \code{\link{control_solarOption}} for details.
#'
#' @return An object of the class `solarOptionPayoff`.
#' @examples
#' # Solar model
#' model <- solarModel$new(spec)
#' model$fit()
#' # Simulate scenarios
#' scenario <- solarScenario(model, from = "2011-01-01", to = "2012-01-01", by = "1 month", nsim = 10, seed = 3)
#'
#' solarOption_scenario(model, scenario)
#' solarOption_historical(model)
#' solarScenario_plot(scenario)
#'
#' @rdname solarOption_scenario
#' @name solarOption_scenario
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_scenario <- function(model, scenario, nmonths = 1:12, put = TRUE, nsim, control_options = control_solarOption()){
  # Control parameters
  nyears <- control_options$nyears
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- scenario$target
  target_bar <- paste0(target, "_bar")

  # Simulated daily Payoffs
  sim <- scenario$sim
  # Filter for control years
  df_payoff <- dplyr::filter(sim, date >= control_options$from & date <= control_options$to)
  # Filter for selected months
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonths)
  # Eventually reduce the number of simulations used
  if (!missing(nsim)) {
    nsim <- min(c(nsim, nrow(df_payoff$data[[1]])))
    df_payoff <- dplyr::mutate(df_payoff, data = purrr::map(data, ~.x[1:nsim,]))
  }
  # Add seasonal mean (strike)
  data <- dplyr::left_join(df_payoff, model$data[,c("date", target_bar, "Ct")], by = "date")
  df_payoff$strike <- data[[target_bar]]*exp(K)
  df_payoff$Ct <- data$Ct
  df_payoff <- tidyr::unnest(df_payoff, cols = "data")
  df_payoff$Rt <- df_payoff[[target]]
  # Compute simulated payoffs
  df_payoff <- df_payoff %>%
    dplyr::mutate(side = ifelse(put, "put", "call"),
                  exercise = dplyr::case_when(side == "put" & Rt <= strike ~ 1,
                                              side == "put" & Rt > strike ~ 0,
                                              side == "call" & Rt <= strike ~ 0,
                                              side == "call" & Rt > strike ~ 1),
                  payoff_sim = dplyr::case_when(side == "put" ~ (strike - Rt)*exercise,
                                                side == "call" ~ (Rt - strike)*exercise),
                  max_payoff = Ct * (1-model$transform$alpha - model$transform$beta),
    ) %>%
    dplyr::group_by(date) %>%
    dplyr::group_by(side, date, Year, Month, Day) %>%
    dplyr::reframe(
      Rt = mean(Rt),
      strike = mean(strike),
      premium = mean(payoff_sim),
      exercise = mean(exercise),
      max_payoff = mean(max_payoff)
    ) %>%
    dplyr::ungroup()
  # Add and compute realized payoff
  data <- dplyr::left_join(scenario$emp, model$data[,c("date", target_bar)], by = "date")
  scenario$emp$Rt <- scenario$emp[[target]]
  scenario$emp$strike <- data[[target_bar]]*exp(K)
  scenario$emp$side <- ifelse(put, "put", "call")
  scenario$emp <- dplyr::mutate(scenario$emp,
                                payoff = dplyr::case_when(
                                  side == "put" ~ (strike - Rt)*ifelse(Rt < strike, 1, 0),
                                  TRUE ~ (Rt - strike)*ifelse(Rt > strike, 1, 0)))

  df_payoff <- dplyr::left_join(df_payoff, scenario$emp[,c("date", "payoff")], by = "date")
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise, max_payoff)
  # Output structure
  solarOptionPayoff(df_payoff, control_options$leap_year)
}

#' Compute the premium given the moments
#'
#' @inheritParams solarOption_historical
#' @param moments Tibble containing the forecasted moments for different days ahead. See the function \code{\link{solarMoments}} for more details.
#' @param portfolio Optional, A list of objects of the class `solarOptionPortfolio`.
#' @param lambda Numeric, Sugeno parameter.
#' @param implvol Numeric, implied volatility.
#'
#' @examples
#' # Model
#' model = solarModel$new(spec)
#' model$fit()
#' # Pricing without portfolio
#' moments <- model$moments$unconditional
#' # Premium
#' premium_Vt <- solarOption_model(model, moments, put = TRUE)
#' # Pricing date
#' t_now <- as.Date("2021-12-31")
#' # Inception date
#' t_init <- as.Date("2022-01-01")
#' # Maturity date
#' t_hor <- as.Date("2022-12-31")
#' # solar options portfolio
#' portfolio <- SoRadPorfolio(model, t_now, t_init, t_hor)
#' # Moments
#' moments <- purrr::map_df(portfolio, ~model$Moments(t_now, .x$t_hor))
#' # Premium
#' premium_Vt <- solarOption_model(model, moments, portfolio, put = TRUE)
#' premium_Vt$payoff_year$premium
#' @rdname solarOption_model
#' @name solarOption_model
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_model <- function(model, moments, portfolio, nmonths = 1:12, lambda = 0, implvol = 1, put = TRUE, control_options = control_solarOption()){
  # Options control parameters
  if (missing(portfolio)) {
    moments <- dplyr::filter(moments, date >= control_options$from & date <= control_options$to)
    moments <- dplyr::filter(moments, Month %in% nmonths)
    moments$S_Y0 <- moments$S_Y0 * implvol
    moments$S_Y1 <- moments$S_Y1 * implvol
    moments$lambda <- lambda
    data <- purrr::map_df(1:nrow(moments), ~solarOption_pricing(moments[.x,], lambda = moments$lambda[.x],
                                                                put = put, link = model$spec$transform$link, control_options = control_options))

  } else {
    moments$S_Y0 <- moments$S_Y0 * implvol
    moments$S_Y1 <- moments$S_Y1 * implvol
    moments$lambda <- lambda
    data <- purrr::map_df(portfolio, ~solarOption_pricing(dplyr::filter(moments, date == .x$t_hor), .x,
                                                          lambda = dplyr::filter(moments, date == .x$t_hor)$lambda,
                                                          put = put, link = model$spec$transform$link, control_options = control_options))
  }
  # Add realized GHI
  data <- dplyr::left_join(data, dplyr::select(model$data, date, Rt = "GHI"), by = "date")
  # Compute realized payoff
  data <- dplyr::mutate(data,
                        payoff = dplyr::case_when(
                          side == "put" ~ (strike - Rt)*ifelse(Rt < strike, 1, 0),
                          TRUE ~ (Rt - strike)*ifelse(Rt > strike, 1, 0)))
  # Reorder variables
  solarOptionPayoff(data, control_options$leap_year)
}

#' Compute the price of a `solarOption`
#'
#' @param sorad An object of the class `solarOption`.
#' @inheritParams solarOption_model
#'
#' @examples
#' # Model
#' model = solarModel$new(spec)
#' model$fit()
#' moments <- filter(model$moments$conditional, Year == 2022)
#' # Pricing without contracts
#' solarOption_pricing(moments[1,])
#' # Pricing with contracts specification
#' sorad <- solarOption$new()
#' sorad$set_contract("2021-12-31", "2022-01-01", "2022-04-20", moments$GHI_bar[1])
#' solarOption_pricing(moments[1,], sorad)
#' solarOption_pricing(moments[1,], sorad, lambda = 0.02)
#' solarOption_pricing(moments[1,], sorad, lambda = -0.02)
#'
#' @rdname solarOption_pricing
#' @name solarOption_pricing
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_pricing <- function(moments, sorad, lambda = 0, put = TRUE, link = "invgumbel", control_options = control_solarOption()){
  # Moments
  df_n <- moments
  # Mixture moments
  comb <- dplyr::tibble(mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), probs = c(df_n$p1, 1-df_n$p1))
  # Density of Yt
  pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
  # Distribution of Yt
  cdf_Yt <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
  # Density and distribution of Rt
  if (lambda != 0){
    # Sugeno-distorted distributions
    pdf_R <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = link)
    cdf_R <- function(x) psolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, cdf_Yt, link = link)
    cdf_Rt <- psugeno(cdf_R, lambda = lambda)
    pdf_Rt <- dsugeno(pdf_R, cdf_R, lambda = lambda)
  } else {
    # Standard distribution of Rt
    pdf_Rt <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = link)
    cdf_Rt <- function(x) psolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, cdf_Yt, link = link)
  }
  # Strike price
  K <- df_n$GHI_bar * exp(control_options$K)
  t_hor <- df_n$date
  if (!missing(sorad)) {
    K <- sorad$strike * exp(control_options$K)
    t_hor <- sorad$t_hor
  }
  # G-function
  e_RT <- function(lower, upper) integrate(function(x) x * pdf_Rt(x),
                                           stop.on.error = FALSE, lower = lower, upper = upper)$value
  # Integration bounds
  bounds <- c(Rt_min = df_n$Ct*(1-df_n$alpha-df_n$beta), K = K, Rt_max = df_n$Ct*(1-df_n$alpha))

  # Option pricing
  if (put) {
    # Probability of exercise
    exercise <- cdf_Rt(K)
    # Option expected value
    premium <- K * exercise - e_RT(bounds[1], bounds[2])
    # Maximum payoff
    max_payoff <- K - bounds[1]
  } else {
    # Probability of exercise
    exercise <- 1 - cdf_Rt(K)
    # Option expected value
    premium <- e_RT(bounds[2], bounds[3]) - K * exercise
    # Maximum payoff
    max_payoff <- bounds[3] - K
  }

  # Select only relevant variables
  df_n <- dplyr::tibble(
    date = t_hor,
    Year = lubridate::year(date),
    Month = lubridate::month(date),
    Day = lubridate::day(date),
    side = ifelse(put, "put", "call"),
    premium = premium,
    exercise = exercise,
    strike = K,
    max_payoff = max_payoff)
  return(df_n)
}

#' Compute Choquet price for a Solar Option
#'
#' @inheritParams solarOption_model
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' moments <- model$moments$unconditional[1:365,]
#' lambda = 0
#' control_options = control_solarOption()
#' solarOption_choquet(model, moments[1:30,], lambda = 0.01)
#' solarOption_model(model, moments[1,])
#' solarOption_choquet(model, moments, lambda = 0.1, put = F)
#' @rdname solarOption_choquet
#' @name solarOption_choquet
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_choquet <- function(model, moments, portfolio, nmonths = 1:12, lambda = 0, implvol = 1, put = TRUE, control_options = control_solarOption()){
  # ASK price (seller)
  lam_neg <- purrr::map_dbl(lambda, ~sugeno_bounds(.x)$neg)
  premium_ask <- solarOption_model(model = model, moments = moments, portfolio, nmonths = nmonths,
                                   lambda = lam_neg, implvol = implvol, put = put, control_options = control_options)
  # BID price (buyer)
  lam_pos <- purrr::map_dbl(lambda, ~sugeno_bounds(.x)$pos)
  premium_bid <- solarOption_model(model = model, moments = moments, portfolio, nmonths = nmonths,
                                   lambda = lam_pos, implvol = implvol, put = put, control_options = control_options)

  df_payoff <- premium_ask$payoff %>%
    dplyr::select(-Rt, -strike, -payoff, -max_payoff) %>%
    dplyr::mutate(side = ifelse(put, "put", "call")) %>%
    dplyr::rename(premium_up = "premium", exercise_up = "exercise") %>%
    dplyr::left_join(premium_bid$payoff, by = c("side", "date", "Year", "Month", "Day")) %>%
    dplyr::rename(premium_dw = "premium", exercise_dw = "exercise") %>%
    dplyr::mutate(premium = (premium_up + premium_dw)/2, exercise = (exercise_up + exercise_dw)/2) %>%
    dplyr::select(side:Day, Rt, strike, premium_up, premium, premium_dw, exercise_up, exercise, exercise_dw, payoff, max_payoff)

  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   #premium = sum(premium),
                   premium = sum((premium_dw+premium_up)/2),
                   premium_dw = sum(premium_dw),
                   premium_up = sum(premium_up),
                   exercise_dw = mean(exercise_dw),
                   exercise = mean(exercise),
                   exercise_up = mean(exercise_up),
                   max_payoff = sum(max_payoff))

  # Aggregation by Month and Day
  df_month_day <- df_payoff %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = mean((premium_dw+premium_up)/2),
                   premium_dw = mean(premium_dw),
                   premium_up = mean(premium_up),
                   exercise_dw = mean(exercise_dw),
                   exercise = mean(exercise),
                   exercise_up = mean(exercise_up),
                   max_payoff = mean(max_payoff))

  # Aggregation by Month
  df_month <- df_year_month %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium_dw = mean(premium_dw),
                   premium = mean(premium),
                   premium_up = mean(premium_up),
                   exercise_dw = mean(exercise_dw),
                   exercise = mean(exercise),
                   exercise_up = mean(exercise_up),
                   max_payoff = mean(max_payoff))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium_dw = sum(premium_dw),
                   premium = sum(premium),
                   premium_up = sum(premium_up),
                   exercise_dw = mean(exercise_dw),
                   exercise = mean(exercise),
                   exercise_up = mean(exercise_up),
                   max_payoff = sum(max_payoff))
  # Structured output
  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionChoquet", "list")
  )
}

#' Compute the first fourth moments and variance of a SoRad
#' @inheritParams solarOption_historical
#' @param transform slot transform of a solarModel
#' @param lambda Sugeno parameter
#' @examples
#' # Solar model
#' model <- solarModel$new(spec)
#' model$fit()
#' moments <- model$moments$unconditional
#' solarOption_moments(moments[1:365,], model$transform)
#'
#' @rdname solarOption_moments
#' @name solarOption_moments
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_moments <- function(moments, transform, lambda = 0, put = TRUE, link = "invgumbel", control_options = control_solarOption()){
  # Moments
  df_tT <- moments
  # Selected payoff function
  if (put) {
    # Put function
    payoff <- function(x, K){(K - x) * ifelse(x > K, 0, 1)}
  } else {
    # Call function
    payoff <- function(x, K){(x - K) * ifelse(x < K, 0, 1)}
  }
  # Initialize the variables
  df_tT$e_Gamma1 <- df_tT$e_Gamma2 <- df_tT$e_Gamma3 <- df_tT$e_Gamma4 <- df_tT$v_Gamma <- df_tT$strike <- 0
  for(i in 1:nrow(df_tT)){
    # Moments at time i
    df_n <- df_tT[i,]
    # Density of Yt
    pdf_Yt <- function(x) dmixnorm(x, c(df_n$M_Y1, df_n$M_Y0), c(df_n$S_Y1, df_n$S_Y0), c(df_n$p1, 1-df_n$p1))
    # Density of Rt
    if (lambda != 0){
      # Sugeno-distorted distributions
      pdf_R <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = link)
      pdf_Rt <- dsugeno(pdf_R, cdf_R, lambda = lambda)
    } else {
      # Standard distribution of Rt
      pdf_Rt <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = link)
    }
    # Strike price
    df_tT$strike[i] <- df_n$GHI_bar * exp(control_options$K)
    # Integration bounds
    bounds <- c(Rt_min = df_n$Ct*(1-df_n$alpha-df_n$beta), Rt_max = df_n$Ct*(1-df_n$alpha))
    # First moment
    df_tT$e_Gamma1[i] <- integrate(function(x) payoff(x, df_tT$strike[i]) * pdf_Rt(x),
                                   stop.on.error = FALSE, lower = bounds[1], upper = bounds[2])$value
    # Second moment
    df_tT$e_Gamma2[i] <- integrate(function(x) payoff(x, df_tT$strike[i])^2 * pdf_Rt(x),
                                   stop.on.error = FALSE, lower = bounds[1], upper = bounds[2])$value
    # Third moment
    df_tT$e_Gamma3[i] <- integrate(function(x) payoff(x, df_tT$strike[i])^3 * pdf_Rt(x),
                                   stop.on.error = FALSE, lower = bounds[1], upper = bounds[2])$value
    # Fourth moment
    df_tT$e_Gamma4[i] <- integrate(function(x) payoff(x, df_tT$strike[i])^4 * pdf_Rt(x),
                                   stop.on.error = FALSE, lower = bounds[1], upper = bounds[2])$value
    # Variance
    df_tT$v_Gamma[i]<- df_tT$e_Gamma2[i] - df_tT$e_Gamma1[i]^2
  }
  # Select only relevant variables
  df_tT <- dplyr::select(df_tT, date:Day, Ct, strike, e_Gamma1:v_Gamma)
  return(df_tT)
}

#' Compute the covariance between two SoRad
#' @rdname solarOption_covariance
#' @name solarOption_covariance
#' @keywords solarOption
#' @note Version 1.0.0.
#' @noRd
solarOption_covariance <- function(t_now, mom_t, mom_T, model){

  t_cond <- mom_t$date
  t_hor <- mom_T$date
  print(paste0("Covariance ", t_cond, " - ", t_hor))

  mom_T_1 <- solarModel_moments_path(model, mom_T, Bj = 1, t_cond)
  mom_T_0 <- solarModel_moments_path(model, mom_T, Bj = 0, t_cond)

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
  cv_given_0 <- sum(data_psi_t$psi_t_0 * data_psi_T_0$psi_j)

  # Pdf
  mu_11 <- c(mom_T$M_YT_11, mom_T$M_Yt_11)
  sigma11 <- matrix(c(mom_T$S_YT_11^2, cv_given_1, cv_given_1, mom_T$S_Yt_11^2), 2, 2, byrow = TRUE)
  pdf_11 <- function(x) mvtnorm::dmvnorm(x, mean = mu_11, sigma = sigma11)
  # Pdf
  mu_01 <- c(mom_T$M_YT_01, mom_T$M_Yt_01)
  sigma01 <- matrix(c(mom_T$S_YT_01^2, cv_given_0, cv_given_0, mom_T$S_Yt_01^2), 2, 2, byrow = TRUE)
  pdf_01 <- function(x) mvtnorm::dmvnorm(x, mean = mu_01, sigma = sigma01)
  # Pdf
  mu_10 <- c(mom_T$M_YT_10, mom_T$M_Yt_10)
  sigma10 <- matrix(c(mom_T$S_YT_10^2, cv_given_1, cv_given_1, mom_T$S_Yt_10^2), 2, 2, byrow = TRUE)
  pdf_10 <- function(x) mvtnorm::dmvnorm(x, mean = mu_10, sigma = sigma10)
  # Pdf
  mu_00 <- c(mom_T$M_YT_00, mom_T$M_Yt_00)
  sigma00 <- matrix(c(mom_T$S_YT_00^2, cv_given_0, cv_given_0, mom_T$S_Yt_00^2), 2, 2, byrow = TRUE)
  pdf_00 <- function(x) mvtnorm::dmvnorm(x, mean = mu_00, sigma = sigma00)
  # Mixture pdf
  pdf_Yt <- function(x) mom_T$p1 * mom_t$p1 * pdf_11(x) + mom_T$p1 * (1-mom_t$p1) * pdf_01(x) + (1-mom_T$p1) * mom_t$p1 * pdf_10(x) + (1-mom_T$p1) * (1-mom_t$p1) * pdf_00(x)
  #pdf_Rt <- function(x) dmvsolarGHI(x, c(mom_T$Ct, mom_t$Ct), c(mom_T$alpha, mom_t$alpha), c(mom_T$beta, mom_t$beta), pdf_Yt)

  payoff <- function(x, K, Ct) (K - transform$iRY(x, Ct))*ifelse(transform$iRY(x, Ct) > K, 0, 1)
  mom_t$GHI_bar <- filter(model$data, date == mom_t$date)$GHI_bar
  mom_T$GHI_bar <- filter(model$data, date == mom_T$date)$GHI_bar

  e_Gamma_t <- integrate(function(x) payoff(x, mom_t$GHI_bar, mom_t$Ct) * dmixnorm(x, c(mom_t$M_Y1, mom_t$M_Y0), c(mom_t$S_Y1, mom_t$S_Y0), c(mom_t$p1, 1-mom_t$p1)),
                         stop.on.error = FALSE, lower = -Inf, upper = Inf)$value
  e_Gamma_T <- integrate(function(x) payoff(x, mom_T$GHI_bar, mom_T$Ct) * dmixnorm(x, c(mom_T$M_Y1, mom_T$M_Y0), c(mom_T$S_Y1, mom_T$S_Y0), c(mom_T$p1, 1-mom_T$p1)),
                         stop.on.error = FALSE, lower = -Inf, upper = Inf)$value
  e_Gamma_t_Gamma_T <- cubature::hcubature(function(x) payoff(x[2], mom_t$GHI_bar, mom_t$Ct)*payoff(x[1], mom_T$GHI_bar, mom_T$Ct)*pdf_Yt(x),
                                           tol = 0.05, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))$integral
  # Covariance
  e_Gamma_t_Gamma_T - e_Gamma_T * e_Gamma_t
}

#' Compute the Greeks of a Solar Option
#'
#' @examples
#' # Model fit
#' model <- solarModel$new(spec)
#' model$fit()
#' solarOption_greeks(model, model$Moments("2023-01-01", "2024-05-01"))
#'
#' @rdname solarOption_greeks
#' @name solarOption_greeks
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_greeks <- function(model, moments, put = TRUE, elasticities = FALSE){
  mom <- moments
  # Clearsky
  Ct <- mom$Ct
  # Strike
  K <- mom$GHI_bar
  # Bounds parameters
  alpha <- mom$alpha
  beta <- mom$beta
  # Mixture means
  mu1 <- mom$M_Y1
  mu0 <- mom$M_Y0
  # Mixture std. deviations
  sd1 <- mom$S_Y1
  sd0 <- mom$S_Y0
  # Mixture probabilities
  p <- mom$p1
  # Extract time of the year
  t <- number_of_day(mom$date)
  # Forecast moments change
  #mom_Th <- model$Moments(moments$psi_j[[1]]$date[1], moments$date+1, quiet = TRUE)
  dmu1_dT <- NA#mom_Th$M_Y1 - mom$M_Y1
  dmu0_dT <- NA#mom_Th$M_Y0 - mom$M_Y0
  dsd1_dT <- NA#mom_Th$S_Y1 - mom$S_Y1
  dsd0_dT <- NA#mom_Th$S_Y0 - mom$S_Y0
  dp_dT <- NA#mom_Th$p1 - mom$p1
  # Forecast moments change
  # mom_th <- model$Moments(moments$psi_j[[1]]$date[1]+1, moments$date, quiet = TRUE)
  dmu1_dt <- NA#mom_th$M_Y1 - mom$M_Y1
  dmu0_dt <- NA#mom_th$M_Y0 - mom$M_Y0
  dsd1_dt <- NA#mom_th$S_Y1 - mom$S_Y1
  dsd0_dt <- NA#mom_th$S_Y0 - mom$S_Y0
  # *************************************************************************
  # Normalization function
  eta <- function(x, alpha, beta, Ct) (1-alpha-x/Ct)/beta
  # Inverse Normalization function
  ieta <- function(x, alpha, beta, Ct) Ct*(1-alpha-beta*x)
  # Link function
  g <- function(x) model$transform$Q_H(x)
  # First derivative link function
  g_prime <- function(x) 1/(x*log(x))
  # Inverse link function
  ig <- function(y)  model$transform$F_H(y)
  # Distribution of Y
  F_Y <- function(x, mu1, mu0, sd1, sd0, p) p*pnorm(x, mu1, sd1) + (1-p)*pnorm(x, mu0, sd0)
  # Density of Y
  f_Y <- function(x, mu1, mu0, sd1, sd0, p) p*dnorm(x, mu1, sd1) + (1-p)*dnorm(x, mu0, sd0)
  # Integral function G_kq
  # G_kq <- function(x, q = 0, mu, sd) integrate(function(y) ((y-mu)^q)*ig(y)*f_Yk(y, mu, sd), upper = Inf, lower = x, subdivisions = 1000L)$value
  Gi_q <- function(x, q = 0, mu, sd) integrate(function(y) ((y-mu)^q)*ig(y)*dnorm(y, mu, sd), lower = -Inf, upper = x, subdivisions = 1000L)$value
  # Integral function G_q
  G_q <- function(x, q = 0, mu1, mu0, sd1, sd0, p) Gi_q(x, q, mu1, sd1) * p + (1-p) * Gi_q(x, q, mu0, sd0)
  # *************************************************************************
  # Standardized variable
  K_prime <- eta(K, alpha, beta, Ct)
  K_Y <- g(K_prime)
  # Probability of exercise
  Phi_Y <- F_Y(K_Y, mu1, mu0, sd1, sd0, p)
  Phi_Y1 <- pnorm(K_Y, mu1, sd1)
  Phi_Y0 <- pnorm(K_Y, mu0, sd0)
  # First derivative of the probability of exercise
  phi_Y <- f_Y(K_Y, mu1, mu0, sd1, sd0, p)
  phi_Y1 <- dnorm(K_Y, mu1, sd1)
  phi_Y0 <- dnorm(K_Y, mu0, sd0)
  # Integral G0
  G1_0 <- Gi_q(K_Y, 0, mu1, sd1)
  G0_0 <- Gi_q(K_Y, 0, mu0, sd0)
  G_0 <- p * G1_0 + (1-p) * G0_0
  # Integral G1
  G1_1 <- Gi_q(K_Y, 1, mu1, sd1)
  G0_1 <- Gi_q(K_Y, 1, mu0, sd0)
  # Integral G2
  G1_2 <- Gi_q(K_Y, 2, mu1, sd1)
  G0_2 <- Gi_q(K_Y, 2, mu0, sd0)
  # Exercise sets
  d1 <- (K_Y-mu1)/sd1
  d0 <- (K_Y-mu0)/sd0
  # Option price
  Pt <- (K - Ct*(1-alpha)) * Phi_Y + beta * Ct * G_0
  # *************************************************************************
  # Greeks
  # Derivative with respect to alpha
  d_Pt_d_alpha <- Ct * Phi_Y
  # Derivative with respect to beta
  d_Pt_d_beta <- Ct * G_0
  # Derivative with respect to strike
  d_Pt_d_strike <- Phi_Y
  # Derivative with respect to C_T
  d_Pt_d_CT <- -(1-alpha) * Phi_Y + beta * G_0
  # Derivative with respect to first mixture mean
  d_Pt_d_M_Y1 <- p * ( (Ct * (1-alpha) - K) * phi_Y1 + beta * Ct * G1_1 / sd1^2)
  # Derivative with respect to second mixture mean
  d_Pt_d_M_Y0 <- (1-p) * ( (Ct * (1 - alpha) - K) * phi_Y0 + beta * Ct * G0_1 / sd0^2)
  # Derivative with respect to first mixture std. deviation
  d_Pt_d_S_Y1 <-  (p/sd1) * ((Ct * (1 - alpha) - K) * phi_Y1 * d1 + beta * Ct * (-G1_0 + G1_2 / sd1^2))
  # Derivative with respect to first mixture std. deviation
  d_Pt_d_S_Y0 <-  ((1-p)/sd0) * ((Ct * (1 - alpha) - K) * phi_Y0 * d0 + beta*Ct * (-G0_0 + G0_2 / sd0^2))
  # Derivative with respect to Mixture probs
  d_Pt_d_pT <- (K - Ct * (1-alpha)) * (Phi_Y1 - Phi_Y0) + beta*Ct * (G1_0 - G0_0)
  # *************************************************************************
  # Derivative of the clear-sky
  c_ <- model$seasonal_model_Ct$coefficients
  # Derivative of extraterrestrial radiation
  d_Hon_dT <- model$seasonal_model_Ct$ssf$Hon(t, model$seasonal_model_Ct$lat, deriv = TRUE)
  # Periodicity
  omega <- (2*base::pi)/model$seasonal_model_Ct$period
  dC_dT <- c_[2]*d_Hon_dT + omega*(c_[3]*cos(omega*t) - c_[4]*sin(omega*t))
  # Derivative of KY with respect to T
  d_KY_dT  <- g_prime(K_prime) * K / beta * (1/Ct^2) * dC_dT
  d_FY1_dT <- (d_KY_dT - dmu1_dT  - d1  * dsd1_dT) * phi_Y1
  d_FY0_dT <- (d_KY_dT - dmu0_dT  - d0  * dsd0_dT) * phi_Y0
  d_FY_dT  <- dp_dT * (Phi_Y1 - Phi_Y0) + d_FY1_dT * p + d_FY0_dT * (1-p)
  # Derivative of G0 with respect to T
  d_G10_dT <- G1_1 / sd1^2 * dmu1_dT + (G1_2 / sd1^3 - G1_0 / sd1) * dsd1_dT
  d_G00_dT <- G0_1 / sd0^2 * dmu0_dT + (G0_2 / sd0^3 - G0_0 / sd0) * dsd0_dT
  d_G0_dT  <- K_prime * phi_Y * d_KY_dT + dp_dT * (G1_0 - G0_0) + d_G10_dT * p + d_G00_dT * (1-p)
  # Total derivative of Pt
  d_Pt_dT  <- dC_dT * (beta * G_0 - (1-alpha) * Phi_Y) + (K - Ct*(1-alpha)) * d_FY_dT + beta * Ct * d_G0_dT
  # *************************************************************************
  # Derivative of FY with respect to t
  d_FY1_dt <- phi_Y1 * (dmu0_dt  - d0  * dsd0_dt)
  d_FY0_dt <- phi_Y0 * (dmu1_dt  - d1  * dsd1_dt)
  d_FY_dt <- p * d_FY1_dt + (1 - p) * d_FY0_dt
  # Derivative of G0 with respect to t
  d_G10_dt <- G1_1 / sd1^2 * dmu1_dt + (G1_2 / sd1^3 - G1_0 / sd1) * dsd1_dt
  d_G00_dt <- G0_1 / sd0^2 * dmu0_dt + (G0_2 / sd0^3 - G0_0 / sd0) * dsd0_dt
  d_G0_dt  <- d_G10_dt * p + d_G00_dt * (1-p)
  # Total derivative
  d_Pt_dt <- (K - Ct * (1-alpha)) * d_FY_dt + beta * Ct * d_G0_dt

  if (elasticities) {
    dplyr::tibble(
      Pt = Pt,
      alpha = d_Pt_d_alpha * (alpha / Pt),
      beta = d_Pt_d_beta * (beta / Pt),
      strike = d_Pt_d_strike * (K / Pt),
      C_T = d_Pt_d_CT * (Ct / Pt),
      M_Y1 = d_Pt_d_M_Y1 * (mu1 / Pt),
      M_Y0 = d_Pt_d_M_Y0 * (mu0 / Pt),
      S_Y1 = d_Pt_d_S_Y1 * (sd1 / Pt),
      S_Y0 = d_Pt_d_S_Y0 * (sd0 / Pt),
      p = d_Pt_d_pT * (p / Pt),
      dT = d_Pt_dT,
      dt = d_Pt_dt
    )
  } else {
    dplyr::tibble(
      Pt = Pt,
      alpha = d_Pt_d_alpha,
      beta = d_Pt_d_beta,
      strike = d_Pt_d_strike,
      C_T = d_Pt_d_CT,
      M_Y1 = d_Pt_d_M_Y1,
      M_Y0 = d_Pt_d_M_Y0,
      S_Y1 = d_Pt_d_S_Y1,
      S_Y0 = d_Pt_d_S_Y0,
      p = d_Pt_d_pT,
      dT = d_Pt_dT,
      dt = d_Pt_dt
    )
  }
}

#' Compute the Greeks of a Solar Option index
#'
#' @examples
#' # Model fit
#' model <- solarModel$new(spec)
#' model$fit()
#'
#' t_now <- as.Date("2024-01-01")
#' t_hor <- as.Date("2024-12-31")
#' dates <- seq.Date(t_now+1, t_hor, 1)
#' mom <- purrr::map_df(dates, ~model$Moments(t_now, .x))
#' solarOption_index_greeks(model, mom, elasticities = TRUE)
#' mom_test <- mom
#' mom_test$beta <- mom_test$beta*1.01
#' solarOption_index_greeks(model, mom_test, elasticities = TRUE)
#'
#' @rdname solarOption_index_greeks
#' @name solarOption_index_greeks
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_index_greeks <- function(model, moments, elasticities = FALSE){

  if (elasticities) {
    # Elasticities
    greeks <- purrr::map_df(1:nrow(moments), ~solarOption_greeks(model, moments[.x,], elasticities = TRUE))
    # Weights proportion of price
    w_E <- greeks[,c("Pt")] / sum(greeks[,c("Pt")])
    # Normalize to allow summation
    idx_deriv <- which(!(colnames(greeks) %in% c("Pt", "date")))
    for(i in idx_deriv) {
      greeks[,i] <- greeks[,i] * w_E
    }
  } else {
    # Sensitivities
    greeks <- purrr::map_df(1:nrow(moments), ~solarOption_greeks(model, moments[.x,], elasticities = FALSE))
  }

  # Extract t_now
  t_now <- head(moments$date, 1) - 1
  # Extract t_hor
  t_hor <- tail(moments$date, 1)
  # Time to maturity in days
  time_to_maturity <- nrow(moments)
  # Summarise the greeks
  greeks <- dplyr::summarise_all(greeks, sum)

  dplyr::bind_cols(
    date = t_hor,
    t_now = t_now,
    h = time_to_maturity,
    greeks
  )
}

