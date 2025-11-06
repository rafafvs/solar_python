#' Control parameters for solar hedging
#'
#' @param n_panels Numeric scalar, number of meters squared of solar panels.
#' @param efficiency Numeric scalar, average electricity produced with 1 \eqn{m^2} of solar panels given 1 kWh/\eqn{m^2} of solar radiation received.
#' @param PUN Numeric scalar, fixed electricity price at which the produced energy is sold.
#' @param tick Numeric scalar, tick for the monetary conversion of the payoff of a solar derivative from kWh/\eqn{m^2} to Euros.
#' @param n_contract Numeric scalar, number of solar derivative contracts bought.
#'
#' @rdname control_solarHedging
#' @name control_solarHedging
#' @keywords control
#' @note Version 1.0.0.
#' @export
control_solarHedging <- function(n_panels = 1, efficiency = 1, PUN = 1, tick = 1, n_contracts = 1){
  structure(
    list(
      n_panels = n_panels,
      efficiency = efficiency,
      PUN = PUN,
      tick = tick,
      n_contracts = n_contracts
    ),
    class = c("control", "list")
  )
}

#' Global Minimum Variance number of contracts
#'
#' Compute the optimal number of contracts, such that the variance of the cash flow of a solar power producer is minimum.
#'
#' @param model An object with the class `solarModel`. See the function \code{\link{solarModel}} for more details.
#' @param moments Tibble containing the forecasted moments for different days ahead. See the function \code{\link{solarMoments}} for more details.
#' @param P0_P Optional numeric scalar, expected value of 1 solar derivative with unitary tick under \eqn{\mathbb{P}}.
#' @param gamma Numeric scalar, risk aversion parameter.
#' @inheritParams solarOption_historical
#' @param control_hedge Named list, control parameters for hedging. See the function \code{\link{control_solarHedging}} for more details.
#'
#' @rdname solarHedging_model
#' @name solarHedging_model
#' @keywords solarHedging
#' @note Version 1.0.0.
#' @export
solarHedging_model <- function(model, moments, P0_P, r_star, gamma = 0.01, put = TRUE, control_options = control_solarOption(), control_hedge = control_solarHedging()){
  # Select only time variable
  data <- dplyr::select(moments, date, Year, Month, Day)
  # Initialize a column for the number of contracts
  data$w <- 0
  # Maximum and minimum of clearness index
  K_min_max <- model$transform$bounds("Kt")
  i <- 1
  .l <- list()
  for(i in 1:nrow(moments)){
    # Extract moments
    mom <- moments[i,]
    # Density of Yt
    pdf_Yt <- function(x) dmixnorm(x, mean = c(mom$M_Y1, mom$M_Y0), sd = c(mom$S_Y1, mom$S_Y0), alpha = c(mom$p1, 1-mom$p1))
    # Maximum and minimum of solar radiation
    RT_min_max <- mom$Ct * K_min_max
    # Density of Rt
    pdf_Rt <- function(x) dsolarGHI(x, mom$Ct, mom$alpha, mom$beta, pdf_Yt)
    # Payoff function
    if (put) {
      # Put function
      Gamma <- function(x, K, Ct) (K - x)*ifelse(x > K, 0, 1)
    } else {
      # Call function
      Gamma <- function(x, K, Ct) (x - K)*ifelse(x < K, 0, 1)
    }
    # Strike price
    K <- mom$GHI_bar * exp(control_options$K)
    # Expected value of Gamma
    e_Gamma <- integrate(function(x) Gamma(x, K, mom$Ct) * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Second moment of Gamma
    e2_Gamma <- integrate(function(x) Gamma(x, K, mom$Ct)^2 * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Variance payoff
    v_Gamma <- e2_Gamma - e_Gamma^2
    # Joint expected value
    e_RT_K_plus <- integrate(function(x) x * K * ifelse(x > K, 0, 1) * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Second moment of RT
    e_RT2_plus <- integrate(function(x) (x^2) * ifelse(x > K, 0, 1) * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Expected value of RT
    e_RT <- integrate(function(x) x * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Expected value of RT2
    e_RT2 <- integrate(function(x) x^2 * pdf_Rt(x), lower = RT_min_max[1], upper = RT_min_max[2])$value
    # Variance of RT
    v_RT <- e_RT2 - e_RT
    # Covariance between RT and Gamma
    cv_RT_Gamma <- e_RT_K_plus - e_RT2_plus - e_RT * e_Gamma

    # Exposure to solar radiation
    w_R <- control_hedge$n_panels * control_hedge$PUN * control_hedge$efficiency

    # Expectation unhedged
    e_Pi_uh <- w_R * e_RT
    # Variance unhedged
    v_Pi_uh <- w_R^2 * v_RT

    P0_Q <- P0_P * (1 + r_star)
    # Minimum variance weights
    w_gmv <- - w_R / control_hedge$tick * (cv_RT_Gamma / v_Gamma)
    # Optimal mean-variance weights
    w_mv <- w_gmv + control_hedge$tick * (P0_P - P0_Q) / (gamma * control_hedge$tick^2 * v_Gamma)

    # Variance hedged
    v_Pi_mv <- v_Pi_uh +  w_mv^2 * control_hedge$tick^2 * v_Gamma + 2 * w_R * w_mv * control_hedge$tick * cv_RT_Gamma
    # Expectation hedged
    e_Pi_mv <- e_Pi_uh + w_mv * control_hedge$tick * (P0_P - P0_Q)

    # Variance hedged
    v_Pi_gmv <- v_Pi_uh +  w_gmv^2 * control_hedge$tick^2 * v_Gamma + 2 * w_R * w_gmv * control_hedge$tick * cv_RT_Gamma
    # Expectation hedged
    e_Pi_gmv <- e_Pi_uh + w_gmv * control_hedge$tick * (P0_P - P0_Q)

    .l[[i]] <- dplyr::bind_cols(data[,i],
                                dplyr::tibble(
                                  e_Pi_uh = e_Pi_uh,
                                  e_Pi_mv = e_Pi_mv,
                                  e_Pi_gmv = e_Pi_gmv,
                                  sd_Pi_uh = sqrt(v_Pi_uh),
                                  sd_Pi_mv = sqrt(v_Pi_mv),
                                  sd_Pi_gmv = sqrt(v_Pi_gmv),
                                  satis_uh = e_Pi_uh - gamma/2 * sd_Pi_uh^2,
                                  satis_mv = e_Pi_mv - gamma/2 * sd_Pi_mv^2,
                                  satis_gmv = e_Pi_gmv - gamma/2 * sd_Pi_gmv^2,
                                  ratio = (satis_mv/satis_uh - 1),
                                  w_gmv = w_gmv,
                                  w_mv = w_mv,
                                  w_R = w_R,
                                  tick = control_hedge$tick,
                                  r_max = gamma * tick * w_gmv * v_Gamma / P0_P
                                ))
  }
  data <- dplyr::bind_rows(.l)
  return(data)
}

#' Compute the optimal number of solar derivative index
#'
#' Compute the optimal number of contracts, such that the variance of the cash flow of a solar power
#' producer with a given setup is minimum.
#'
#' @param scenarios Object with the class `solarScenario`. See the function \code{\link{solarScenario}} for more details.
#' @inheritParams solarHedging_model
#'
#' @rdname solarHedging_scenarios
#' @name solarHedging_scenarios
#' @keywords solarHedging
#' @note Version 1.0.0.
#' @export
solarHedging_scenarios <- function(scenarios, P0_P, r_star = 0, gamma = 0.0001, put = TRUE, control_options = control_solarOption(), control_hedge = control_solarHedging()){
  # 0) Compute the payoff of solar options
  payoff_sim <- scenarios$sim %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::mutate(n = number_of_day(date)) %>%
    dplyr::mutate(strike = GHI_bar * exp(control_options$K))
  if (put) {
    payoff_sim <- payoff_sim %>%
      dplyr::mutate(
        exercise = ifelse(strike > GHI, 1, 0),
        payoff = (strike - GHI) * exercise)
  } else {
    payoff_sim <- payoff_sim %>%
      dplyr::mutate(
        exercise = ifelse(strike < GHI, 1, 0),
        payoff = (GHI-strike) * exercise)
  }
  # ****************************************
  # 1) Compute the moments of the yearly sum
  payoff_sim <- payoff_sim %>%
    dplyr::group_by(Year, Month, Day) %>%
    dplyr::mutate(premium = mean(payoff)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Year, nsim) %>%
    dplyr::summarise(Gamma = sum(payoff),
                     strike = sum(strike),
                     Rt = sum(GHI),
                     V0 = sum(premium))
  # ****************************************
  # 2.A) Expectations
  ## Expectation of RT
  e_RT <- mean(payoff_sim$Rt)
  ## Expectation of Gamma
  e_Gamma <- mean(payoff_sim$Gamma)
  # 2.B) Variances
  ## Variance of RT
  v_RT <- var(payoff_sim$Rt)
  ## Variance of Gamma
  v_Gamma <- var(payoff_sim$Gamma)
  # 2.C) Joint expectations RT and Gamma
  e_RT_Gamma <- mean(payoff_sim$Rt * payoff_sim$Gamma)
  # 2.D) Covariances RT and Gamma
  cv_RT_Gamma <- e_RT_Gamma - e_RT * e_Gamma
  # 2.E) Correlation RT and Gamma
  cr_RT_Gamma <- cv_RT_Gamma / sqrt(v_RT * v_Gamma)
  # ****************************************
  # Total exposure to solar radiation
  w_R <- control_hedge$n_panels * control_hedge$PUN * control_hedge$efficiency
  # 3.A) Expectation and variance (unhedged)
  e_Pi_uh <- w_R * e_RT
  v_Pi_uh <- w_R^2 * v_RT
  # 3.B) Simulated cash flows (unhedged)
  Pi_uh <- w_R * payoff_sim$Rt
  # ****************************************
  # 4) Compute price with premium charged by the seller
  P0_P <- ifelse(missing(P0_P), mean(payoff_sim$Gamma), P0_P)
  P0_Q <- P0_P * exp(r_star)
  # ****************************************
  # 4.A) Global minimum variance number of contracs
  w_gmv <- - w_R / control_hedge$tick * (cv_RT_Gamma / v_Gamma)
  # Maximum premium tolerated
  r_max <- gamma * control_hedge$tick * w_gmv * v_Gamma / P0_P
  # 4.B) Expectation and variance (hedged GMV)
  e_Pi_gmv <- e_Pi_uh + w_gmv * control_hedge$tick * (P0_P - P0_Q)
  v_Pi_gmv <- v_Pi_uh +  w_gmv^2 * control_hedge$tick^2 * v_Gamma + 2 * w_R * w_gmv * control_hedge$tick * cv_RT_Gamma
  # 4.C) Simulated cash flows (hedged GMV)
  Pi_gmv <- Pi_uh + w_gmv * control_hedge$tick * (payoff_sim$Gamma - P0_Q)
  # ****************************************
  # 5.A) Optimal mean-variance number of contracs
  w_mv <- w_gmv + control_hedge$tick * (P0_P - P0_Q) / (gamma * control_hedge$tick^2 * v_Gamma)
  # 5.B) Expectation and variance (hedged MV)
  e_Pi_mv <- e_Pi_uh + w_mv * control_hedge$tick * (P0_P - P0_Q)
  v_Pi_mv <- v_Pi_uh +  w_mv^2 * control_hedge$tick^2 * v_Gamma + 2 * w_R * w_mv * control_hedge$tick * cv_RT_Gamma
  # 5.C) Simulated cash flows (hedged MV)
  Pi_mv <- Pi_uh + w_mv * control_hedge$tick * (payoff_sim$Gamma - P0_Q)
  # ****************************************

  dplyr::tibble(
    # Option tick
    tick = control_hedge$tick,
    # Maximum price
    P0_P = P0_P,
    # Maximum price
    P0_Q = P0_P * (1 + r_max),
    # Premium charged
    r_star = r_star,
    # Risk adversion
    gamma = gamma,
    # Maximum premium tolerated
    r_max = r_max,
    # Optimal number of contracts
    w_gmv = w_gmv,
    w_mv = w_mv,
    w_R = w_R,
    # Expectations
    e_Pi_uh = e_Pi_uh,
    e_Pi_mv = e_Pi_mv,
    e_Pi_gmv = e_Pi_gmv,
    # Std. deviations
    sd_Pi_uh = sqrt(v_Pi_uh),
    sd_Pi_mv = sqrt(v_Pi_mv),
    sd_Pi_gmv = sqrt(v_Pi_gmv),
    # Mean-variance satisfactions
    satis_uh = e_Pi_uh - gamma/2 * sd_Pi_uh^2,
    satis_mv = e_Pi_mv - gamma/2 * sd_Pi_mv^2,
    satis_gmv = e_Pi_gmv - gamma/2 * sd_Pi_gmv^2,
    # Ratio MV vs Unhedged (should be > 1)
    ratio = (satis_mv/satis_uh - 1),
    # Simulated cash flows
    Pi_uh = list(Pi_uh),
    Pi_mv = list(Pi_mv),
    Pi_gmv = list(Pi_gmv),
    # Correlations
    cr_RT_Gamma = cr_RT_Gamma,
    # Controls parameters
    control_options = list( control_options),
    control_hedge = list( control_hedge)
  )
}



