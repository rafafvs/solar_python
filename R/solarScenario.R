#' Simulate multiple scenarios
#'
#' Simulate multiple scenarios of solar radiation with a `solarModel` object.
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param from character, start Date for simulations in the format `YYYY-MM-DD`.
#' @param to character, end Date for simulations in the format `YYYY-MM-DD`.
#' @param by character, steps for multiple scenarios, e.g. `1 day` (day-ahead simulations), `15 days`, `1 month`, `3 months`, ecc.
#' For each step are simulated `nsim` scenarios.
#' @param nsim integer, number of simulations.
#' @param theta numeric, shift parameter for the mixture.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' scen <- solarScenario(model, "2005-01-10", "2020-01-01", theta = 0, nsim = 4, by = "1 year")
#' # Plot
#' solarScenario_plot(scen, nsim = 3)
#' # Solar Option
#' solarOption_scenario(model, scen)
#' solarOption_historical(model)
#'
#' @rdname solarScenario
#' @name solarScenario
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario <- function(model, from = "2010-01-01", to = "2011-01-01", by = "1 year", theta = 0, nsim = 1, seed = 1, quiet = FALSE){

  idx_date <- seq.Date(as.Date(from), as.Date(to), by = by)
  scenarios <- list()
  df_emp <- dplyr::tibble()
  n_scenario <- length(idx_date)
  j <- 2
  for(j in 2:n_scenario){
    if (!quiet) {
      # To report progress
      pb <- txtProgressBar(min = 1,            # Minimum value of the progress bar
                           max = n_scenario,   # Maximum value of the progress bar
                           style = 3,          # Progress bar style (also available style = 1 and style = 2)
                           width = 50,         # Progress bar width. Defaults to getOption("width")
                           char = "#")
      setTxtProgressBar(pb, j)
    }
    simSpec <- solarScenario_spec(model, from = idx_date[j-1], to = idx_date[j]-1, theta = theta, exclude_known = TRUE, quiet = TRUE)
    simSpec <- solarScenario_residuals(simSpec, nsim = nsim, seed = seed)
    simSpec <- solarScenario_filter(simSpec)
    df_emp <- dplyr::bind_rows(df_emp, simSpec$emp)
    scenarios[[j]] <- dplyr::bind_cols(seed = seed, dplyr::bind_rows(simSpec$simulations))
    seed <- seed + j - 1
  }
  if (!quiet) close(pb)

  simSpec$emp <- df_emp[!duplicated(df_emp),]
  simSpec$simulations <- scenarios
  return(as_solarScenario(simSpec, by))
}

#' Specification of a solar scenario
#'
#' @inheritParams solarScenario
#' @param exclude_known when true the two starting points (equals for all the simulations) will be excluded from the output.
#'
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' simSpec <- solarScenario_spec(model)
#'
#' @rdname solarScenario_spec
#' @name solarScenario_spec
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario_spec <- function(model, from = "2010-01-01", to = "2010-12-31", theta = 0, exclude_known = FALSE, quiet = FALSE){
  # Extract information
  data <- model$data
  place <- model$place
  # Maximum number of lags to start
  i_start <- max(c(model$ARMA$order, model$GARCH_model$order)) + 1
  # Initial date
  from <- as.Date(from)
  # End date
  to <- as.Date(to)
  # Selected columns
  cols_emp <- c("date", "n", "Year", "Month", "Day", "GHI", "clearsky",
                 "Xt", "Yt", "Yt_tilde", "Yt_tilde_hat", "eps", "eps_tilde", "sigma", "u_tilde", "B", "z1", "z2")
  cols_sim <- c(cols_emp, "Ct", "Yt_bar", "GHI_bar", "sigma_bar", "Yt_tilde_uncond", "sigma_uncond", "mu1", "mu2", "sd1", "sd2", "p1")

  # Initialize a dataset
  max_date_from <- max(data$date)
  max_date_to <- max_date_from - i_start
  if (max_date_to >= to) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)) & date <= to)
    df_emp <- dplyr::bind_cols(place = place, df_emp)
  } else if (max_date_to >= from & max_date_from >= from) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)))
    df_new_emp <- dplyr::tibble(date = seq.Date(max(df_emp$date) + 1, to, by = "1 day"))
    df_emp <- dplyr::bind_rows(df_emp, df_new_emp)
    df_emp <- dplyr::mutate(df_emp,
                            Year = lubridate::year(date),
                            Month = lubridate::month(date),
                            Day = lubridate::day(date))
    df_emp$n <- solarr::number_of_day(df_emp$date)
    # Add seasonal variables
    df_emp <- df_emp[, cols_emp]
    df_emp <- dplyr::left_join(df_emp, model$seasonal_data, by = c("Month", "Day", "n"))
  } else {
    msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
    if (!quiet) warning(msg)
    return(NULL)
  }

  # Initialize simulation dataset
  df_sim <- df_emp[, cols_sim]
  # Initialize lambda
  df_sim$theta <- theta
  # Filter df_emp to be in [from - to] dates
  df_emp <- df_emp[, cols_emp]
  if (exclude_known) {
    df_emp <- dplyr::filter(df_emp, date >= from & date <= to)
  }

  # Output structure
  structure(
    list(
      # Base dataset for simulations
      sim = df_sim,
      # Dataset with empirical data
      emp = df_emp,
      # Model place
      place = model$place,
      # Coordinates
      coords = model$spec$coords,
      # Model target
      target = model$spec$target,
      # Model transform
      transform = model$transform,
      # Monthly probability
      p_up = function(nmonth){model$NM_model$coefficients$p1[nmonth]},
      # AR model
      ARMA = model$ARMA,
      # GARCH model
      GARCH_model = model$GARCH,
      # Number of lags for simulations
      i_start = i_start,
      # Other info
      exclude_known = exclude_known,
      from = from,
      to = to,
      seed = NA,
      nsim = NA,
      quiet = quiet,
      # List to store residuals and simulations
      residuals = list(),
      simulations = list()
    ),
    class = c("solarScenarioSpec", "list")
  )
}

#' Simulate residuals for a `solarScenario_spec`
#'
#' @inheritParams solarScenario
#' @param simSpec object with the class `solarScenario_spec`. See the function \code{\link{solarScenario_spec}} for details.
#'
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' simSpec <- solarScenario_spec(model, from = "2010-01-01", to = "2010-01-01")
#' simSpec <- solarScenario_residuals(simSpec, nsim = 10)
#' @rdname solarScenario_residuals
#' @name solarScenario_residuals
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario_residuals <- function(simSpec, nsim = 1, seed = 1){
  # Random seed
  set.seed(seed)
  # Initialize the sequence of dates
  # Initialize a dataset for storing simulations
  df_sim <- dplyr::tibble(date = simSpec$sim$date[-c(1:(simSpec$i_start-1))]) %>%
    dplyr::mutate(Year = lubridate::year(date), Month = lubridate::month(date)) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    tidyr::nest(data_B = date, data_X = date)

  m <- 1
  # Simulating scenarios
  for(m in 1:nrow(df_sim)) {
    # Number of days of the month
    n <- df_sim$n[m]
    # Complete Normal simulation
    sim_x12 <- mvtnorm::rmvnorm(n*nsim, mean = rep(0, 2), sigma = diag(c(1,1)))
    colnames(sim_x12) <- c("x1_1", "x2_1")
    sim_x12 <- dplyr::as_tibble(sim_x12)
    # 2) Bernoulli simulation
    # sim_B <- rbinom(n*nsim, 1, prob = simSpec$p_up(m))
    sim_B <- bindata::rmvbin(n*nsim, margprob = simSpec$p_up(df_sim$Month[m]))
    colnames(sim_B) <- simSpec$place
    sim_B <- dplyr::as_tibble(sim_B)
    # Slit the mixture simulations
    sim_x1 <- dplyr::select(sim_x12, dplyr::contains("x1_"))
    sim_x2 <- dplyr::select(sim_x12, dplyr::contains("x2_"))
    sim_X <- sim_x1
    for(i in 1:nrow(sim_B)){
      for(j in 1:ncol(sim_B)){
        sim_X[i,j] <- ifelse(sim_B[i,j] == 1, sim_x1[i,j], sim_x2[i,j])
      }
    }
    colnames(sim_B) <- simSpec$place
    colnames(sim_X) <- simSpec$place
    # Structure binomial simulations
    sim_B <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_B)
    sim_B <- tidyr::nest(dplyr::group_by(sim_B, scenario))
    sim_B <- dplyr::mutate(sim_B, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))
    # Structure normal simulations
    sim_X <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_X)
    sim_X <- dplyr::group_by(sim_X, scenario) %>% tidyr::nest()
    sim_X <- dplyr::mutate(sim_X, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))
    # Store data
    df_sim$data_B[[m]] <- dplyr::bind_cols(df_sim$data_B[[m]], sim_B)
    df_sim$data_X[[m]] <- dplyr::bind_cols(df_sim$data_X[[m]], sim_X)
  }
  # Structure binomial simulations
  df_sim_X <- dplyr::bind_rows(df_sim$data_X) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest()%>%
    dplyr::rename(X = "data")
  # Structure normal simulations
  df_sim_B <- dplyr::bind_rows(df_sim$data_B) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest() %>%
    dplyr::rename(B = "data")
  # Create a unique dataset
  df_sim <- dplyr::left_join( df_sim_X, df_sim_B, by = "nsim")
  df_sim <- dplyr::bind_cols(seed = seed, df_sim)
  simSpec$residuals <- df_sim
  simSpec$nsim <- nsim
  simSpec$seed <- seed
  return(simSpec)
}

#' Simulate trajectories from a a `solarScenario_spec`
#'
#' @inheritParams solarScenario
#' @inheritParams solarScenario_residuals
#'
#' @examples
#' model <- Bologna
#' simSpec <- solarScenario_spec(model, from = "2023-01-01", to = "2023-12-31")
#' simSpec <- solarScenario_residuals(simSpec, nsim = 1, seed = 1)
#' simSpec <- solarScenario_filter(simSpec)
#' # Empiric data
#' df_emp <- simSpec$emp
#' # First simulation
#' df_sim <- simSpec$simulations[[1]]
#' ggplot()+
#' geom_line(data = df_emp, aes(date, GHI))+
#' geom_line(data = df_sim, aes(date, GHI), color = "red")
#'
#' @rdname solarScenario_filter
#' @name solarScenario_filter
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario_filter <- function(simSpec){

  if (purrr::is_empty(simSpec$residuals)) {
    stop("The slot `simSpec$residuals` is empty! Consider running `simSpec <- solarScenario_residuals(simSpec)` before!")
  }
  # Number of lags to start
  i_start <- simSpec$i_start
  # ARMA order
  arma_order <- simSpec$ARMA$order
  # GARCH order
  garch_order <- simSpec$GARCH_model$order
  # Number of simulations
  nsim <- nrow(simSpec$residuals)
  j <- 1
  simulations <- list()
  for(j in 1:nsim){
    # Initialize dataset for storing the simulation
    df_sim <- simSpec$sim
    # Add simulated residuals
    df_sim$z <- 0
    df_sim$z[(i_start):nrow(df_sim)] <- simSpec$residuals[j,]$X[[1]][, simSpec$place][[1]]
    df_sim$B[(i_start):nrow(df_sim)] <- simSpec$residuals[j,]$B[[1]][, simSpec$place][[1]]
    # Verbose message
    if (!simSpec$quiet) message("Simulation: ", j, "/", nsim, " (", round(j/nsim*100, 4), " %) \r", appendLF = FALSE)
    # Routine
    i <- i_start
    for(i in i_start:nrow(df_sim)){
      # ARCH component
      eps_tilde <- c(1)
      if (garch_order[1] > 0){
        eps_tilde <- df_sim$eps_tilde[(i-1):(i-garch_order[1])]
      }
      # GARCH component
      sigma_t <- c(1)
      if (garch_order[2] > 0){
        sigma_t <- df_sim$sigma[(i-1):(i-garch_order[2])]
      }
      # Simulated GARCH next step standard deviation (sigma)
      df_sim$sigma[i] <- sqrt(simSpec$GARCH_model$next_step(eps_tilde, sigma_t^2))

      # AR component
      Yt_tilde <- c()
      if (arma_order[1] > 0){
        Yt_tilde <- df_sim$Yt_tilde[(i-1):(i-arma_order[1])]
      }
      # MA component
      eps <- c()
      if (arma_order[2] > 0){
        eps <- df_sim$eps[(i-1):(i-arma_order[2])]
      }
      # Simulated ARMA next step (Yt_tilde)
      df_sim$Yt_tilde_hat[i] <- simSpec$ARMA$next_step(c(Yt_tilde, eps), 1)[1,1]
      df_sim$z[i] <- df_sim$z[i] + df_sim$theta[i]
      # Simulated monthly normal mixture (u_tilde)
      df_sim$u_tilde[i] <- (df_sim$mu1[i] + df_sim$sd1[i]*df_sim$z[i])*df_sim$B[i] + (df_sim$mu2[i] + df_sim$sd2[i]*df_sim$z[i])*(1-df_sim$B[i])
      # Simulated deseasonalized residuals (eps_tilde)
      df_sim$eps_tilde[i] <- df_sim$sigma[i] * df_sim$u_tilde[i]
      # Simulated ARMA residuals (eps)
      df_sim$eps[i] <- df_sim$eps_tilde[i] * df_sim$sigma_bar[i] * df_sim$sigma_uncond[i]
      # Simulated deseasonalized series (Yt_tilde)
      df_sim$Yt_tilde[i] <- df_sim$Yt_tilde_hat[i] + df_sim$eps[i]
      # Simulated transformed variable (Yt)
      df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i] + df_sim$Yt_tilde_uncond[i]
    }
    # Simulated cloudiness ratio (Xt)
    df_sim$Xt <- simSpec$transform$iX_prime(simSpec$transform$iY(df_sim$Yt))
    # Simulated solar radiation (GHI)
    df_sim[[simSpec$target]] <- simSpec$transform$iX(df_sim$Xt, df_sim$Ct)
    # Remove redundant variables
    df_sim <- dplyr::select(df_sim, -mu1, -mu2, -sd1, -sd2, -p1, -Ct, -Yt_bar, -sigma_bar, -Yt_tilde_hat, -Yt_tilde_uncond, -sigma_uncond)
    # Remove initial values
    if (simSpec$exclude_known) {
      df_sim <- dplyr::filter(df_sim, date >= simSpec$from & date <= simSpec$to)
    }
    # Store simulations
    simulations[[j]] <- dplyr::bind_cols(nsim = j, df_sim)
  }
  # Add simulations
  simSpec$simulations <- simulations
  return(simSpec)
}

#' Extract and structure simulations from a `solarScenarioSpec`
#'
#' @inheritParams solarScenario_residuals
#' @param by Optional character. Represent the steps used for multiple scenarios.
#' @rdname as_solarScenario
#' @name as_solarScenario
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
as_solarScenario <- function(simSpec, by) {

  if (purrr::is_empty(simSpec$simulations)) {
    stop("The slot `simSpec$simulations` is empty! Consider running `simSpec <- solarScenario_filter(simSpec)` before!")
  }
  df_sim <- dplyr::bind_rows(simSpec$simulations) %>%
    dplyr::group_by(date, Year, Month, Day) %>%
    tidyr::nest() %>%
    dplyr::ungroup()

  structure(
    list(
      emp = simSpec$emp,
      sim = df_sim,
      target = simSpec$target,
      by = ifelse(missing(by), NA_character_, by)
    ),
    class = c("solarScenario", "list")
  )
}

#' Compute the Value at Risk from simulated values
#'
#' @param scenarios An object of the class `solarScenario`
#' @param alpha Confidence level for the VaR
#'
#' @examples
#' model <- Bologna
#' scen <- solarScenario(model, "2016-01-01", "2017-01-01", nsim = 10, by = "1 month")
#' solarScenario_VaR(scen, 0.05)
#' @rdname solarScenario_VaR
#' @name solarScenario_VaR
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario_VaR <- function(scenarios, alpha = 0.05){
  # Compute the VaR
  simulated <- scenarios$sim %>%
    tidyr::unnest(cols = c("data")) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(VaR_alpha = quantile(GHI, probs = alpha)) %>%
    dplyr::select(date, VaR_alpha) %>%
    dplyr::ungroup() %>%
    unique()
  # Empiric values
  data <- dplyr::left_join(simulated, scenarios$emp, by = "date")
  # Violations of VaR
  data$et <- ifelse(data$GHI <= data$VaR_alpha, 1, 0)
  # Select only relevant variables
  data <- dplyr::select(data, date, Year, Month, Day, VaR_alpha, et)
  # Output structure
  structure(
    list(
      data = data,
      alpha = alpha
    )
  )
}


