#' Compute monthly moments for mixture with 16 components
#' @param nb_models a list with 4 neighborhoods models
#' @param nmonths numeric vector, reference months. Default is `1:12`, i.e. all the months.
#' @param weights numeric, vector of 4 weights.
#' @noRd
#' @keywords spatialModel internal
#' @export
spatialModel_combinations = function(nb_models, nmonths = 1:12, weights = rep(1/4, 4), nobs.min = 5){
  # Weights in vector notation ensuring normalization
  w <- matrix(weights/sum(weights), ncol = 1)
  # Possible combinations in a grid with 4 reference locations
  possible_combinations <- data.frame(
    State = 1:16,
    Month = rep(1, 16),
    Location_A = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    Location_B = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0),
    Location_C = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
    Location_D = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  )
  # Add observations weights
  data_A <- dplyr::bind_cols(nb_models[[1]]$data)
  data_B <- dplyr::bind_cols(nb_models[[2]]$data)
  data_C <- dplyr::bind_cols(nb_models[[3]]$data)
  data_D <- dplyr::bind_cols(nb_models[[4]]$data)
  # Filter only for train data with non-zero weight
  df_A <- dplyr::filter(data_A, isTrain & weights != 0)
  df_B <- dplyr::filter(data_B, isTrain & weights != 0)
  df_C <- dplyr::filter(data_C, isTrain & weights != 0)
  df_D <- dplyr::filter(data_D, isTrain & weights != 0)
  # Select only relevant variables
  df_A <- dplyr::select(df_A, date, B_A = "B", u_A = "u_tilde")
  df_B <- dplyr::select(df_B, date, B_B = "B", u_B = "u_tilde")
  df_C <- dplyr::select(df_C, date, B_C = "B", u_C = "u_tilde")
  df_D <- dplyr::select(df_D, date, B_D = "B", u_D = "u_tilde")
  # Create a unique dataset
  data <- dplyr::left_join(df_A, df_B, by = "date") %>%
    dplyr::left_join(df_C, by = "date") %>%
    dplyr::left_join(df_D, by = "date") %>%
    dplyr::mutate(Month = lubridate::month(date))

  # Main loop
  m <- 1
  combinations <- list()
  for(m in nmonths){
    # Initialize the dataset
    combinations[[m]] <- dplyr::as_tibble(possible_combinations)
    combinations[[m]]$probs <- 0
    combinations[[m]]$mean <- 0
    combinations[[m]]$sd <- 0
    combinations[[m]]$n <- 0
    combinations[[m]]$cv <- as.list(1:16)
    # Filter dataset for reference month
    df <- dplyr::filter(data, Month == m)
    s <- 1
    for(s in 1:16){
      # Extract the states
      states <- unlist(possible_combinations[s,][-c(1:2)])
      # Compute joint probabilities
      combinations[[m]]$probs[s] <- mean(df$B_A == states[1] & df$B_B == states[2] & df$B_C == states[3] & df$B_D == states[4], na.rm = TRUE)
      # Extract the realized joint series
      df_u <- dplyr::filter(df, B_A == states[1] & B_B == states[2] & B_C == states[3] & B_D == states[4])
      df_u <- dplyr::select(df_u, dplyr::contains("u"))
      # Number of observations for the estimates
      combinations[[m]]$n[s] <- nrow(df_u)
      # If the number of cases is too low set probability to zero
      if (nrow(df_u) <= nobs.min) {
        combinations[[m]]$probs[s] <- 0
        combinations[[m]]$sd[s] <- 1
        next
      }
      # Variance-covariance matrix
      combinations[[m]]$cv[[s]] <- cov(df_u)
      # Expected value for the i-th state
      combinations[[m]]$mean[s] <- sum(dplyr::summarise_all(df_u, mean)*w)
      combinations[[m]]$sd[s] <- sqrt(t(w) %*% combinations[[m]]$cv[[s]] %*% w)[1]
    }
    # Normalize the probabilities
    combinations[[m]]$probs <- combinations[[m]]$probs/sum(combinations[[m]]$probs)
    combinations[[m]]$Month <- m
    colnames(combinations[[m]]) <- c("State", "Month", names(nb_models), "probs", "mean", "sd","n", "cv", "mu2")
  }
  dplyr::bind_rows(combinations)
}

#' Perform the bilinear interpolation for a target variable.
#' @param nb_models a list with 4 neighborhoods models
#' @param target character, name of the target variable to interpolate.
#' @param n number of neighborhoods to use for interpolation.
#' @param weights numeric, vector of weights of with length equal to n.
#' @param day_date date for interpolation, if missing all the available dates will be used.
#' @noRd
#' @keywords spatialModel internal
#' @export
spatialModel_interpolator = function(nb_models, target = "GHI", n = 4, weights = rep(1/n, n), day_date){
  # Initialize the output dataset
  interp_data <- dplyr::tibble(x3 = lubridate::as_date(NA), x4 = NA, x5 = FALSE)
  colnames(interp_data) <- c("date", target, "interpolated")
  # Neighborhoods models
  interp_data <- as.list(interp_data)
  if (!missing(day_date)) {
    data <- purrr::map(nb_models, ~dplyr::filter(.x$data, date %in% as.Date(day_date)))
    data <- purrr::map(data, ~.x[, c("date", target)])
  } else {
    data <- purrr::map(nb_models, ~.x$data[, c("date", target)])
  }
  df <- data[[1]]
  for(i in 2:n){
    df <- dplyr::left_join(df, data[[i]], by = "date")
  }
  interp_data$date <- df$date
  interp_data[[target]] <- 0
  for(i in 1:n){
    interp_data[[target]] <- interp_data[[target]] + df[,i+1][[1]]*weights[i]
  }
  interp_data <- dplyr::bind_cols(interp_data)
  interp_data$interpolated <- TRUE
  return(interp_data)
}

#' @noRd
#' @keywords spatialModel internal
#' @export
spatialModel_simulate_filter <- function(model, from = "2010-01-01", to = "2010-12-31", u, B, exclude_known = FALSE, quiet = FALSE){
  # AR parameters
  AR_model_Yt <- model$AR_model_Yt
  arOrder <- model$control$mean.model$arOrder
  # GARCH(p,q) model
  GARCH <- model$GARCH
  # Intercept
  GARCH$omega <- GARCH$coef[names(GARCH$coef) == "omega"]
  # Arch parameters
  GARCH$alpha <- GARCH$coef[stringr::str_detect(names(GARCH$coef), "alpha")]
  archOrder <- max(c(length(GARCH$alpha), 1))
  # Garch parameters
  GARCH$beta <- GARCH$coef[stringr::str_detect(names(GARCH$coef), "beta")]
  garchOrder <- max(c(length(GARCH$beta), 1))
  # Garch next step function
  GARCH_next_step <- GARCH_pq_next_step(GARCH$omega, GARCH$alpha, GARCH$beta)
  # Initialize empirical and simulated data
  sim_data <- solarModel_simulate_data(model, from, to, quiet)
  # Number of lags to consider
  i_start <- sim_data$i_start

  # Filter df_emp to be in [from - to] dates
  if (exclude_known) {
    sim_data$emp <- dplyr::filter(sim_data$emp, date >= from & date <= to)
  }

  j <- 1
  # Initialize dataset for storing the simulation
  df_sim <- sim_data$sim
  # Bernoulli jump
  df_sim$B <- B
  # Simulate Normal mixture (ut)
  df_sim$z <- u
  i <- i_start
  for(i in i_start:nrow(df_sim)){
    # Simulated GARCH standard deviation
    df_sim$sigma[i] <- GARCH_next_step(df_sim$eps_tilde[(i-archOrder):(i-1)], df_sim$sigma[(i-garchOrder):(i-1)])
    # Simulated Yt_tilde
    df_sim$Yt_tilde_hat[i] <- predict(AR_model_Yt, newdata = df_sim[(i-i_start):i,])[i_start]
    # Simulated normal mixture
    df_sim$u_tilde[i] <- (df_sim$mu_up[i] + df_sim$sd_up[i]*df_sim$z[i])*df_sim$B[i] + (1-df_sim$B[i])*(df_sim$mu_dw[i] + df_sim$sd_dw[i]*df_sim$z[i])
    # Simulated standardized monthly residuals
    df_sim$u[i] <- df_sim$sigma_m[i]*df_sim$u_tilde[i]
    # Simulated standardized residuals
    df_sim$eps_tilde[i] <- df_sim$sigma[i]*df_sim$u[i]
    # Simulated AR residuals
    df_sim$eps[i] <- df_sim$eps_tilde[i]*df_sim$sigma_bar[i]
    # Simulated Yt_tilde
    df_sim$Yt_tilde[i] <- df_sim$Yt_tilde_hat[i] + df_sim$eps[i]
    # Simulated Yt
    df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i] + df_sim$Yt_tilde_uncond[i]
    # Simulated Xt
    df_sim$Xt[i] <- model$transform$iY(df_sim$Yt[i])
    # Simulated GHI
    df_sim[[model$target]][i] <- model$transform$GHI(df_sim$Xt[i], df_sim$Ct[i])
  }
  # Remove redundant variables
  df_sim <- dplyr::select(df_sim, -mu_up, -mu_dw, -sd_up, -sd_dw, -p_up, -Ct, -Yt_bar, -sigma_bar, -sigma_m, -Yt_tilde_hat, -Yt_tilde_uncond, -z)
  # Remove initial values
  if (exclude_known) {
    df_sim <- dplyr::filter(df_sim, date >= from & date <= to)
  }

  structure(
    list(
      sim = df_sim,
      emp = sim_data$emp
    ),
    class = c("solarModelSimulation", "list")
  )
}

#' @noRd
#' @keywords spatialModel internal
#' @export
spatialModel_simulate_residuals <- function(sc, places, from = "2022-01-01", to = "2022-12-31", nsim = 10, seed = 1){

  seq_date <- seq.Date(as.Date(from)-2, as.Date(to), by = "1 day")

  df_sim <- dplyr::tibble(date = seq_date) %>%
    dplyr::mutate(Year = lubridate::year(date), Month = lubridate::month(date)) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    tidyr::nest(data_B = date, data_X = date)

  set.seed(seed)
  m <- 1
  for(m in 1:nrow(df_sim)) {
    # Extract correlation matrices
    models_cr <- sc$get(places, nmonth = df_sim$Month[m])
    # Number of simulations
    n <- df_sim$n[m]
    # Correlation matrix for mixture components
    cr_X <- models_cr$cr_X
    # 1) Complete Normal simulation
    sim_x12 <- mvtnorm::rmvnorm(n*nsim, mean = rep(0, ncol(cr_X)), sigma = cr_X)
    colnames(sim_x12) <- colnames(cr_X)
    sim_x12 <- dplyr::as_tibble(sim_x12)
    # 2) Bernoulli simulation
    sim_B <- bindata::rmvbin(n*nsim, margprob = models_cr$margprob, sigma = models_cr$sigma_B)
    colnames(sim_B) <- colnames(models_cr$sigma_B)
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
    colnames(sim_B) <- colnames(models_cr$sigma_B)
    colnames(sim_X) <- colnames(models_cr$sigma_B)

    sim_B <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_B)
    sim_B <- dplyr::group_by(sim_B, scenario) %>% tidyr::nest()
    sim_B <- dplyr::mutate(sim_B, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))

    sim_X <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_X)
    sim_X <- dplyr::group_by(sim_X, scenario) %>% tidyr::nest()
    sim_X <- dplyr::mutate(sim_X, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))

    df_sim$data_B[[m]] <- dplyr::bind_cols(df_sim$data_B[[m]], sim_B)
    df_sim$data_X[[m]] <- dplyr::bind_cols(df_sim$data_X[[m]], sim_X)
  }

  df_sim_X <- dplyr::bind_rows(df_sim$data_X) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest()

  df_sim_B <- dplyr::bind_rows(df_sim$data_B) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest()

  structure(
    list(
      x = dplyr::bind_cols(seed = seed, df_sim_X),
      B = dplyr::bind_cols(seed = seed, df_sim_B)
    )
  )
}
