#' Specification of a solar scenario
#'
#' @param sm `spatialModel` object
#' @param sc `spatialCorrelation` object
#' @param places target places
#' @param from character, start Date for simulations in the format `YYYY-MM-DD`.
#' @param to character, end Date for simulations in the format `YYYY-MM-DD`.
#' @param exclude_known when true the two starting points (equals for all the simulations) will be excluded from the output.
#' @param quiet logical
#' @rdname spatialScenario_spec
#' @name spatialScenario_spec
#' @export
spatialScenario_spec <- function(sm, sc, places, from = "2010-01-01", to = "2010-01-31", exclude_known = FALSE, quiet = FALSE){

  if (missing(places)) {
    places <- sc$places
  } else {
    places <- sc$.__enclos_env__$private$check_places(places)
  }

  simSpec <- purrr::map(places, ~solarScenario_spec(sm$gridModel(.x), from = from, to = to, theta = 0,
                                                    exclude_known = exclude_known, quiet = TRUE))

  locations <- purrr::map_df(places, ~dplyr::bind_rows(sm$gridModel(.x)$coords))
  locations <- dplyr::bind_cols(place = places, locations)

  names(simSpec) <- places
  structure(
    list(
      spec = simSpec,
      locations = locations,
      from = as.Date(from),
      to = as.Date(to),
      sc = sc,
      places = places,
      nsim = NA,
      residuals = FALSE,
      filter = FALSE,
      quiet = quiet
    ),
    class = c("spatialScenarioSpec", "list")
  )
}

#' Simulate residuals from a a `spatialScenario_spec`
#'
#' @param simSpec object with the class `spatialScenario_spec`. See the function \code{\link{spatialScenario_spec}} for details.
#' @param nsim integer, number of simulations.
#' @param seed scalar integer, starting random seed.
#' @rdname spatialScenario_residuals
#' @name spatialScenario_residuals
#' @export
spatialScenario_residuals <- function(simSpec, nsim = 1, seed = 1){

  # Correlations matrices
  sc <- simSpec$sc
  # Places
  places <- simSpec$places

  # Random seed
  set.seed(seed)
  # Create a sequence of dates
  dates <- seq.Date(simSpec$from, simSpec$to, by = "1 day")
  # Initialize a dataset for storing simulations
  df_sim <- dplyr::tibble(date = dates) %>%
    dplyr::mutate(Year = lubridate::year(date), Month = lubridate::month(date)) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    tidyr::nest(data_B = date, data_X = date)

  # Simulating scenarios
  for(m in 1:nrow(df_sim)) {
    # Number of days of the month
    n <- df_sim$n[m]
    # Extract correlation matrices
    models_cr <- sc$get(places, nmonth = df_sim$Month[m])
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
    tidyr::nest()%>%
    dplyr::rename(X = "data")

  df_sim_B <- dplyr::bind_rows(df_sim$data_B) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest() %>%
    dplyr::rename(B = "data")

  df_sim <- dplyr::left_join(df_sim_X, df_sim_B, by = "nsim")
  df_sim <- dplyr::bind_cols(seed = seed, df_sim)

  for(i in 1:length(simSpec$spec)){
    place <- simSpec$spec[[i]]$place
    simSpec$spec[[i]]$residuals <- dplyr::mutate(df_sim,
                                                 X = purrr::map(X, ~dplyr::select(.x, date, any_of(place))),
                                                 B = purrr::map(B, ~dplyr::select(.x, date, any_of(place))))
  }

  simSpec$residuals <- TRUE
  simSpec$nsim <- nsim

  return(simSpec)
}

#' Simulate trajectories from a `spatialScenario_spec`
#'
#' @inheritParams spatialScenario_residuals
#'
#' @rdname spatialScenario_filter
#' @name spatialScenario_filter
#' @export
spatialScenario_filter <- function(simSpec){

  custom_msg <- function(place, n_models, i) {
    msg_1 <- paste0("Filtering place (", place, ")")
    msg_2 <- paste0(" ", i, "/", n_models, " ")
    msg_3 <- paste0(" (", format(i/n_models*100, digits = 3), "%)")
    paste0(msg_1, msg_2, msg_3)
  }

  n_models <- length(simSpec$spec)
  for(i in 1:n_models){
    place <- simSpec$spec[[i]]$place
    if (!simSpec$quiet) message(custom_msg(place, n_models, i), "\r", appendLF = FALSE)
    simSpec$spec[[i]] <- solarScenario_filter(simSpec$spec[[i]])
  }
  simSpec$filter <- TRUE
  return(simSpec)
}


