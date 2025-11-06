#' Compute the Value at Risk and Expected Shortfall of a SolarModel
#'
#' @param model solarModel
#' @param alpha Numeric vector of confidence levels. Allows for more than one alpha.
#' @param ci Numeric scalar, confidence levels used to state if the Null is rejected or not on VaR tests.
#' @param ES Logical, when `TRUE` the expected shortfall will be also computed for each alpha.
#' @param type Numeric, type of evaluation, `full` on the complete data, `train` on the train data, `test` on the test data.
#'
#' @rdname solarModel_VaR
#' @name solarModel_VaR
#' @export
solarModel_VaR <- function(model, alpha = 0.05, ci = 0.05, ES = FALSE, type = "full"){
  # Link function
  link <- model$spec$transform$link
  # Conditional moments
  moments <- model$moments$conditional
  # Type of VaR
  type <- match.arg(type, choices = c("full", "train", "test"))
  # Filter moments
  if (type == "train") {
    moments <- moments[model$data$isTrain,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  # Realized solar radiation
  Rt <- dplyr::filter(model$data, date %in% moments$date)$GHI
  # Initialization of the dataset
  moments$VaR <- NA
  moments$GHI <- Rt
  if (ES) {
    moments$ES <- NA
  }
  # Base dataset
  data <- dplyr::select(moments, date, GHI)
  # Number of observations
  N <- nrow(moments)
  # Number of VaRs
  k <- length(alpha)

  # 1) Compute VaR, Violations and ES
  # Initialize the matrices to store the output
  VaR_alpha <- Viol_VaR_alpha <- ES_alpha <- matrix(0, nrow = N, ncol = k)
  # Assign standard names
  colnames(VaR_alpha) <- paste0("VaR_", alpha)
  colnames(Viol_VaR_alpha) <- paste0("e_", alpha)
  if (ES) {
    ES_alpha <- matrix(0, nrow = N, ncol = k)
    colnames(ES_alpha) <- paste0("ES_", alpha)
  }
  # Iterate over the sequence of dates
  for(t in 1:N){
    # Moments
    df_n <- moments[t,]
    # Computation of Expected Shortfall or not
    # Quantile of Yt
    cdf_Y <- function(x) pmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
    # Quantile function of Rt
    Q_R <- function(alpha) qsolarGHI(alpha, df_n$Ct, df_n$alpha, df_n$beta, cdf_Y, link = link)
    # Value at Risk with probability `ci`
    VaR_alpha[t,] <- Q_R(alpha)
    # Violation of the VaR
    Viol_VaR_alpha[t,] <- ifelse(Rt[t] <= VaR_alpha[t,], 1, 0)
    if (ES) {
      # Quantile of Yt
      pdf_Y <- function(x) dmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
      # Quantile function of Rt
      f_R <- function(alpha) dsolarGHI(alpha, df_n$Ct, df_n$alpha, df_n$beta, pdf_Y, link = link)
      # Expected shortfall
      ES_R <- function(VaR, alpha) integrate(function(x) x*f_R(x), lower = df_n$Ct*(1-df_n$alpha-df_n$beta),
                                             upper = VaR, subdivisions = 50L, stop.on.error = FALSE)$value/alpha
      # Expected Shortfall at probability `ci`
      ES_alpha[t,] <- purrr::map2_dbl(VaR_alpha[t,], alpha, ~ES_R(.x, .y))
    }
  }

  # 2) Perform the tests on the violations
  # Tests on the VaR
  VaR_tests <- purrr::map_df(1:length(alpha), ~VaR_test(Viol_VaR_alpha[,.x], alpha[.x], ci = ci))

  # 3) Summarise the results
  # Initialize a matrix to store the VaR
  VaR_alpha_emp <- dplyr::bind_rows(VaR_alpha[1,])
  VaR_alpha_emp[1,] <- NA
  if (ES) {
    # Initialize a matrix to store the model's Expected Shortfall
    ES_alpha_mod <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_mod[1,] <- NA
    # Initialize a matrix to store the empiric Expected Shortfall
    ES_alpha_emp <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_emp[1,] <- NA
  }
  # Evaluate Empiric VaR and Expected Shortfalls
  for(i in 1:k){
    # Index of the violations
    idx_violations <- which(Viol_VaR_alpha[,i] == 1)
    if (!purrr::is_empty(idx_violations)){
      VaR_alpha_emp[1,i] <- sum(Viol_VaR_alpha[,i] == 1)/nrow(Viol_VaR_alpha)
      if (ES) {
        ES_alpha_mod[1,i] <- mean(ES_alpha[idx_violations,i])
        ES_alpha_emp[1,i] <- mean(Rt[idx_violations])
      }
    }
  }

  # Standard output
  output <- structure(
    list(
      VaR = dplyr::bind_cols(data, VaR_alpha),
      Viol = dplyr::bind_cols(data, Viol_VaR_alpha),
      VaR_emp = VaR_alpha_emp,
      VaR_test = VaR_tests
    )
  )
  if (ES) {
    output$ES <- dplyr::bind_cols(data, ES_alpha)
    output$ES_tot <- dplyr::bind_cols(Type = c("Model", "Empiric"),
                                      dplyr::bind_rows(ES_alpha_mod, ES_alpha_emp))
  }
  return(output)
}

#' Compute the Value at Risk and Expected Shortfall of a SolarMixture
#'
#' @param model solarMixture
#' @param alpha Numeric vector of confidence levels. Allows for more than one alpha.
#' @param ci Numeric scalar, confidence levels used to state if the Null is rejected or not on VaR tests.
#' @param ES Logical, when `TRUE` the expected shortfall will be also computed for each alpha.
#' @param type Numeric, type of evaluation, `full` on the complete data, `train` on the train data, `test` on the test data.
#'
#' @rdname solarMixture_VaR
#' @name solarMixture_VaR
#' @export
solarMixture_VaR <- function(solarMix, date, x, alpha = 0.05, ci = 0.05, ES = FALSE){
  # Create a base dataset
  data <- tibble(date = date, x = x)
  # Number of observations
  N <- length(date)
  # Number of VaRs
  k <- length(alpha)
  # 1) Compute VaR, Violations and ES
  VaR_alpha <- solarMix$VaR(date, alpha)
  # Violations of the VaR
  Viol_VaR_alpha <- VaR_alpha
  for(i in 1:k){
    Viol_VaR_alpha[,i] <- ifelse(x < VaR_alpha[,1], 1, 0)
  }
  colnames(Viol_VaR_alpha) <- paste0("e_", alpha)
  # Expected shortfall
  if (ES) {
    ES_alpha <- sm$ES(date, alpha)
  }

  # 2) Perform the tests on the violations
  # Tests on the VaR
  VaR_tests <- purrr::map_df(1:length(alpha), ~VaR_test(Viol_VaR_alpha[,.x], alpha[.x], ci = ci))

  # 3) Summarise the results
  # Initialize a matrix to store the VaR
  VaR_alpha_emp <- dplyr::bind_rows(VaR_alpha[1,])
  VaR_alpha_emp[1,] <- NA
  if (ES) {
    # Initialize a matrix to store the model's Expected Shortfall
    ES_alpha_mod <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_mod[1,] <- NA
    # Initialize a matrix to store the empiric Expected Shortfall
    ES_alpha_emp <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_emp[1,] <- NA
  }
  # Evaluate Empiric VaR and Expected Shortfalls
  for(i in 1:k){
    # Index of the violations
    idx_violations <- which(Viol_VaR_alpha[,i] == 1)
    if (!purrr::is_empty(idx_violations)){
      VaR_alpha_emp[1,i] <- sum(Viol_VaR_alpha[,i] == 1)/nrow(Viol_VaR_alpha)
      if (ES) {
        ES_alpha_mod[1,i] <- mean(ES_alpha[idx_violations,i])
        ES_alpha_emp[1,i] <- mean(x[idx_violations])
      }
    }
  }
  # Standard output
  output <- structure(
    list(
      VaR = bind_cols(data, VaR_alpha),
      Viol = bind_cols(data, Viol_VaR_alpha),
      VaR_emp = VaR_alpha_emp,
      VaR_test = VaR_tests
    )
  )
  if (ES) {
    output$ES <- bind_cols(data, ES_alpha)
    output$ES_tot <- dplyr::bind_cols(Type = c("Model", "Empiric"),
                                      dplyr::bind_rows(ES_alpha_mod, ES_alpha_emp))
  }
  return(output)
}




