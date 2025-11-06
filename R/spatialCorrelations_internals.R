#' Make a matrix positive semidefined
#'
#' The matrix is decomposed using a spectral decomposition.
#' Then the negative eigenvalues are imputed with `neg_values` and the original matrix is constructed again.
#'
#' @param x symmetric matric
#' @param neg_values numeric
#' @keywords spatialCorrelations internals
#' @noRd
#' @export
makeSemiPositive <- function(x, neg_values = 1e-5){
  mat <- x
  dec <- eigen(x)
  e <- dec$vectors
  lam <- dec$values
  if (any(lam < 0)) {
    lam[lam < 0] <- neg_values
    mat <- e %*% diag(lam) %*% t(e)
  }
  attr(mat, "index_neg_values") <- which(dec$values < 0)
  attr(mat, "original_values") <- dec$values[dec$values < 0]
  return(mat)
}

#' Check admissibility of common probabilities matrix
#'
#' @param commonprob symmetric matrix of joint probabilities.
#' @param check.commonprob logical
#' @param simulvals `bindata::SimulVals`
#' @keywords spatialCorrelations internals
#' @noRd
#' @export
interpolate_commonprob <- function(commonprob, check.commonprob = TRUE, simulvals = NULL){

  if (is.null(simulvals)) {
    simulvals <- bindata::SimulVals
  }
  # Check marginal probabilities
  if (check.commonprob) {
    check <- bindata::check.commonprob(commonprob)
    if (!check) {
      cat(attr(check, "message"), sep = "\n")
      message("Matrix commonprob not admissible. Run check.commonprob(commonprob) for more details.")
    }
  }
  margprob <- diag(commonprob)
  for (m in 1:(ncol(commonprob) - 1)) {
    for (n in (m + 1):nrow(commonprob)) {
      x <- cbind(margprob[m], margprob[n], as.numeric(dimnames(simulvals)[[3]]))
      y <- e1071::interpolate(x, simulvals)
      # Check admissibility and impute errors
      if (commonprob[m, n] > max(y)) {
        print(paste0("Error: (", n, " ", m, ")", " substitued from ", commonprob[m, n], " to ", max(y)))
        commonprob[m, n] <- max(y)
      } else if (commonprob[m, n] < min(y)){
        print(paste0("Error: (", n, " ", m, ")", " substitued from ", commonprob[m, n], " to ", min(y)))
        commonprob[m, n] <- min(y)
      }
    }
    message(m, "/", (ncol(commonprob) - 1), "\r", appendLF = FALSE)
  }
  # Check if any elements is NAs
  if (any(is.na(commonprob))) {
    stop("Some elements in commonprob are NAs... margprob and commonprob are not compatible?")
  }
  return(commonprob)
}

#' spatialCorrelation object
#'
#' @noRd
#' @export
spatialCorrelation_mixture <- function(models, nmonths){

  if (missing(nmonths)){
    filter_data <- function(data) dplyr::filter(data, isTrain & weights != 0)
  } else {
    filter_data <- function(data) dplyr::filter(data, Month %in% nmonths & isTrain & weights != 0)
  }

  # Number of models
  n_models <- length(models)
  # Extract a dataset containing all U and B variables
  data <- filter_data(models[[1]]$data)
  # Rename the variables
  df_u <- dplyr::tibble(date = data$date, u1 = data$u_tilde, B1 = data$B, x1_1 = ((u1-data$mu1)/data$sd1)*B1, x2_1 = ((u1-data$mu2)/data$sd2)*(1-B1))
  df_u <- dplyr::select(df_u, -B1)
  # Create a unique dataset
  for(i in 2:n_models){
    message(i, "/", n_models, "\r", appendLF = FALSE)
    data <- filter_data(models[[i]]$data)
    data <- dplyr::tibble(date = data$date, u1 = data$u_tilde, B1 = data$B, x1 = ((u1-data$mu1)/data$sd1)*B1, x2 = ((u1-data$mu2)/data$sd2)*(1-B1))
    data <- dplyr::select(data, -B1)
    colnames(data) <- c("date", paste0(c("u", "x1_", "x2_"), i))
    df_u <- dplyr::left_join(df_u, data, by = "date")
  }
  # Compute correlations
  message("Computing correlations... ")
  if(missing(nmonths)){
    df_u <- na.omit(df_u)[,-1]
    # Split the datasets
    data_12 <- dplyr::select(df_u, dplyr::contains("x"))
    message("Computing mixture correlations: X1-X1 | X2-X2 | X1-X2 | X2-X1 |")
    # Store the results
    cr_list <- cor(data_12)
  } else {
    m <- 1
    cr_list <- list()
    df_u$Month <- lubridate::month(df_u$date)
    df_u <- na.omit(df_u)[,-1]
    for(m in nmonths){
      message("Computing correlation for Month: ", m, "...")
      # Split the datasets
      data_12 <- dplyr::select(dplyr::filter(df_u, Month == m), dplyr::contains("x"))
      message("Computing mixture correlations: X1-X1 | X2-X2 | X1-X2 | X2-X1 |")
      # Store the results
      cr_list[[m]] <- cor(data_12)
    }
    names(cr_list) <- lubridate::month(nmonths, label = TRUE)
  }
  # Unique names as models names
  return(cr_list)
}

#' spatialCorrelation object
#'
#' @noRd
#' @export
spatialCorrelation_binprobs <- function(models, nmonths = 1:12){
  # Number of models
  n_models <- length(models)
  # Extract a dataset containing all U and B variables
  data <- models[[1]]$data
  data <- dplyr::filter(data, Month %in% nmonths & isTrain & weights != 0)
  # Rename the variables
  df_B <- dplyr::tibble(date = data$date, u1 = data$u_tilde, B1 = data$B)
  # Create a unique dataset
  for(i in 2:n_models){
    message(i, "/", n_models, "\r", appendLF = FALSE)
    data <- models[[i]]$data
    data <- dplyr::filter(data, Month %in% nmonths & isTrain & weights != 0)
    data <- dplyr::tibble(date = data$date, u = data$u_tilde, B = data$B)
    colnames(data) <- c("date", paste0(c("u", "B"), i))
    df_B <- dplyr::left_join(df_B, data, by = "date")
  }
  df_B <- na.omit(df_B)
  # Add month variable
  df_B$Month <- lubridate::month(df_B$date)
  # Compute correlations
  cr_list <- list()
  cr_list$sigma <- cr_list$commonprob <- cr_list$cr_B <- list()
  m <- 1
  for(m in nmonths){
    message("Computing correlation for Month: ", m, "...")
    # Split the datasets
    data_B <- dplyr::select(dplyr::filter(df_B, Month == m), dplyr::contains("B"))
    message("Computing mixture correlations...", appendLF = FALSE)
    cr_list$cr_B[[m]] <- cor(data_B)
    message("Done!", appendLF = TRUE)
    message("Computing marginal probabilities B...", appendLF = FALSE)
    # Marginal probabilities
    cr_list$commonprob[[m]] <- diag(apply(data_B, 2, mean, na.rm = TRUE))
    message("Done!", appendLF = TRUE)
    message("Computing joint probabilities of event [Bi = 1 & Bj = 1]...")
    # Common probabilities
    for(i in 1:(n_models-1)){
      for(j in (i+1):n_models){
        pairwise_probs_ij <-  mean(data_B[,i] == 1 & data_B[,j] == 1, na.rm = TRUE)
        cr_list$commonprob[[m]][i,j] <- cr_list$commonprob[[m]][j,i] <- pairwise_probs_ij
      }
      if (i == ncol(data_B)) {
        message("Done!")
      } else {
        message(i, "/", ncol(data_B), "\r", appendLF = FALSE)
      }
    }
    message("Adjusting commonprobs...")
    # Convert commonprob into a matrix sigma to improve simulations
    cr_list$commonprob[[m]] <- interpolate_commonprob(cr_list$commonprob[[m]], check.commonprob = FALSE)
    message("Done!", appendLF = TRUE)
    message("Computing implied sigma from commonprob...", appendLF = FALSE)
    # Adjust implied matrix sigma
    cr_list$sigma[[m]] <- bindata::commonprob2sigma(cr_list$commonprob[[m]])
    cr_list$sigma[[m]] <- makeSemiPositive(cr_list$sigma[[m]], neg_values = 1e-3)
    message("Done!", appendLF = TRUE)
    # Unique names as models names
    colnames(cr_list$cr_B[[m]]) <- rownames(cr_list$cr_B[[m]]) <- names(models)
    colnames(cr_list$commonprob[[m]]) <- rownames(cr_list$commonprob[[m]]) <- names(models)
    colnames(cr_list$sigma[[m]]) <- rownames(cr_list$sigma[[m]]) <- names(models)
  }
  names(cr_list$commonprob) <- lubridate::month(nmonths, label = TRUE)
  names(cr_list$cr_B) <- lubridate::month(nmonths, label = TRUE)
  names(cr_list$sigma) <- lubridate::month(nmonths, label = TRUE)
  return(cr_list)
}

