#' Fit a grid of solarModels
#'
#' @param spec Specification
#' @param place Reference location
#' @param arOrder Numeric, maximum AR order.
#' @param maOrder Numeric, maximum MA order.
#' @param archOrder Numeric, maximum ARCH order.
#' @param garchOrder Numeric, maximum GARCH order.
#' @examples
#' spec <- solarModel_spec$new()
#' models <- solarModels_grid(spec, "Bologna", 1,1,1,1)
#' models[which.min(models$L),]
#'
#' @rdname solarModels_grid
#' @name solarModels_grid
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModels_grid <- function(spec, arOrder = 2, maOrder = 2, archOrder = 1, garchOrder = 1, QMLE = FALSE){

  AR_order <- paste0("ARMA(", 1:arOrder)
  MA_order <- paste0(",", 0:maOrder, ")")
  grid <- expand.grid(x = AR_order, y = MA_order)
  for(i in 1:nrow(grid)){
    grid[i,3] <- paste0(grid[i,1], grid[i,2])
  }
  grid_AR_MA <- grid[,3]

  ARCH_order <- paste0("GARCH(", 0:archOrder)
  GARCH_order <- paste0(",", 0:garchOrder, ")")
  grid <- expand.grid(x = ARCH_order, y = GARCH_order)
  for(i in 1:nrow(grid)){
    if (grid[i,1] == "GARCH(0" & grid[i,2] != ",0)"){
      next
    }
    grid[i,3] <- paste0(grid[i,1], grid[i,2])
  }
  grid_ARCH_GARCH <- na.omit(grid[,3])

  grid <- expand.grid(x = grid_AR_MA, y = grid_ARCH_GARCH)
  for(i in 1:nrow(grid)){
    grid[i,3] <- paste0(grid[i,1], "-", grid[i,2])
  }
  grid_models <- grid[,3]
  # Extract model orders
  order <- stringr::str_extract_all(grid_models, "[0-9],[0-9]")
  order <- purrr::map(order, ~as.numeric(unlist(stringr::str_split(.x, ","))))

  models <- list()
  models_name <- c()
  for(i in 1:length(order)){
    print(paste0("Fitting: ", i, "/", length(order)))
    # Specification
    spec$set_mean.model(arOrder = order[[i]][1], maOrder = order[[i]][2])
    if (order[[i]][3] == 0 & order[[i]][4] == 0){
      spec$set_variance.model(archOrder = order[[i]][3], garchOrder = order[[i]][4], garch_variance = FALSE)
    } else {
      spec$set_variance.model(archOrder = order[[i]][3], garchOrder = order[[i]][4], garch_variance = TRUE)
    }
    models_name[i] <- paste0("ARMA(", order[[i]][1], ", ", order[[i]][2],")-GARCH(", order[[i]][3], ", ", order[[i]][4], ")")
    print(models_name[i])
    # Initialize the model
    model <- solarModel$new(spec)
    # Model fit
    model$fit()
    print(model$loglik)
    if (QMLE) {
      model <- solarModel_QMLE(model)
      print(model$loglik)
    }
    # Fit result
    models[[i]] <- model$clone(TRUE)
  }
  names(models) <- models_name
  return(models)
}

#' Compute the AIC and BIC of a solarModel object
#'
#' @param model solarmodel
#'
#' @rdname solarModel_AIC_BIC
#' @name solarModel_AIC_BIC
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_AIC_BIC <- function(model, target = "GHI", type = c("train", "test", "full")){
  moments <- model$moments$conditional
  if (type == "train") {
    moments <-  moments[model$data$isTrain,]
  } else if (type == "test"){
    moments <-  moments[!model$data$isTrain,]
  }
  # Log-likelihood
  L <- model$logLik(moments, target = "GHI")
  n <- sum(!is.infinite(L))
  L <- sum(L[!is.infinite(L)])
  # Model's parameters
  params <- unlist(model$coefficients[-c(1,2)])
  # Remove zero params
  params <- params[params!=0]
  # Number of parameters (minus omega)
  k <- length(params) - 1
  # AIC and BIC
  AIC <- -2 * L + 2 * k
  BIC <- log(n) * k - 2 * L
  # Model name
  model_name <- paste0("ARMA", "(", model$ARMA$order[1], ", ", model$ARMA$order[2],
                       ")-GARCH", "(", model$GARCH$order[1], ", ", model$GARCH$order[2], ")")

  dplyr::tibble(
    Place = model$place,
    Model = model_name,
    ARMA = list(ARMA = model$ARMA$order),
    GARCH = list(GARCH = model$GARCH$order),
    Spec = list(model$spec$clone(TRUE)),
    L = L,
    k = k,
    n = n,
    AIC = AIC,
    BIC = BIC
  )
}

#' Select the Best Model
#'
#' @param spec specification
#'
#' @rdname solarModel_selection
#' @name solarModel_selection
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModels_selection <- function(spec, arOrder = 2, maOrder = 2, archOrder = 1, garchOrder = 1){
  # Logit grid
  spec_logis <- spec$clone(TRUE)
  spec_logis$.__enclos_env__$private$..transform$link <- "logis"
  models_logis <- solarModels_grid(spec_logis, arOrder = arOrder, maOrder = maOrder, archOrder = archOrder, garchOrder = garchOrder)
  # Probit grid
  spec_norm <- spec$clone(TRUE)
  spec_norm$.__enclos_env__$private$..transform$link <- "norm"
  models_norm <- solarModels_grid(spec_norm, arOrder = arOrder, maOrder = maOrder, archOrder = archOrder, garchOrder = garchOrder)
  # Gumbel grid
  spec_invgumb <- spec$clone(TRUE)
  spec_invgumb$.__enclos_env__$private$..transform$link <- "invgumbel"
  models_invgumb <- solarModels_grid(spec_invgumb, arOrder = arOrder, maOrder = maOrder, archOrder = archOrder, garchOrder = garchOrder)
  # AIC / BIC (train)
  train_logis <- purrr::map_df(models_logis, ~solarModel_AIC_BIC(.x, type = "train"))
  train_norm <- purrr::map_df(models_norm, ~solarModel_AIC_BIC(.x, type = "train"))
  train_invgumb <- purrr::map_df(models_invgumb, ~solarModel_AIC_BIC(.x, type = "train"))
  # Select the best model with Train data
  AIC_BIC <- dplyr::bind_rows(
    dplyr::bind_cols(link = "logis", train_logis),
    dplyr::bind_cols(link = "norm", train_norm),
    dplyr::bind_cols(link = "invgumbel", train_invgumb)
  )
  # AIC on train
  best_models_AIC <- dplyr::bind_rows(
    head(AIC_BIC  %>%
           dplyr::filter(link == "logis") %>%
           dplyr::arrange(AIC), n = 1),
    head(AIC_BIC  %>%
           dplyr::filter(link == "norm") %>%
           dplyr::arrange(AIC), n = 1),
    head(AIC_BIC  %>%
           dplyr::filter(link == "invgumbel") %>%
           dplyr::arrange(AIC), n = 1)
  )
  # BIC on train
  best_models_BIC <- dplyr::bind_rows(
    head(AIC_BIC  %>%
           dplyr::filter(link == "logis") %>%
           dplyr::arrange(BIC), n = 1),
    head(AIC_BIC  %>%
           dplyr::filter(link == "norm") %>%
           dplyr::arrange(BIC), n = 1),
    head(AIC_BIC  %>%
           dplyr::filter(link == "invgumbel") %>%
           dplyr::arrange(BIC), n = 1)
  )
  # Best models overall
  smallest_AIC <- head(dplyr::arrange(AIC_BIC, AIC), n = 1)
  smallest_BIC <- head(dplyr::arrange(AIC_BIC, BIC), n = 1)
  # AIC / BIC (test)
  test_logis <- purrr::map_df(models_logis, ~solarModel_AIC_BIC(.x, type = "test"))
  test_norm <- purrr::map_df(models_norm, ~solarModel_AIC_BIC(.x, type = "test"))
  test_invgumb <- purrr::map_df(models_invgumb, ~solarModel_AIC_BIC(.x, type = "test"))
  # Evaluate the choice of the best model on test data
  AIC_BIC_test <- dplyr::bind_rows(
    dplyr::bind_cols(link = "logis", test_logis),
    dplyr::bind_cols(link = "norm", test_norm),
    dplyr::bind_cols(link = "invgumbel", test_invgumb)
  )
  best_models_AIC_test <- dplyr::bind_rows(
    dplyr::filter(AIC_BIC_test, link == best_models_AIC$link[1] & Model == best_models_AIC$Model[1]),
    dplyr::filter(AIC_BIC_test, link == best_models_AIC$link[2] & Model == best_models_AIC$Model[2]),
    dplyr::filter(AIC_BIC_test, link == best_models_AIC$link[3] & Model == best_models_AIC$Model[3])
  )
  best_models_BIC_test <- dplyr::bind_rows(
    dplyr::filter(AIC_BIC_test, link == best_models_BIC$link[1] & Model == best_models_BIC$Model[1]),
    dplyr::filter(AIC_BIC_test, link == best_models_BIC$link[2] & Model == best_models_BIC$Model[2]),
    dplyr::filter(AIC_BIC_test, link == best_models_BIC$link[3] & Model == best_models_BIC$Model[3])
  )

  # Total models
  models <- list(
    logis = models_logis,
    norm = models_norm,
    invgumb = models_invgumb,
    train = list(
      logis = train_logis,
      norm = train_norm,
      invgumb = train_invgumb
    ),
    AIC_BIC_train = AIC_BIC,
    best_AIC = best_models_AIC,
    best_BIC = best_models_BIC,
    test = list(
      logis = test_logis,
      norm = test_norm,
      invgumb = test_invgumb
    ),
    AIC_BIC_test = AIC_BIC_test,
    best_AIC_test = best_models_AIC_test,
    best_BIC_test = best_models_BIC_test,
    best_all_AIC = smallest_AIC,
    best_all_BIC = smallest_BIC
  )
  return(models)
}

