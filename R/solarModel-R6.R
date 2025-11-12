#' Solar Model in R6 Class
#'
#' The `solarModel` class is an implementation of a comprehensive solar model that includes fitting seasonal models,
#' detecting outliers, performing transformations, and applying time-series models such as AR and GARCH. This model
#' is specifically designed to predict solar radiation data, and it uses seasonal and Gaussian Mixture models to
#' capture the underlying data behavior.
#'
#' @description
#' The `solarModel` class allows for the step-by-step fitting and transformation of solar radiation data, from clear sky
#' models to GARCH models for residual analysis. It utilizes various private and public methods to fit the seasonal
#' clearsky model, compute risk drivers, detect outliers, and apply time-series models.
#'
#' @examples
#'
#' # Model specification
#' spec <- solarModel_spec$new()
#' spec$set_mean.model(arOrder = 1, maOrder = 1)
#' spec$specification("Bologna")
#' spec
#' # Model fit
#' Bologna <- solarModel$new(spec)
#' Bologna$fit()
#' Bologna
#' # save(spec, file = "data/Bologna.RData")
#'
#' # Extract and update the parameters
#' model <- Bologna$clone(TRUE)
#' params <- model$coefficients
#' model$update(params)
#' model$filter()
#'
#' # Fit a model with the realized clear sky
#' spec$control$stochastic_clearsky <- TRUE
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' # Model fit
#' model$fit()
#'
#' # Fit a model for the clearsky
#' spec_Ct <- spec
#' spec_Ct$control$stochastic_clearsky <- FALSE
#' spec_Ct$target <- "clearsky"
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' # Model fit
#' model$fit()
#'
#' @rdname solarModel
#' @name solarModel
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
#self <- Bologna$.__enclos_env__$self
#private <- Bologna$.__enclos_env__$private
solarModel <- R6::R6Class("solarModel",
                          # ====================================================================================================== #
                          #                                             Public slots
                          # ====================================================================================================== #
                          public = list(
                            #' @description
                            #' Initialize a `solarModel`
                            #' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
                            initialize = function(spec){
                              # Seasonal data by month and day for an year with 366 days
                              dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "1 day")
                              seasonal_data <- dplyr::tibble(date = dates)
                              seasonal_data <- dplyr::mutate(seasonal_data,
                                                             Month = lubridate::month(date),
                                                             Day = lubridate::day(date),
                                                             n = number_of_day(date))
                              seasonal_data <- dplyr::select(seasonal_data, -date)
                              seasonal_data <- dplyr::arrange(seasonal_data, Month, Day)
                              # **************************************************** #
                              # Private components
                              private$..data <- dplyr::select(spec$data, -H0)
                              private$..data$loglik <- 0
                              private$..seasonal_data <- seasonal_data
                              private$..monthly_data <- dplyr::tibble(Month = 1:12, Yt_tilde_uncond = 0, sigma_uncond = 1)
                              private$..loglik <- NA
                              private$..spec <- spec$clone(TRUE)
                              private$..transform <- solarTransform$new(alpha = 0, beta = 1, link = private$..spec$transform$link)
                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Initialize and fit a \code{\link{solarModel}} object given the specification contained in `$control`.
                            fit = function(){
                              # 1) Clearsky
                              self$fit_seasonal_model_Ct()
                              # 2) Risk-driver
                              self$compute_risk_drivers()
                              # 3) Solar transform
                              self$fit_transform()
                              # 4) Seasonal mean
                              self$fit_seasonal_model_Yt()
                              #    Center the mean to be exactly equal to zero
                              self$fit_monthly_mean()
                              # 5) ARMA model
                              self$fit_ARMA()
                              # 6) Seasonal variance
                              self$fit_seasonal_variance()
                              # 7) GARCH variance
                              self$fit_GARCH()
                              # 8) Mixture model
                              self$fit_NM_model()
                              self$update_NM_classification()
                              # Update the moments
                              self$update_moments()
                              # Update the log-likelihoods
                              self$update_logLik()
                            },
                            #' @description
                            #' Initialize and fit a \code{\link{seasonalClearsky}} model given the specification contained in `$control`.
                            fit_seasonal_model_Ct = function(){
                              # Control
                              control <- self$spec$clearsky
                              # Arguments
                              data <- dplyr::filter(private$..data, isTrain & weights != 0)
                              # **************************************************** #
                              # Initialize a seasonal model for clear sky radiation
                              seasonal_model_Ct <- seasonalClearsky$new(control = control)
                              # Fit the seasonal model
                              seasonal_model_Ct$fit(data[[self$spec$target]],
                                                    data[["date"]],
                                                    self$spec$coords$lat,
                                                    data[["clearsky"]])
                              # **************************************************** #
                              # Private components
                              # Add seasonal clear sky max to seasonal data
                              private$..data[["Ct"]] <- seasonal_model_Ct$predict(newdata = self$data)
                              # Store the model
                              private$..seasonal_model_Ct <- seasonal_model_Ct$clone(deep = TRUE)
                            },
                            #' @description
                            #' Compute the risk drivers and impute the observation that are greater or equal to the clear sky level.
                            compute_risk_drivers = function(){
                              # Arguments
                              target <- self$spec$target
                              control <- self$spec
                              data <- self$data
                              transform <- self$transform
                              # **************************************************** #
                              if (control$stochastic_clearsky) {
                                # Risk driver
                                data$Xt <- transform$X(data[[target]], data$clearsky)
                                # Detect and impute outliers
                                outliers <- clearsky_outliers(data[["GHI"]], data$clearsky, date = data$date, quiet = control$quiet)
                                # Update GHI
                                data[["GHI"]] <- outliers$x
                                # Risk driver
                                data[["Xt"]] <- transform$X(data[["GHI"]], data$clearsky)
                              } else {
                                # Detect and impute outliers
                                outliers <- clearsky_outliers(data[["GHI"]], data$Ct, date = data$date, quiet = control$quiet)
                                # Update GHI
                                data[["GHI"]] <- outliers$x
                                # Add computed risk driver in data
                                data[["Xt"]] <- transform$X(data[["GHI"]], data$Ct)
                              }
                              # **************************************************** #
                              # Add computed risk driver in data
                              private$..data[["Xt"]] <- data$Xt
                              # Store outliers data
                              private$outliers <- outliers
                              # private$..data$weights[outliers$index] <- 0
                            },
                            #' @description
                            #' Fit the parameters of the \code{\link{solarTransform}} object.
                            fit_transform = function(){
                              # Arguments
                              control <- self$spec$transform
                              data <- dplyr::filter(private$..data, isTrain & weights != 0)
                              outliers <- private$outliers
                              # **************************************************** #
                              # Transformation parameters
                              params <- self$transform$fit(data$Xt, control$threshold, control$min_pos, control$max_pos)
                              # Update transform parameters
                              private$..transform$update(params$alpha, params$beta)
                              data <- private$..data
                              # Rescale the minimum to avoid that extreme values
                              # breaks the ARMA-GARCH routines (especially the minimum)
                              idx_Xt_min <- which(data$Xt <= params$Xt_min)
                              data$Xt[idx_Xt_min] <- params$Xt_min * (1 + control$delta)
                              # Rescale also the maximum
                              idx_Xt_max <- which(data$Xt >= params$Xt_max)
                              data$Xt[idx_Xt_max] <- params$Xt_max * (1 - control$delta)
                              # Store the index of the imputed values
                              outliers$index_type$transform <- c(idx_Xt_min, idx_Xt_max)
                              # Update outliers index and dates
                              outliers$index <- unique(c(outliers$index, outliers$index_type$transform))
                              outliers$date <- data$date[outliers$index]
                              # **************************************************** #
                              # Compute the transformed variable
                              data$Yt <- self$transform$Y(self$transform$X_prime(data[["Xt"]]))
                              # Add Yt to private data
                              private$..data[["Yt"]] <- data$Yt
                              # Store the updated outliers
                              private$outliers <- outliers
                            },
                            #' @description
                            #' Fit a \code{\link{seasonalModel}} the transformed variable (`Yt`) and compute deseasonalized series (`Yt_tilde`).
                            fit_seasonal_model_Yt = function(){
                              # Arguments
                              target <- self$spec$target
                              control <- self$spec$seasonal.mean
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # Custom formula
                              formula_Yt <- ifelse(control$include.intercept, "Yt ~ 1", "Yt ~ -1")
                              formula_Yt <- ifelse(control$include.trend, paste0(formula_Yt, " + t"), formula_Yt)
                              # Initialize the seasonal model for Yt
                              seasonal_model_Yt <- seasonalModel$new(order = control$order, period = control$period)
                              seasonal_model_Yt$fit(as.formula(formula_Yt), data = data_train)
                              # Compute Yt_bar
                              data$Yt_bar <- seasonal_model_Yt$predict(newdata = data)
                              # Compute Yt_tilde
                              data$Yt_tilde <- data$Yt - data$Yt_bar
                              # Update parameters names
                              coefs_names <- c()
                              orig_names <- seasonal_model_Yt$model$coefficients_names
                              if(control$include.intercept) {
                                coefs_names <- c("a_0")
                                orig_names <- orig_names[-c(1)]
                              }
                              if (control$include.trend) {
                                coefs_names <- c(coefs_names, "t")
                                orig_names <- orig_names[-c(1)]
                               }
                              if (length(control$order) > 1 || control$order > 0) {
                                coefs_names <- c(coefs_names, paste0("a_", orig_names))
                              }
                              # Update parameters name inside the R6 object
                              seasonal_model_Yt$.__enclos_env__$private$..model$coefficients_names <- coefs_names
                              names(seasonal_model_Yt$.__enclos_env__$private$..std.errors) <- coefs_names
                              # **************************************************** #
                              # Private components
                              # Store seasonal model for Yt
                              private$..seasonal_model_Yt <- seasonal_model_Yt$clone(deep = TRUE)
                              # Add Yt_tilde to data
                              private$..data[["Yt_tilde"]] <- data$Yt_tilde
                              # Add Yt_bar to data
                              private$..data[["Yt_bar"]] <- self$seasonal_model_Yt$predict(newdata = data)
                              # Add GHI_bar to data
                              private$..data[[paste0(target, "_bar")]] <- self$transform$iRY(self$data$Yt_bar, self$data$Ct)
                            },
                            #' @description
                            #' Correct the deseasonalized series (`Yt_tilde`) by subtracting its monthly mean (`Yt_tilde_uncond`).
                            fit_monthly_mean = function(){
                              control <- self$spec
                              # Unconditional mean for Yt_tilde
                              if (control$seasonal.mean$monthly.mean) {
                                data <- private$..data
                                # Train data
                                train_data <- dplyr::filter(data, isTrain & weights != 0)
                                # Compute monthly unconditional mean
                                monthly_mean <- train_data %>%
                                  dplyr::group_by(Month) %>%
                                  dplyr::summarise(Yt_tilde_uncond = mean(Yt_tilde)) %>%
                                  dplyr::ungroup()
                                # Add unconditional mean to the dataset
                                data <- dplyr::left_join(data, monthly_mean, by = c("Month"))
                                # Update Yt_tilde
                                data$Yt_tilde <- data$Yt_tilde - data$Yt_tilde_uncond
                                # **************************************************** #
                                # Private components
                                # Add unconditional mean to monthly data
                                private$..monthly_data[["Yt_tilde_uncond"]] <- monthly_mean$Yt_tilde_uncond
                                # Add Yt_tilde to data
                                private$..data[["Yt_tilde"]] <- data$Yt_tilde
                              }
                            },
                            #' @description
                            #' Fit an AR model (`Yt_tilde`) and compute AR residuals (`eps`).
                            fit_ARMA = function(){
                              # Arguments
                              control <- self$spec$mean.model
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain)# & weights != 0)
                              # Initialize an AR model
                              ARMA <- ARMA_modelR6$new(control$arOrder, control$maOrder, control$include.intercept)
                              # Fitted AR model
                              ARMA$fit(data_train$Yt_tilde)
                              # Fitted Yt_tilde
                              data$Yt_tilde_hat <- ARMA$filter(data$Yt_tilde)
                              # Fitted residuals
                              data$eps <- data$Yt_tilde - data$Yt_tilde_hat
                              # **************************************************** #
                              # Private components
                              # Store ARMA model
                              private$..ARMA <- ARMA$clone(TRUE)
                              # Add fitted Yt_tilde
                              private$..data[["Yt_tilde_hat"]] <- data$Yt_tilde_hat
                              # Add fitted residuals
                              private$..data[["eps"]] <- data$eps
                            },
                            #' @description
                            #' Fit a \code{\link{seasonalModel}} on AR squared residuals (`eps`) and compute deseasonalized residuals `eps_tilde`.
                            fit_seasonal_variance = function(){
                              # Control
                              control <- self$spec$seasonal.variance
                              # Dataset
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # Squared residuals
                              data_train$eps2 <- data_train$eps^2
                              # Custom formula
                              formula <- ifelse(control$include.trend, "eps2 ~ 1 + t", "eps2 ~ 1")
                              # Initialize the seasonal model for eps2
                              seasonal_variance <- seasonalModel$new(order = control$order, period = control$period)
                              seasonal_variance$fit(formula = as.formula(formula), data = data_train)
                              # Fitted seasonal standard deviation
                              private$..data[["sigma_bar"]] <- sqrt(seasonal_variance$predict(newdata = data))
                              # Compute standardized residuals
                              private$..data[["eps_tilde"]] <- private$..data[["eps"]] / private$..data[["sigma_bar"]]
                              # Update parameters names
                              coefs_names <- c("c_0")
                              orig_names <- seasonal_variance$model$coefficients_names[-1]
                              if (control$include.trend) {
                                coefs_names <- c(coefs_names, "t")
                                orig_names <- orig_names[-c(1)]
                              }
                              if (control$order > 0) {
                                coefs_names <- c(coefs_names, paste0("c_", orig_names))
                              }
                              # Update parameters name inside the R6 object
                              seasonal_variance$.__enclos_env__$private$..model$coefficients_names <- coefs_names
                              names(seasonal_variance$.__enclos_env__$private$..std.errors) <- coefs_names
                              # **************************************************** #
                              # Store seasonal variance model
                              private$..seasonal_variance <- seasonal_variance$clone(deep = TRUE)
                              # **************************************************** #
                              # Compute monthly corrective variance to ensure unitary monthly variance
                              self$fit_monthly_variance()
                              # Correction of the parameters to ensure unitary variance
                              self$correct_seasonal_variance()
                            },
                            #' @description
                            #' Correct the standardized series (`eps_tilde`) by subtracting its monthly mean (`sigma_uncond`).
                            fit_monthly_variance = function(){
                              # Arguments
                              control <- self$spec$seasonal.variance
                              data <- private$..data
                              # **************************************************** #
                              # Unconditional variance for eps_tilde
                              if (control$monthly.mean) {
                                # Train data
                                train_data <- dplyr::filter(data, isTrain & weights != 0)
                                # Compute monthly unconditional mean
                                monthly_data <- train_data %>%
                                  dplyr::group_by(Month) %>%
                                  dplyr::summarise(sigma_uncond = sd(eps_tilde, na.rm = TRUE)) %>%
                                  dplyr::ungroup()
                                # Add unconditional variance to monthly data
                                private$..monthly_data[["sigma_uncond"]] <- monthly_data$sigma_uncond
                                # Updated standardized residuals
                                private$..data[["eps_tilde"]] <- private$..data[["eps"]] / (self$data$sigma_bar * self$data$sigma_uncond)
                              }
                            },
                            #' @description
                            #' Correct the parameters of the seasonal variance to ensure a unitary variance
                            correct_seasonal_variance = function(){
                              # Control
                              control <- self$spec$seasonal.variance
                              # Dataset with all the parameters
                              data <- self$data
                              # **************************************************** #
                              # Initialize a slot to store monthly correction
                              private$..seasonal_variance$extra_params$correction <- 1
                              # Correction to ensure unitary variance
                              if (control$correction) {
                                seasonal_variance <- private$..seasonal_variance$clone(TRUE)
                                # Update train data
                                data_train <- dplyr::filter(data, isTrain & weights != 0)
                                # Store correction factor
                                seasonal_variance$extra_params$correction <- var(data_train$eps / (data_train$sigma_bar * data_train$sigma_uncond))
                                # Correct parameters to ensure unitary variance
                                std.errors <- seasonal_variance$std.errors
                                seasonal_variance$update(seasonal_variance$coefficients * seasonal_variance$extra_params$correction)
                                seasonal_variance$update_std.errors(seasonal_variance$extra_params$correction * std.errors)
                                # Update sigma_bar into the seasonal data
                                private$..data[["sigma_bar"]] <- sqrt(seasonal_variance$predict(newdata = private$..data))
                                # Updated standardized residuals
                                private$..data[["eps_tilde"]] <- private$..data[["eps"]] / (private$..data[["sigma_bar"]] * self$data$sigma_uncond)
                                # Update seasonal variance model
                                private$..seasonal_variance <- seasonal_variance$clone(TRUE)
                              }
                            },
                            #' @description
                            #' Fit a `GARCH` model on the deseasonalized residuals (`eps_tilde`).
                            #' Compute the standardized (`u`) and monthly deseasonalized residuals (`u_tilde`).
                            fit_GARCH = function(){
                              # Arguments
                              control <- self$spec
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # GARCH specification
                              GARCH_spec <- control$variance.model
                              # Initialize a GARCH model
                              GARCH_model <- sGARCH$new(GARCH_spec$archOrder, GARCH_spec$garchOrder)
                              # Control for garch variance
                              if (control$garch_variance) {
                                # Fit the model
                                GARCH_model$fit(data_train$eps_tilde, data_train$weights)
                                # Fitted std. deviation on true eps_tilde
                                data$sigma <- sqrt(GARCH_model$filter(data$eps_tilde))
                                # Fitted standardized residuals
                                data$u_tilde <- data$eps_tilde / data$sigma
                              } else {
                                # Fitted std. deviation
                                data$sigma <- 1
                                # Not compute standardized residuals
                                data$u_tilde <- data$eps_tilde
                              }
                              # **************************************************** #
                              # Store GARCH parameters
                              private$..GARCH <- GARCH_model$clone(TRUE)
                              # Fitted GARCH std. deviation residuals
                              private$..data[["sigma"]] <- data$sigma
                              # Fitted standardized residuals
                              private$..data[["u_tilde"]] <- data$u_tilde
                            },
                            #' @description
                            #' Initialize and fit a `solarMixture` object.
                            fit_NM_model = function(){
                              # Arguments
                              control <- self$spec$mixture.model
                              outliers <- private$outliers
                              # Train data
                              data_train <- dplyr::filter(private$..data, isTrain & weights != 0)
                              # **************************************************** #
                              # Gaussian Mixture fitted only on train data
                              NM_model <- solarMixture$new(components = 2, abstol = control$abstol,
                                                           maxit = control$maxit, maxrestarts = control$maxrestarts)
                              # Match moments
                              if (control$match.expectation | control$match.variance) {
                                # Target empirical moments
                                target <- data_train %>%
                                  group_by(Month) %>%
                                  summarise(mu_target = mean(u_tilde),
                                            var_target = var(u_tilde))
                                # Match or not empirical moments
                                if (!control$match.empiric) {
                                  target$mu_target <- 0
                                  target$var_target <- 1
                                }
                                # Match or not expectation
                                if (!control$match.expectation) {
                                  target$mu_target <- NA
                                }
                                # Match or not variance
                                if (!control$match.variance) {
                                  target$var_target <- NA
                                }
                                NM_model$fit(x = data_train$u_tilde, date = data_train$date, weights = data_train$weights,
                                             mu_target = target$mu_target, var_target = target$var_target,
                                             method = control$method)

                              } else {
                                NM_model$fit(x = data_train$u_tilde, date = data_train$date, weights = data_train$weights,
                                             method = control$method)
                              }
                              # **************************************************** #
                              # Private components
                              # Store Gaussian Mixture parameters
                              private$..NM_model <- NM_model$clone(TRUE)
                            },
                            #' @description
                            #' Update the parameters inside object
                            #' @param params List of parameters. See the slot `$coefficients` for a template.
                            update = function(params){
                              if (!missing(params)) {
                                if (!is.list(params)) {
                                  params <- solarModel_match_params(params, self$coefficients)
                                }
                                # Update transform parameters
                                private$..transform$update(params$params$alpha, params$params$beta)
                                # Update clear sky model
                                private$..seasonal_model_Ct$update(unlist(params$seasonal_model_Ct))
                                # Update seasonal mean model
                                private$..seasonal_model_Yt$update(unlist(params$seasonal_model_Yt))
                                # Update AR mean model
                                private$..ARMA$update(unlist(params$ARMA))
                                # Update seasonal variance model
                                private$..seasonal_variance$update(unlist(params$seasonal_variance))
                                # Update GARCH variance model
                                private$..GARCH$update(unlist(params$GARCH))
                                # ***************** Update Gaussian Mixture model *****************
                                means <- cbind(mu1 = unlist(params$NM_mu_up), mu2 = unlist(params$NM_mu_dw))
                                sd <- cbind(sd1 = unlist(params$NM_sd_up), sd2 = unlist(params$NM_sd_dw))
                                p <- cbind(p1 = unlist(params$NM_p_up), p2 = 1 - unlist(params$NM_p_up))
                                private$..NM_model$update(means = means, sd = sd, p = p)
                                # Set Log-likelihood to NA
                                private$..loglik <- NA
                              } else {
                                message("`params` is missing nothing to update!")
                              }
                            },
                            #' @description
                            #' Update the moments inside object
                            update_moments = function(){
                              # Update conditional moments
                              private$..moments$conditional <- solarMoments_conditional(self$data, control_model = self$spec)
                              # Update unconditional moments
                              private$..moments$unconditional <- solarMoments_unconditional(self$data, ARMA = self$ARMA, GARCH = self$GARCH)
                            },
                            #' @description
                            #' Update the log-likelihood inside object
                            update_logLik = function(){
                              # Update the log-likelihoods
                              private$..data[["loglik"]] <- self$logLik()
                              # Update total log-likelihood
                              private$..loglik <- sum(self$data$loglik)
                            },
                            #' @description
                            #' Update the clear sky and risk drivers
                            update_risk_drivers = function(){
                              # Update clear sky
                              private$..data[["Ct"]] <- self$seasonal_model_Ct$predict(newdata = self$data)
                              # Update risk driver
                              self$compute_risk_drivers()
                              # Update solar transform and compute Yt
                              self$fit_transform()
                            },
                            #' @description
                            #' Update the classification of the Bernoulli random variable.
                            #' @param filter Logical, when `TRUE` before the classification will be runned the command
                            #' ` self$NM_model$filter()` to update the mixture classification.
                            update_NM_classification = function(filter = FALSE){
                              if (filter) {
                                # Update mixture classification
                                private$..NM_model$filter()
                              }
                              # Complete data
                              data <- private$..data
                              # **************************************************** #
                              # Ensure that no columns with B name are included
                              data <- data[, !(colnames(data) %in% c("B", "z1","z2"))]
                              # Classify the series
                              df_m <- list()
                              for(nmonth in 1:12){
                                df_m[[nmonth]] <- dplyr::filter(data, Month == nmonth)
                                df_nm <- private$..NM_model$model[[nmonth]]$classify(df_m[[nmonth]]$u_tilde)
                                df_m[[nmonth]]$B <- df_nm$B1
                                df_m[[nmonth]]$z1 <- df_nm$z1
                                df_m[[nmonth]]$z2 <- df_nm$z2
                              }
                              # Add the data
                              df_m <- dplyr::select(dplyr::bind_rows(df_m), date, B, z1, z2)
                              data <- dplyr::left_join(data, df_m, by = "date")
                              # **************************************************** #
                              # Private components
                              # Update fitted series of Bernoulli
                              private$..data[["B"]] <- data$B
                              # Update fitted series of Gaussians components
                              private$..data[["z1"]] <- data$z1
                              private$..data[["z2"]] <- data$z2
                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Filter the time series when new parameters are supplied in the method `$update(params)`.
                            #' @param fit Logical, when `TRUE`, if in the model's specification, the monthly mean and variances will be re estimated and the seasonal variance corrected
                            #' such that the total variance of the deseasonalized residuals is zero.
                            #' @return Update the slots `$data`, `$seasonal_data`, `$monthly_data`
                            filter = function(fit = TRUE){
                              # Arguments
                              control <- self$spec
                              target <- self$spec$target
                              # **************************************************** #
                              # Update seasonal mean of Yt
                              private$..data[["Yt_bar"]] <- self$seasonal_model_Yt$predict(newdata = private$..data)
                              # Update seasonal mean of target variable
                              private$..data[[paste0(target, "_bar")]] <- self$transform$iRY(private$..data[["Yt_bar"]], self$data[["Ct"]])
                              # Update Yt_tilde
                              private$..data[["Yt_tilde"]] <- private$..data[["Yt"]] - private$..data$Yt_bar
                              # Fit the corrective mean
                              if (fit) {
                                self$fit_monthly_mean()
                              }
                              # Update Yt_tilde_hat
                              private$..data[["Yt_tilde_hat"]] <- self$ARMA$filter(private$..data[["Yt_tilde"]])
                              # Update ARMA residuals
                              private$..data[["eps"]] <- private$..data[["Yt_tilde"]] - private$..data[["Yt_tilde_hat"]]
                              # **************************************************** #
                              # Update seasonal std. deviation
                              private$..data[["sigma_bar"]] <- sqrt(self$seasonal_variance$predict(newdata = private$..data))
                              # Update standardized residuals
                              private$..data[["eps_tilde"]] <- private$..data[["eps"]] / private$..data[["sigma_bar"]]
                              if (fit) {
                                # Compute corrective monthly variance
                                self$fit_monthly_variance()
                                # Correct the seasonal variance
                                self$correct_seasonal_variance()
                              }
                              # **************************************************** #
                              # Update Garch standard deviation
                              if (self$spec$garch_variance) {
                                private$..data[["sigma"]] <- sqrt(self$GARCH$filter(private$..data[["eps_tilde"]]))
                              } else {
                                private$..data[["sigma"]] <- 1
                              }
                              # Fitted standardized residuals
                              private$..data[["u_tilde"]] <- private$..data[["eps_tilde"]] / private$..data[["sigma"]]
                            },
                            #' @description
                            #' Compute the conditional moments
                            #' @param t_now Character date. Today date.
                            #' @param t_hor Character date. Horizon date.
                            #' @param theta Numeric, shift parameter for the mixture.
                            #' @param quiet Logical for verbose messages.
                            Moments = function(t_now, t_hor, theta = 0, quiet = FALSE){
                              purrr::map_df(t_hor, ~solarMoments(t_now, .x, self$data, self$ARMA,
                                                                       self$GARCH, self$NM_model, self$transform, theta = theta, quiet = quiet))
                            },
                            #' @description
                            #' Value at Risk for a `solarModel`
                            #' @param moments moments dataset
                            #' @param t_now Character date. Today date.
                            #' @param t_hor Character date. Horizon date.
                            #' @param theta Numeric, shift parameter for the mixture.
                            #' @param ci Confidence interval (one tail).
                            VaR = function(moments, t_now, t_hor, theta = 0, ci = 0.05){
                              # Moments
                              if (missing(moments)) {
                                t_seq <- seq.Date(as.Date(t_now), as.Date(t_hor), 1)[-1]
                                moments <- purrr::map_df(t_seq, ~self$Moments(t_now, .x, theta, quiet = FALSE))
                              }
                              # Add realized GHI
                              moments <- dplyr::left_join(moments, self$data[,c("date", "GHI")], by = "date")
                              # Initialize the VaR
                              moments$VaR <- NA
                              for(i in 1:nrow(moments)){
                                # Moments
                                df_n <- moments[i,]
                                # Quantile of Yt
                                cdf_Y <- function(x) pmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
                                # Value at Risk with probability `ci`
                                # moments$VaR[i] <- qsolarGHI(ci, df_n$Ct, df_n$alpha, df_n$beta, cdf_Yt)
                                moments$VaR[i] <- qsolarGHI(ci, df_n$Ct, df_n$alpha, df_n$beta, cdf_Y, link = self$spec$transform$link)
                              }
                              # Violations of VaR
                              moments$et <- ifelse(moments$GHI < moments$VaR, 1, 0)
                              # Select only relevant variables
                              moments <- dplyr::select(moments, date, Year, Month, Day, GHI, VaR, et)
                              return(moments)
                            },
                            #' @description
                            #' Compute the log-likelihood of the model and update the slot `$loglik`.
                            #' @param moments Dataset containing the moments to use for computation.
                            #' @param target Character. Target variable to use "Yt" or "GHI".
                            #' @param quasi Logical, when `TRUE` is computed the pseudo-likelihood with Gaussian link.
                            logLik = function(moments, target = "Yt", quasi = FALSE){
                              # Default argument
                              if (missing(moments)) {
                                moments <- self$moments$conditional
                              }
                              # Add Yt and weights
                              moments <- dplyr::left_join(moments, private$..data[, c("date", target, "weights", "isTrain")], by = "date")
                              # Compute the quasi and true log-likelihood
                              moments$loglik <- 0
                              # Density: solar radiation
                              if (target == "GHI") {
                                for(i in 1:nrow(moments)){
                                  df_n <- moments[i,]

                                  if (df_n$weights == 0 & df_n$isTrain) {
                                    moments$loglik[i] <- 0
                                    next
                                  }
                                  if (quasi) {
                                    # Normal density
                                    pdf_quasi <- function(x) dnorm(x, df_n$e_Yt, df_n$sd_Yt)
                                    # Quasi Log-likelihood
                                    moments$loglik[i] <- log(dsolarGHI(df_n$GHI, df_n$Ct, df_n$alpha, df_n$beta, pdf_quasi, link = self$spec$transform$link))
                                  } else {
                                    # True mixture density
                                    pdf_Yt <- function(x) dmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
                                    # True Log-likelihood
                                    moments$loglik[i] <- log(dsolarGHI(df_n$GHI, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = self$spec$transform$link))
                                  }
                                }
                              } else {
                                if (quasi){
                                  # Standardize the time series
                                  moments$z <- (moments$Yt - moments$e_Yt) / moments$sd_Yt
                                  # Pseudo log-likelihood
                                  moments$loglik <- log(dnorm(moments$z) / moments$sd_Yt)
                                } else {
                                  # Standardize the time series into its components
                                  moments$z1 <- (moments$Yt - moments$M_Y1) / moments$S_Y1
                                  moments$z0 <- (moments$Yt - moments$M_Y0) / moments$S_Y0
                                  # True mixture log-likelihood
                                  moments$loglik <- log(dnorm(moments$z1) / moments$S_Y1 * moments$p1 + dnorm(moments$z0) / moments$S_Y0 * (1 - moments$p1))
                                  sum(moments$loglik)
                                  nrow(moments)
                                }
                              }
                              # **************************************************** #
                              return(moments$loglik)
                            },
                            #' @description
                            #' Print method for `solarModel` class.
                            print = function(){
                              # Complete data specifications
                              data <- self$spec$dates$data
                              # Train data specifications
                              train <- self$spec$dates$train
                              train$perc <- format(train$perc*100, digits = 4)
                              # Test data specifications
                              test <- self$spec$dates$test
                              test$perc <- format(test$perc*100, digits = 4)
                              cat(paste0("--------------------- ", "solarModel", " (", "\033[1;35m", self$place, "\033[0m", ") ", "--------------------- \n"))
                              cat(paste0("Model: ", "ARMA", "(", self$ARMA$order[1], ", ", self$ARMA$order[2], ")-GARCH", "(", self$GARCH$order[1], ", ", self$GARCH$order[2], ")\n"))
                              cat(paste0("Target: ", self$spec$target, " \n Coordinates: (Lat: ", self$spec$coords$lat, ", Lon: ", self$spec$coords$lon, ", Alt: ", self$spec$coords$alt, ") \n"))
                              cat(paste0(" Dates: ", data$from, " - ", data$to, "\n Observations: ", data$nobs, "\n"))
                              cat(paste0("---------------------------------------------------------------\n"))
                              cat(paste0("Train dates: ", train$from, " - ", train$to, " (", train$nobs, " points ~ ", train$perc, "%)", "\n"))
                              cat(paste0(" Test dates: ", test$from, " - ", test$to, " (", test$nobs, " points ~ ", test$perc, "%)", "\n"))
                              cat(paste0("---------------------------------------------------------------\n"))
                              cat(paste0("Log-Likelihood: ", format(self$loglik, digits = 8), "\n"))
                              cat(paste0("Empiric Mixture Parameters: ", ifelse(is.null(self$NM_model$use_empiric), "FALSE",
                                                                                self$NM_model$use_empiric), "\n"))
                              cat(paste0("Interpolated: ", private$interpolated, "\n"))
                              cat(paste0("Version: ", private$version, "\n"))
                            }
                          ),
                          # ====================================================================================================== #
                          #                                             Private slots
                          # ====================================================================================================== #
                          private = list(
                            version = "1.0.1",
                            ..data = NA,
                            ..seasonal_data = NA,
                            ..monthly_data = NA,
                            ..spec = NA,
                            ..transform = NA,
                            outliers = NA,
                            ..seasonal_model_Ct = NA,
                            ..seasonal_model_Yt = NA,
                            ..ARMA = list(),
                            ..seasonal_variance = NA,
                            ..GARCH = list(),
                            ..NM_model = list(),
                            ..loglik = NA,
                            ..hessian = NA,
                            ..jacobian = NA,
                            ..moments = list(conditional = NA, unconditional = NA),
                            interpolated = FALSE
                          ),
                          # ====================================================================================================== #
                          #                                             Active slots
                          # ====================================================================================================== #
                          active = list(
                            #' @field place Character, optional name of the location considered.
                            place = function(){
                              private$..spec$place
                            },
                            #' @field model_name Character, model's name.
                            model_name = function(){
                              paste0(self$spec$transform$link, "-ARMA(", self$ARMA$order[1], ", ", self$ARMA$order[2], ")",
                                     "GARCH(", self$GARCH$order[1], ", ", self$GARCH$order[2], ")")
                            },
                            #' @field data A data frame with the fitted data, and the seasonal and monthly parameters.
                            data = function(){
                              # Seasonal data
                              seasonal_data <- dplyr::select(self$seasonal_data, -n)
                              dplyr::left_join(private$..data, seasonal_data, by = c("Month", "Day"))
                            },
                            #' @field seasonal_data A data frame containing seasonal and monthly parameters.
                            seasonal_data = function(){
                              dplyr::left_join(private$..seasonal_data, self$monthly_data, by = "Month")
                            },
                            #' @field monthly_data A data frame that contains monthly parameters.
                            monthly_data = function(){
                              if (!purrr::is_empty(self$NM_model)) {
                                dplyr::left_join(private$..monthly_data, self$NM_model$coefficients, by = "Month")
                              } else {
                                private$..monthly_data
                              }
                            },
                            #' @field loglik The log-likelihood computed on train data.
                            loglik = function(){
                              private$..loglik
                            },
                            #' @field spec A list with the specification that govern the behavior of the model's fitting process.
                            spec = function(){
                              private$..spec
                            },
                            #' @field location A data frame with coordinates of the location considered.
                            location = function(){
                              dplyr::bind_cols(place = self$place, target = self$spec$target, dplyr::bind_rows(self$spec$coords))
                            },
                            #' @field transform A \code{\link{solarTransform}} object with the transformation functions applied to the data.
                            transform = function(){
                              private$..transform
                            },
                            #' @field seasonal_model_Ct The fitted model for clear sky radiation, used for predict the maximum radiation available.
                            seasonal_model_Ct = function(){
                              private$..seasonal_model_Ct
                            },
                            #' @field seasonal_model_Yt The fitted seasonal model for the target variable.
                            seasonal_model_Yt = function(){
                              private$..seasonal_model_Yt
                            },
                            #' @field ARMA The fitted ARMA model for the target variable.
                            ARMA = function(){
                              private$..ARMA
                            },
                            #' @field seasonal_variance The fitted model for seasonal variance.
                            seasonal_variance = function(){
                              private$..seasonal_variance
                            },
                            #' @field GARCH A model object representing the GARCH model fitted to the residuals.
                            GARCH = function(){
                              private$..GARCH
                            },
                            #' @field NM_model A model object representing the Gaussian Mixture model fitted to the standardized residuals.
                            NM_model = function(){
                              private$..NM_model
                            },
                            #' @field moments Get a list containing the conditional and unconditional moments.
                            moments = function(){
                              moments <- private$..moments
                              # Extra data to add
                              data_extra <- dplyr::select(self$data, date, GHI_bar, Ct, p1)
                              data_extra <- dplyr::mutate(data_extra, alpha = self$transform$alpha, beta = self$transform$beta)
                              # Conditional moments
                              if (length(moments$conditional) != 1) {
                                moments$conditional <- dplyr::left_join(moments$conditional, data_extra, by = c("date"))
                              }
                              # Unconditional moments
                              if (length(moments$unconditional) != 1) {
                                moments$unconditional <- dplyr::left_join(moments$unconditional, data_extra, by = c("date"))
                              }
                              return(moments)
                            },
                            #' @field coefficients Get the model parameters as a named list.
                            coefficients = function(){
                              # 1. Clear sky seasonal model
                              coefs_names <- c()
                              seasonal_model_Ct <- as.list(self$seasonal_model_Ct$coefficients)
                              # 2. Seasonal model Yt
                              seasonal_model_Yt <- as.list(self$seasonal_model_Yt$coefficients)
                              # 3. AR model
                              ARMA <- as.list(self$ARMA$coefficients)
                              # 4. Seasonal variance model
                              seasonal_variance <- as.list(self$seasonal_variance$coefficients)
                              # 5. GARCH variance model
                              GARCH <- as.list(self$GARCH$coefficients)
                              # 6. Gaussian mixture model
                              NM_mu_up <- self$NM_model$coefficients$mu1
                              names(NM_mu_up) <- paste0("mu_up_", 1:12)
                              NM_mu_dw <- self$NM_model$coefficients$mu2
                              names(NM_mu_dw) <- paste0("mu_dw_", 1:12)
                              NM_sd_up <- self$NM_model$coefficients$sd1
                              names(NM_sd_up) <- paste0("sd_up_", 1:12)
                              NM_sd_dw <- self$NM_model$coefficients$sd2
                              names(NM_sd_dw) <- paste0("sd_dw_", 1:12)
                              NM_p_up <- self$NM_model$coefficients$p1
                              names(NM_p_up) <- paste0("p_", 1:12)
                              # Output list of parameter for each model
                              params <- list(
                                location = self$location,
                                params = dplyr::bind_cols(alpha = self$transform$alpha, beta = self$transform$beta),
                                seasonal_model_Ct = dplyr::bind_cols(seasonal_model_Ct),
                                seasonal_model_Yt = dplyr::bind_cols(seasonal_model_Yt),
                                ARMA = dplyr::bind_cols(ARMA),
                                seasonal_variance = dplyr::bind_cols(seasonal_variance),
                                GARCH = dplyr::bind_cols(GARCH),
                                NM_mu_up = dplyr::bind_rows(NM_mu_up),
                                NM_mu_dw = dplyr::bind_rows(NM_mu_dw),
                                NM_sd_up = dplyr::bind_rows(NM_sd_up),
                                NM_sd_dw = dplyr::bind_rows(NM_sd_dw),
                                NM_p_up = dplyr::bind_rows(NM_p_up)
                              )
                              # Eventually add correlation (stochastic clearsky models)
                              check_rho_up <- suppressWarnings(!is.null(self$NM_model$rho_up[1]))
                              if (check_rho_up) {
                                NM_rho_up <- self$NM_model$coefficients$rho_up
                                names(NM_rho_up) <- paste0("rho_up_", 1:12)
                                params <- append(params, values = list(rho_up = dplyr::bind_rows(NM_rho_up)))
                              }
                              # Eventually add correlation (stochastic clearsky models)
                              check_rho_dw <- suppressWarnings(!is.null(self$NM_model$coefficients$rho_dw[1]))
                              if (check_rho_dw) {
                                NM_rho_dw <- self$NM_model$rho_dw
                                names(NM_rho_dw) <- paste0("rho_dw_", 1:12)
                                params <- append(params, values = list(rho_dw = dplyr::bind_rows(NM_rho_dw)))
                              }
                              return(params)
                            },
                            #' @field var_theta Variance-covariance matrix of the parameters with robust std. errors.
                            var_theta = function(){
                              H <- private$..hessian
                              J <- private$..jacobian
                              solve(H) %*% crossprod(J, J) %*% solve(H)
                            }
                          )
)

