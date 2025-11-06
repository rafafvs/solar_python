#' `spatialKringing` object
#'
#' @rdname spatialKringing
#' @name spatialKringing
#' @export
spatialKringing <- R6::R6Class("spatialKringing",
                               public = list(
                                 #' @field quiet Logical
                                 quiet = FALSE,
                                 #' @description
                                 #' Initialize a `spatialKringing` object
                                 #' @param data dataset with spatial parameters and lon, lat.
                                 #' @param params_names Names of the parameters to fit.
                                 #' @param vg_models an optional list of kernelRegression models already fitted.
                                 #' @param grid description
                                 #' @param sample List of parameter used as sample.
                                 #' @param quiet description
                                 initialize = function(data, params_names, grid, vg_models, quiet = FALSE){
                                   # Dataset containing all the parameters
                                   private$..data <- data
                                   # Set quiet parameter
                                   self$quiet <- quiet
                                   # Names of all the parameters
                                   private$params_names <- params_names
                                   if (missing(vg_models)) {
                                     # Initialize a list for all the models
                                     private$..models <- as.list(rep(NA, length(private$params_names)))
                                     # Name the list with the parameter's names
                                     names(private$..models) <- private$params_names
                                   } else {
                                     if (any(names(vg_models) != private$params_names)){
                                       if (!self$quiet) warning("The parameters in `models` do not match the parameters names!")
                                     } else {
                                       # Update the models with the given list
                                       private$..models <- vg_models
                                     }
                                   }
                                   # Extract a row to fill with the fitted parameters
                                   sample <- data[1,]
                                   sample$place <- ""
                                   sample$lat <- NA_integer_
                                   sample$lon <- NA_integer_
                                   # Initialize the sample dataset to store the predictions
                                   private$sample <- sample
                                   # Store the spatial grid
                                   private$grid <- grid$clone(TRUE)
                                 },
                                 #' @description
                                 #' Fit a `kernelRegression` object for a parameter or a group of parameters.
                                 #' @param params list of parameters names to fit. When missing all the parameters will be fitted.
                                 fit = function(params){
                                   if (missing(params)) {
                                     params <- private$params_names
                                   } else {
                                     params <- match.arg(params, choices = private$params_names, several.ok = TRUE)
                                   }
                                   # number of parameters to fit
                                   n.params <- length(params)
                                   for(i in 1:n.params){
                                     param_name <- params[i]
                                     if (!self$quiet) message("Fitting param:", param_name, " ", i, "/", n.params, " parameters...")
                                     data_params_sp <- self$data[, c("lat", "lon", param_name)]
                                     colnames(data_params_sp) <- c("lat", "lon", "par")
                                     sp::coordinates(data_params_sp) <- ~lon+lat
                                     vg_emp <- gstat::variogram(par ~ 1, data_params_sp)  # residuals have no trend
                                     vg_model <- gstat::fit.variogram(vg_emp, model = vgm(model = "Sph"))
                                     # Store model
                                     private$..models[[param_name]] <- vg_model
                                   }
                                 },
                                 #' @description
                                 #' Predict all the parameters for a specified location.
                                 #' @param lat Numeric vector, latitudes in degrees.
                                 #' @param lon Numeric vector, longitudes in degrees.
                                 #' @param n Integer, number of neighborhoods to consider for interpolation.
                                 predict = function(lat, lon, n = 4){
                                   predictions <- list()
                                   coords <- dplyr::tibble(lat = lat, lon = lon)
                                   for(k in 1:nrow(coords)){
                                     # Extract reference latitude and longitude
                                     lat_k <- coords$lat[k]
                                     lon_k <- coords$lon[k]
                                     # Define prediction location (s_0)
                                     location_k <- data.frame(lon = lon_k, lat = lat_k)
                                     sp::coordinates(location_k) <- ~lon+lat

                                     if(!self$quiet) message("Fitting parameters for latitude: ", lat_k, " and longitude ", lon_k)
                                     # Detect closest points
                                     ngb <- private$grid$neighborhoods(lat_k, lon_k)
                                     # Extract a row to fill with the fitted parameters
                                     data_params_fit <- dplyr::mutate(private$sample, place = "", lat = lat_k, lon = lon_k)
                                     # Kringing variances
                                     data_params_fit_var <- data_params_fit
                                     # If the points is known skip interpolation
                                     if (nrow(ngb) == 1){
                                       # Predicted parameter
                                       data_params_fit <- dplyr::filter(self$data, place == ngb$id)
                                       data_params_fit_var <- data_params_fit
                                       data_params_fit_var[,-c(1:3)] <- 0
                                       # Store the predictions
                                       predictions[[k]] <- list()
                                       predictions[[k]]$mean <- data_params_fit
                                       predictions[[k]]$variance <- data_params_fit_var
                                       next
                                     }
                                     # Extract the known parameters
                                     data_params_ngb <- dplyr::filter(self$data, place %in% ngb$id)
                                     # Add coordinates
                                     sp::coordinates(data_params_ngb) <- ~lon+lat

                                     coords_ngb <- sp::coordinates(data_params_ngb)
                                     coords_k <- sp::coordinates(location_k)
                                     for(param_name in private$params_names){
                                       if(!self$quiet) message("Fitting parameter ", param_name, "...\r", appendLF = FALSE)
                                       vg_model <- private$..models[[param_name]]
                                       # Fill kriging matrix K
                                       K <- matrix(0, n+1, n+1)
                                       d <- numeric(n+1)
                                       for(i in 1:n){
                                         for(j in 1:n){
                                           h <- sp::spDistsN1(coords_ngb, coords_ngb[i,], longlat = FALSE)[j]
                                           K[i, j] <- gstat::variogramLine(vg_model, dist_vector = h)$gamma
                                         }
                                         K[i, n+1] <- 1
                                         K[n+1, i] <- 1
                                         # Right-hand side
                                         h0 <- sp::spDistsN1(coords_ngb, coords_k, longlat = FALSE)[i]
                                         d[i] <- gstat::variogramLine(vg_model, dist_vector = h0)$gamma
                                       }
                                       d[n+1] <- 1
                                       # Solve for weights and Lagrange multiplier
                                       sol <- solve(K, d)
                                       ngb$kringing <- sol[1:n]
                                       # Predicted parameter
                                       data_params_fit[[param_name]] <- sum(data_params_ngb@data[[param_name]] * ngb$kringing)
                                       data_params_fit_var[[param_name]] <- (ngb$kringing %*% d[1:n] + d[n])[1]
                                     }
                                     # Store the predictions
                                     predictions[[k]] <- list()
                                     predictions[[k]]$mean <- data_params_fit
                                     predictions[[k]]$variance <- data_params_fit_var
                                   }
                                   e_params <- purrr::map_df(predictions, ~.x$mean)
                                   v_params <- purrr::map_df(predictions, ~.x$variance)
                                   predictions <- list(mean = e_params, variance = v_params)
                                   return(predictions)
                                 }
                               ),
                               private = list(
                                 ..data = NA,
                                 params_names = NA,
                                 ..models = NA,
                                 grid = NA,
                                 sample = NA
                               ),
                               active = list(
                                 #' @field models list of `kernelRegression` objects
                                 models = function(){
                                   private$..models
                                 },
                                 #' @field data dataset with the parameters used for fitting
                                 data = function(){
                                   private$..data
                                 }
                               ))
