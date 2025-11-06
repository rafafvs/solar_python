#' Spatial Grid
#'
#' Create a grid from a range of latitudes and longitudes.
#'
#' @return a tibble with two columns `lat` and `lon`.
#'
#' @examples
#' # Initialize a spatial grid
#' sp <- spatialGrid$new()
#'
#' # Create an equally spaced grid
#' sp$make_grid(c(43.1, 44), c(9.2, 12.3), c(0.1, 0.1))
#' sp
#' # Check known point
#' sp$is_known_point(49.95, 12.15)
#' sp$is_known_point(43.9, 12)
#'
#' # Check if a point is inside the bounds
#' sp$is_inside_bounds(44.8, 10.9)
#'
#' # Extract a point
#' sp$neighborhoods(43.9, 12.1)
#'
#' # Extract its neighborhoods
#' sp$neighborhoods(43.95, 12.15)
#'
#' @rdname spatialGrid
#' @name spatialGrid
#' @export
spatialGrid <- R6::R6Class("spatialGrid",
                           public = list(
                             #' @field weights Weighting function for the distance.
                             weights = NA,
                             #' @description
                             #' Initialize a spatial grid
                             #' @param weights Weighting function for the distance.
                             initialize = function(weights = IDW(2)){
                               # Weighting function
                               self$weights <- weights
                             },
                             #' @description
                             #' Set a spatial grid
                             #' @param grid Tibble with column `id`, `lat` and `lon`.
                             set_grid = function(grid){
                               # Extract unique latitudes
                               lat <- as.numeric(unique(grid$lat))
                               lat <- lat[order(lat)]
                               # Extract unique longitudes
                               lon <- as.numeric(unique(grid$lon))
                               lon <- lon[order(lon)]
                               # Store latitudes bounds
                               private$lat_min <- min(lat)
                               private$lat_max <- max(lat)
                               private$lat_by <- abs(diff(lat)[1])
                               # Sequence of latitudes
                               private$lats = lat
                               # Longitudes bounds
                               private$lon_min <- min(lon)
                               private$lon_max <- max(lon)
                               private$lon_by <- abs(diff(lon)[1])
                               # Sequence of longitudes
                               private$lons = lon
                               # Assign a unique ID for each place
                               if(is.null(grid[["id"]])){
                                 cli::cli_alert_warning("The column `id` is missing in the grid. Assigned new ID.")
                                 labels = paste0("ID_", 1:nrow(grid))
                                 grid <- dplyr::bind_cols(id = labels, grid)
                               }
                               private$..grid <- grid
                             },
                             #' @description
                             #' Create a spatial grid
                             #' @param lat Numeric vector, from which is extracted the minimum and maximum for latitude.
                             #' @param lon Numeric vector, from which is extracted the minimum and maximum for longitude.
                             #' @param by Numeric vector, the first element is used to establish the distance between two latitudes in the grid.
                             #' The second element (if present) is used to establish the distance between two longitudes in the grid.
                             #' @param digits Integer scalar, number of digits for latitudes and longitudes.
                             make_grid = function(lat, lon, by, digits = 5){
                               # Latitudes bounds
                               private$lat_min <- min(lat)
                               private$lat_max <- max(lat)
                               private$lat_by <- by[1]
                               # Longitudes bounds
                               private$lon_min <- min(lon)
                               private$lon_max <- max(lon)
                               private$lon_by <- ifelse(length(by) == 1, by[1], by[2])

                               # Sequence of latitudes
                               private$lats = round(seq(private$lat_min, private$lat_max, by = private$lat_by), digits = digits)
                               # Sequence of longitudes
                               private$lons = round(seq(private$lon_min, private$lon_max, by = private$lon_by), digits = digits)
                               # Matrix of latitudes and longitudes
                               grid_lat_lon <- dplyr::as_tibble(expand.grid(private$lats, private$lons))
                               # Standard column names
                               colnames(grid_lat_lon) <- c("lat", "lon")
                               # Assign a unique ID for each place
                               labels = paste0("ID_", 1:nrow(grid_lat_lon))
                               grid_lat_lon <- dplyr::bind_cols(id = labels, grid_lat_lon)
                               # Store the spatial grid
                               private$..grid <- grid_lat_lon
                             },
                             #' @description
                             #' Check if a location is inside the bounds of the grid.
                             #' @param lat Numeric vector, reference latitudes.
                             #' @param lon Numeric vector, reference longitudes.
                             #' @return `TRUE` when the point is inside the limits and `FALSE` otherwise.
                             is_inside_bounds = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               condition <- rep(TRUE, nrow(coords))
                               for(i in 1:nrow(coords)) {
                                 condition[i] <- coords$lon[i] >= private$lon_min & coords$lon[i] <= private$lon_max
                                 condition[i] <- condition[i] & coords$lat[i] >= private$lat_min & coords$lat[i] <= private$lat_max
                               }
                               return(condition)
                             },
                             #' @description
                             #' Check if a location is a known point inside the grid.
                             #' @param lat Numeric vector, reference latitudes.
                             #' @param lon Numeric vector, reference longitudes.
                             #' @return `TRUE` when the location is known and `FALSE` otherwise.
                             is_known_point = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               condition <- rep(NA, nrow(coords))
                               for(i in 1:nrow(coords)) {
                                 condition_ <- dplyr::filter(self$grid, lat == coords$lat[i] & lon == coords$lon[i])
                                 condition[i] <- nrow(condition_) != 0
                               }
                               return(condition)
                             },
                             #' @description
                             #' Return the id and coordinates of a location inside the grid.
                             #' @param lat Numeric vector, reference latitudes.
                             #' @param lon Numeric vector, reference longitudes.
                             known_point = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               output <- dplyr::tibble()
                               for(i in 1:nrow(coords)) {
                                 if (self$is_known_point(coords$lat[i], coords$lon[i])) {
                                   df_known_point <- dplyr::filter(self$grid, lat == coords$lat[i] & lon == coords$lon[i])
                                   output <- dplyr::bind_rows(output, df_known_point)
                                 }
                               }
                               return(output)
                             },
                             #' @description
                             #' Find the n-closest neighborhoods of a location.
                             #' @param lat Numeric scalar, reference latitude.
                             #' @param lon Numeric scalar, reference longitude.
                             #' @param n number of neighborhoods
                             neighborhoods = function(lat, lon, n = 4){
                               if (self$is_known_point(lat, lon)) {
                                 nb <- self$known_point(lat, lon)
                                 # Add target latitude and longitude
                                 nb$lat_E <- lat
                                 nb$lon_E <- lon
                                 # No interpolation if it is a known point
                                 nb$wgt <- 1
                               } else {
                                 # Grid of locations
                                 nb <- self$grid
                                 # Add target latitude and longitude
                                 nb$lat_E <- lat
                                 nb$lon_E <- lon
                                 # Compute distances from target
                                 nb <- dplyr::mutate(nb, dist = havDistance(lat, lon, lat_E, lon_E))
                                 # Arrange by ascending distances
                                 nb <- head(dplyr::arrange(nb, dist), n = n)
                                 # Compute the weights
                                 nb$wgt <- self$weights(nb$dist)
                                 # Normalize the weights
                                 nb$wgt <- nb$wgt/sum(nb$wgt)
                               }
                               # Output
                               return(nb)
                             },
                             #' @description
                             #' Method print for a `spatialGrid` object.
                             print = function(){
                               cat(paste0("--------------------- ", "spatial Grid", "--------------------- \n"))
                               cat(paste0("Latitudes: ", "(", "\033[1;35m", round(private$lat_min, 3), "\033[0m", " - ", "\033[1;35m", round(private$lat_max, 3), "\033[0m", "; by ", round(private$lat_by, 3), ")", "\n"))
                               cat(paste0("Longitudes: ", "(", "\033[1;35m", round(private$lon_min, 3), "\033[0m", " - ", "\033[1;35m", round(private$lon_max, 3), "\033[0m", "; by ", round(private$lon_by, 3), ")", "\n"))
                               cat(paste0("Number of points: ", "\033[1;35m", nrow(self$grid), "\033[0m", "\n"))
                             }
                           ),
                           private = list(
                             ..grid = NA,
                             lats = NA,
                             lons = NA,
                             lat_min = NA,
                             lat_max = NA,
                             lon_min = NA,
                             lon_max = NA,
                             lon_by = NA,
                             lat_by = NA
                           ),
                           active = list(
                             #' @field grid A tibble with the spatial grid.
                             grid = function(){
                               private$..grid
                             }
                           )
)

#' Haversine distance
#'
#' Compute the Haversine distance between two points.
#'
#' @param lat_1 Numeric vector, latitudes of first location.
#' @param lon_1 Numeric vector, longitudes of first location.
#' @param lat_2 Numeric vector, latitudes of second location.
#' @param lon_2 Numeric vector, longitudes of second location.
#'
#' @examples
#' havDistance(43.3, 12.1, 43.4, 12.2)
#' havDistance(c(43.35, 43.35), c(12.15, 12.1), c(43.4, 44.5), c(12.2, 13.4))
#' @name havDistance
#' @rdname havDistance
#' @return Vector of distances in kilometers.
#' @export
havDistance <- function(lat_1, lon_1, lat_2, lon_2){
  r <- 6371 # Radius of the earth in kilometres
  dist <- c()
  # Solar functions
  sf <- seasonalSolarFunctions$new()
  # Compute distances
  for(i in 1:length(lat_1)){
    # Coordinates
    coord_A <- c(lat_1[i], lon_1[i])
    coord_B <- c(lat_2[i], lon_2[i])
    # Convert in radiant
    point_A <- sf$radiant(coord_A)
    point_B <- sf$radiant(coord_B)
    # Latitudes
    phi <- c(point_A[1], point_B[1])
    delta_phi <- phi[2] - phi[1]
    # Longitudes
    lam <- c(point_A[2], point_B[2])
    delta_lam <- lam[2] - lam[1]
    # Compute distance
    dist[i] <- 2*r*asin(sqrt(0.5*(1 - cos(delta_phi) + cos(phi[1])*cos(phi[2])*(1-cos(delta_lam)))))
  }
  return(dist)
}


#' Inverse Distance Weighting Functions
#'
#' Return a distance weighting function
#'
#' @param beta parameter used in exponential and power functions.
#' @param d0 parameter used only in exponential function.
#' @details When the parameter `d0` is not specified the function returned will be of power type otherwise of exponential type.
#'
#' @examples
#' # Power weighting
#' IDW_pow <- IDW(2)
#' IDW_pow(c(2, 3,10))
#' IDW_pow(c(2, 3,10), normalize = TRUE)
#' # Exponential weighting
#' IDW_exp <- IDW(2, d0 = 5)
#' IDW_exp(c(2, 3,10))
#' IDW_exp(c(2, 3,10), normalize = TRUE)
#' @export
IDW <- function(beta, d0){
  # Power function
  IDW.pow <- function(beta) {
    function(d, normalize = FALSE){
      w <- 1/(d^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Exponential function
  IDW.exp <- function(beta, d0) {
    function(d, normalize = FALSE){

      w <- exp(-(d/d0)^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Output function
  if (missing(d0)) {
    IDW.pow(beta)
  } else {
    IDW.exp(beta, d0)
  }
}
