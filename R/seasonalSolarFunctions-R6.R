#' Solar seasonal functions
#'
#' @examples
#' dates <- seq.Date(as.Date("2022-01-01"), as.Date("2022-12-31"), 1)
#' # Seasonal functions object
#' sf <- seasonalSolarFunctions$new()
#'
#' # Adjustment parameter
#' sf$B(number_of_day(dates))
#'
#' # Time adjustment in minutes
#' sf$E(dates)
#'
#' # Declination
#' sf$declination(dates)
#'
#' # Solar constant
#' sf$Gsc
#'
#' # Solar constant adjusted
#' sf$Gon(dates)
#'
#' # Extraterrestrial radiation
#' sf$Hon(dates, 43)
#'
#' # Number of hours of sun
#' sf$sun_hours(dates, 43)
#'
#' # Sunset hour angle
#' sf$sunset_hour_angle(dates, 43)
#'
#' sf$solar_time("2022-01-01 12:00", 11, 10)
#' sf$hour_angle("2022-01-01 14:00", 11, 15)
#' sf$incidence_angle("2022-06-01 21:00", 31, 12, lon_st = 15, beta = 0, gamma = 0)
#' sf$azimut_angle("2022-01-01 14:00", 30, 17, lon_st = 15)
#'
#' @name seasonalSolarFunctions
#' @rdname seasonalSolarFunctions
#' @note Version 1.0.0
#' @references Duffie, Solar Engineering of Thermal Processes Fourth Edition.
#' @export
seasonalSolarFunctions <- R6::R6Class("seasonalSolarFunctions",
                                public = list(
                                  #' @field legal_hour Logical, when `TRUE` the clock time will be corrected for the legal hour.
                                  legal_hour = TRUE,
                                  #' @description
                                  #' Initialize a `seasonalSolarFunctions` object
                                  #' @param method character, method type for computations. Can be `cooper` or `spencer`.
                                  #' @param legal_hour Logical, when `TRUE` the clock time will be corrected for the legal hour.
                                  initialize = function(method = "spencer", legal_hour = TRUE){
                                    # Method of computation
                                    private$method_ <- match.arg(method, choices = c("spencer", "cooper"))
                                    self$legal_hour <- legal_hour
                                  },
                                  #' @description
                                  #' Extract or update the method used for computations.
                                  #' @param x character, method type. Can be `cooper` or `spencer`.
                                  #' @return When `x` is missing it return a character containing the method that
                                  #' is actually used.
                                  update_method = function(x){
                                    if (missing(x)) {
                                      return(private$method_)
                                    } else {
                                      # Old method
                                      old_method <- private$method_
                                      # New method
                                      private$method_ <- match.arg(x, choices = c("spencer", "cooper"))
                                      return(private$method_)
                                    }
                                  },
                                  #' @description
                                  #' Seasonal adjustment parameter.
                                  #' @param n number of the day of the year
                                  #' @details The function implement Eq. 1.4.2 from Duffie (4th edition), i.e.
                                  #' \deqn{B(n) = \frac{2\pi}{365} n}
                                  B = function(n){
                                    (2*base::pi*n)/365
                                  },
                                  #' @description
                                  #' Convert angles in radiant into an angles in degrees.
                                  #' @param x numeric vector, angles in radiant.
                                  #' @details The function computes:
                                  #' \deqn{\frac{x 180}{\pi}}
                                  degree = function(x){
                                    x*180/base::pi
                                  },
                                  #' @description
                                  #' Convert angles in degrees into an angles in radiant
                                  #' @param x numeric vector, angles in degrees.
                                  #' @details The function computes:
                                  #' \deqn{\frac{x \pi}{180}}
                                  radiant = function(x){
                                    x*base::pi/180
                                  },
                                  #' @description
                                  #' Compute the time adjustment in minutes.
                                  #' @param n number of the day of the year
                                  #' @details The function implement Eq. 1.5.3 from Duffie (4th edition), i.e.
                                  #' \deqn{E = 229.2(0.000075 + 0.001868 \cos(B) - 0.032077\sin(B) - 0.014615\cos(2B) - 0.04089\sin(2B))}
                                  #' @return The time adjustment in minutes.
                                  E = function(n){
                                    n <- number_of_day(n)
                                    B <- self$B(n-1)
                                    # Time adjustment in minutes
                                    E <- 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B))
                                    # Time adjustment in seconds
                                    return(lubridate::dminutes(E))
                                  },
                                  #' @description
                                  #' Compute the angle in the degree given a certain altitude in meters.
                                  #' @param alt Numeric, altitude in meters.
                                  elevation = function(alt){
                                    Re <- 6371*1000
                                    alpha <- acos(Re/(Re+alt))
                                    self$degree(alpha)
                                  },
                                  #' @description
                                  #' Compute the solar time from a clock time.
                                  #' @param x datetime, clock hour.
                                  #' @param lon longitude of interest in degrees.
                                  #' @param lon_st longitude of the Local standard meridian in degrees.
                                  #' @param tz Character, reference time zone.
                                  #' @details The function implement Eq. 1.5.2 from Duffie (4th edition), i.e.
                                  #' \deqn{solartime = clocktime + 4 (lon-lon_{st}) + E(n)}
                                  #' @return A datetime object
                                  solar_time = function(x, lon, lon_st = 15, tz = "Europe/Rome"){
                                    # Convert in a datetime
                                    date_h <- as.POSIXlt(x, tz = tz)
                                    # Adjust the date for legal hour
                                    if (FALSE) {
                                      start_legal_hour <- as.Date(paste0(lubridate::year(date_h[1]), "-03-10"))
                                      end_legal_hour <- as.Date(paste0(lubridate::year(date_h[1]), "-10-30"))
                                      #   - from 27 Mar - 30 Oct (TRUE)
                                      #   - from 30 Oct - 27 Mar (FALSE)
                                      is_legal_hour <- (as.Date(date_h) >= start_legal_hour) & (as.Date(date_h) <= end_legal_hour)
                                      date_h <- dplyr::if_else(is_legal_hour, date_h - lubridate::dhours(1), date_h)
                                    }
                                    # Solar hour for the selected day-time
                                    date_h + lubridate::dminutes(4*(lon-lon_st)) + self$E(date_h)
                                  },
                                  #' @description
                                  #' Compute the solar hour for a specific clock time.
                                  #' @param LST datetime, true solar time.
                                  #' @return Hours
                                  solar_hour = function(LST){
                                    # Local solar time (sum Hours + minutes separately)
                                    lubridate::hour(LST) + lubridate::minute(LST)/60 + lubridate::second(LST)/3600
                                  },
                                  #' @description
                                  #' Compute the solar angle for a specific hour of the day.
                                  #' @param LST datetime, true solar time.
                                  #' @details The function implement Eq. 1.42 from Comini (2013), i.e.
                                  #' \deqn{\omega = 15 (solarhour - 12)}
                                  #' where the "solarhour" is expressed in hours.
                                  #' @return An angle in degrees
                                  hour_angle = function(LST){
                                    # Local solar time (sum Hours + minutes separately)
                                    solar_hour <- self$solar_hour(LST)
                                    # Solar angle in degrees
                                    15*(solar_hour - 12)
                                  },
                                  #' @description
                                  #' Compute the incidence angle
                                  #' @param LST datetime, true solar time.
                                  #' @param lat latitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param beta altitude
                                  #' @param gamma orientation
                                  #' @return An angle in degrees
                                  incidence_angle = function(LST, lat, alt = 0, beta = 0, gamma = 0){
                                    # Altitude of the surface in radiant
                                    beta = self$radiant(beta)
                                    # Orientation of the surface in radiant
                                    gamma = self$radiant(gamma)
                                    # Latitude in radiant
                                    phi = self$radiant(lat)
                                    # Solar hour angle
                                    omega = self$radiant(self$hour_angle(LST))
                                    # Declination in radiant
                                    delta = self$radiant(self$declination(LST))
                                    # Components
                                    T_ = sin(delta)*(sin(phi)*cos(beta) - cos(phi)*sin(beta)*cos(gamma))
                                    U_ = cos(delta)*(cos(phi)*cos(beta) + sin(phi)*sin(beta)*cos(gamma))
                                    V_ = cos(delta)*sin(beta)*sin(gamma)
                                    # Cosine of the angle of incidence
                                    cos_theta_z = T_ + U_*cos(omega) + V_*sin(omega)
                                    cos_theta_z[cos_theta_z < 0] <- 0
                                    # Check if is night
                                    #omega_s <- self$sunset_hour_angle(as.Date(LST), lat, alt)
                                    #is_night <- ifelse(omega < -omega_s | omega > omega_s, 1, 0)
                                    #cos_theta_z = cos_theta_z * (1 - is_night)
                                    # Angle of incidence
                                    #self$degree(acos(cos_theta_z))
                                    cos_theta_z
                                  },
                                  #' @description
                                  #' Compute the solar azimuth angle for a specific time of the day.
                                  #' @param LST datetime, true solar time.
                                  #' @param lat latitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param beta altitude
                                  #' @param gamma orientation
                                  #' @details The function implement Eq. 1.6.6 from Duffie (4th edition), i.e.
                                  #' \deqn{\gamma_s = sign(\omega) \left|\cos^{-1}\left( \frac{\cos \theta_z \sin \phi - \sin \delta}{\sin \theta_z \cos \phi} \right) \right|}
                                  #' @return The solar azimut angle in degrees
                                  azimut_angle = function(LST, lat, alt, beta = 0, gamma = 0){
                                    # Latitude in radiant
                                    phi = self$radiant(lat)
                                    # Declination in radiant
                                    delta = self$radiant(self$declination(LST))
                                    # Solar hour angle
                                    omega = self$radiant(self$hour_angle(LST))
                                    # The angle of incidence is the zenith angle of the sun
                                    theta_z = self$radiant(self$incidence_angle(LST, lat, alt, beta = 0, gamma = 0))
                                    # Azimut angle
                                    numer <- (cos(theta_z)*sin(phi) - sin(delta))
                                    denom <- sin(theta_z)*cos(phi)
                                    ratio <- pmin(pmax(numer / denom, -1), 1)
                                    gamma_s <- sign(omega)*abs(acos(angle))
                                    self$degree(gamma_s)
                                  },
                                  #' @description
                                  #' Compute the solar constant adjusted for the day of the year.
                                  #' @param n number of the day of the year.
                                  #' @param deriv Logical, when `TRUE` will return the first derivative with respect to time.
                                  #' @details When method is `cooper` the function implement Eq. 1.4.1a from Duffie (4th edition), i.e.
                                  #' \deqn{G_{o,n} = G_{sc} (1 + 0.033\cos(B))}
                                  #' otherwise when it is `spencer` it implement Eq. 1.4.1b from Duffie (4th edition):
                                  #' \deqn{G_{o,n} = G_{sc} (1.000110 + 0.034221\cos(B) + 0.001280\sin(B) + 0.000719\cos(2B) + 0.000077\sin(2B))}
                                  #' When `deriv = TRUE` it will be returned the derivatives with respect to time. When the method is `cooper`:
                                  #' \deqn{\frac{\partial G_{o,n}}{\partial n} = - G_{sc} \frac{2\pi}{365} 0.033 \sin(B))}
                                  #' Otherwise if it is `spencer`:
                                  #' \deqn{\frac{\partial G_{o,n}}{\partial n} = G_{sc} \frac{2\pi}{365} (-0.034221\sin(B) + 0.001280\cos(B) - 0.001438\sin(2B) + 0.000154\cos(2B))}
                                  #' @return The solar constant in \eqn{W/m^2} for the day n.
                                  Gon = function(n, deriv = FALSE){
                                    n <- number_of_day(n)
                                    Gsc <- private$..Gsc
                                    if (private$method_ == "cooper") {
                                      B <- self$B(n)
                                      if (!deriv) {
                                        Gn <- Gsc * (1 + 0.033 * cos(B))
                                      } else {
                                        Gn <- -Gsc * (2*base::pi/365) * 0.033 * sin(B)
                                      }
                                    } else if (private$method_ == "spencer") {
                                      B <- self$B(n-1)
                                      if (!deriv) {
                                        Gn <- Gsc * (1.000110 + 0.034221*cos(B) + 0.001280*sin(B) + 0.000719*cos(2*B) + 0.000077*sin(2*B))
                                      } else {
                                        Gn <- Gsc * (2*base::pi/365) * (-0.034221*sin(B) + 0.001280*cos(B) - 0.001438*sin(2*B) + 0.000154*cos(2*B))
                                      }
                                    }
                                    return(Gn)
                                  },
                                  #' @description
                                  #' Compute solar declination in degrees.
                                  #' @param n number of the day of the year
                                  #' @param deriv Logical, when `TRUE` will return the first derivative with respect to time.
                                  #' @details When method is `cooper` the function implement Eq. 1.6.1a from Duffie (4th edition), i.e.
                                  #' \deqn{\delta(n) = 23.45 \sin \left(\frac{2 \pi (284 + n)}{365}\right)}
                                  #' otherwise when it is `spencer` it implement Eq. 1.6.1b from Duffie (4th edition):
                                  #' \deqn{\delta(n) = \frac{180}{\pi}(0.006918 - 0.399912\cos(B) + 0.070257\sin(B) - 0.006758\cos(2B) + 0.000907\sin(2B) - 0.002697\cos(3B) + 0.00148\sin(3B))}
                                  #' When `deriv = TRUE` it will be returned the derivatives with respect to time. When the method is `cooper`:
                                  #' \deqn{\frac{\partial \delta}{\partial n}(n) = 23.45 \frac{2\pi}{365} \cos \left(\frac{2 \pi (284 + n)}{365}\right)}
                                  #' otherwise when the method is `spencer`:
                                  #' \deqn{\frac{\partial \delta}{\partial n}(n) = \frac{360}{365}(0.399912\sin(B) + 0.070257\cos(B) + 0.013516\sin(2B) + 0.001814\cos(2B) + 0.008091\sin(3B) + 0.00444\cos(3B))}
                                  #' @return The solar declination in degrees.
                                  declination = function(n, deriv = FALSE){
                                    n <- number_of_day(n)
                                    if (private$method_ == "cooper") {
                                      if (!deriv) {
                                        declination <- 23.45*sin(2*base::pi*(284 + n)/365)
                                      } else {
                                        declination <- 23.45*cos(2*base::pi*(284+n)/365) * (2*base::pi/365)
                                      }
                                    } else if (private$method_ == "spencer") {
                                      B <- self$B(n-1)
                                      if (!deriv) {
                                        declination <- (180/base::pi)*(0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - 0.006758*cos(2*B) + 0.000907* sin(2*B) - 0.002697*cos(3*B) + 0.00148*sin(3*B))
                                      } else {
                                        declination <- (360/365) * (0.399912*sin(B) + 0.070257*cos(B) + 0.013516*sin(2*B) + 0.001814*cos(2*B) + 0.008091*sin(3*B) + 0.00444*cos(3*B))
                                      }
                                    }
                                    return(declination)
                                  },
                                  #' @description
                                  #' Compute the solar extraterrestrial radiation
                                  #' @param n number of the day of the year
                                  #' @param lat latitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param deriv Logical, when `TRUE` will return the first derivative with respect to time.
                                  #' @details The function implement Eq. 1.10.3 from Duffie (4th edition):
                                  #' \deqn{H_{on} = G_{on} \frac{24 \times 3600}{\pi} (\cos(lat) \cos(\delta) \sin(\omega_s) + \frac{\pi}{180}\sin(lat) \sin(\delta))}
                                  #' @return Extraterrestrial radiation on an horizontal surface in kilowatt hour for meters squared for day.
                                  Hon = function(n, lat, alt, deriv = FALSE){
                                    # Latitude in radiant
                                    phi <- self$radiant(lat)
                                    # Declination in radiant
                                    delta <- self$radiant(self$declination(n))
                                    # Sunset hour angle in degrees
                                    omega_s <- self$radiant(self$sunset_hour_angle(n, lat, alt))
                                    # Scale factor
                                    A <- (24*3600)/(base::pi*3600000)
                                    # Extraterrestrial radiation in daily joules/m^2 per day
                                    Gon <- self$Gon(n)
                                    # Intern angle
                                    B_n <- (cos(phi)*cos(delta)*sin(omega_s) + omega_s*sin(phi)*sin(delta))
                                    if (!deriv) {
                                      Hon <- A * Gon * B_n
                                    } else {
                                      # Derivative of the declination in radiants
                                      d_delta <- self$radiant(self$declination(n, deriv = TRUE))
                                      # Derivative of sunset hour angle in radiants
                                      d_omega_s <- self$radiant(self$sunset_hour_angle(n, lat, alt, deriv = TRUE))
                                      # Derivative of B_n with respect to declination
                                      d_B_n_d_delta <- -sin(delta) * cos(phi) * sin(omega_s) + omega_s * cos(delta) * sin(phi)
                                      # Derivative of B_n with respect to sunset hour angle
                                      d_B_n_d_omega <- cos(delta) * cos(phi) * cos(omega_s) + sin(delta) * sin(phi)
                                      # Total derivative of B_n wrt to time
                                      d_B_n_d_n <- d_B_n_d_delta * d_delta + d_B_n_d_omega * d_omega_s
                                      # Derivative of solar constant wrt to time
                                      d_Gon_dn <- self$Gon(n, deriv = TRUE)
                                      # Total derivative (scaled)
                                      Hon <- A * (B_n * d_Gon_dn + Gon * d_B_n_d_n)
                                    }
                                    # Convert in (kilowatt hour)/m^2 per day
                                    return(Hon)
                                  },
                                  #' @description
                                  #' Compute solar angle at sunset in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param deriv Logical, when `TRUE` will return the first derivative with respect to time.
                                  #' @details The function implement Eq. 1.6.10 from Duffie (4th edition), i.e.
                                  #' \deqn{\omega_s = \cos^{-1}(-\tan(\delta(n))\tan(\phi))}
                                  #' When altitude is not missing it will implement a generalized version with altitude, i.e.
                                  #' \deqn{\omega_s = \cos^{-1}\left(\frac{\sin H - \sin\delta \sin\phi}{\cos\phi \cos \delta}\right)}
                                  #' @return The sunset hour angle in degrees.
                                  sunset_hour_angle = function(n, lat, alt, deriv = FALSE){
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    # Declination from degrees to radiant
                                    declination <- self$radiant(self$declination(n))
                                    # Version without altitude
                                    if (missing(alt)) {
                                      # Compute argument of acos
                                      arg <- -tan(declination) * tan(phi)
                                      # Return the angle
                                      if (!deriv) {
                                        # Safety clipping
                                        arg <- pmin(pmax(arg, -1), 1)
                                        # Sunset hour angle in degrees
                                        omega_s <- self$degree(acos(arg))
                                      } else {
                                        # Derivative of the declination
                                        d_declination <- self$declination(n, deriv = TRUE)
                                        # Derivative of the argument
                                        d_arg <- - tan(phi) / cos(declination)^2
                                        # Final derivative
                                        omega_s <- -(1 / sqrt(1 - arg^2)) * d_arg * d_declination
                                      }
                                      return(omega_s)
                                    }
                                    # Altitude from meters to radiant
                                    alt <- self$radiant(self$elevation(alt))
                                    # Compute argument of acos
                                    arg <- (sin(alt) - sin(declination)*sin(phi))/(cos(phi)*cos(declination))
                                    # Return the angle
                                    if (!deriv) {
                                      # Safety clipping
                                      arg <- pmin(pmax(arg, -1), 1)
                                      # Sunset hour angle in degrees
                                      omega_s <- self$degree(acos(arg))
                                    } else {
                                      # Return the first derivative
                                      d_declination <- self$declination(n, deriv = TRUE)
                                      # Derivative of the argument
                                      d_arg <- (sin(alt) * sin(declination) - sin(phi)) / (cos(phi) * cos(declination)^2)
                                      # Final derivative
                                      omega_s <- -(1 / sqrt(1 - arg^2)) * d_arg * d_declination
                                    }
                                    return(omega_s)
                                  },
                                  #' @description
                                  #' Compute number of sun hours for a day n.
                                  #' @param n number of the day of the year.
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @details The function implement Eq. 1.6.11 from Duffie (4th edition), i.e.
                                  #' \deqn{\frac{2}{15} \omega_s}
                                  sun_hours = function(n, lat, alt){
                                    # Number of sun hours
                                    sun_hours <- self$sunset_hour_angle(n, lat, alt)*(2/15)
                                    return(lubridate::dhours(sun_hours))
                                  },
                                  #' @description
                                  #' Compute solar altitude in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @details The function computes
                                  #' \deqn{\sin^{-1}(-\sin(\delta(n))\sin(\phi) + \cos(\delta(n))\cos(\phi))}
                                  solar_altitude = function(n, lat){
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    declination <- self$radiant(self$declination(n))
                                    # Safety clipping
                                    arg <- sin(declination)*sin(phi) + cos(declination)*cos(phi)
                                    arg <- pmin(pmax(arg, -1), 1)
                                    self$degree(asin(arg))
                                  },
                                  #' @description
                                  #' Compute the solar angle for a latitude in different dates.
                                  #' @param x datetime, clock hour.
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @param lon Numeric, longitude of interest in degrees.
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param lon_st Numeric, longitude of the Local standard meridian in degrees
                                  #' @param beta Numeric angle, inclination of the solar panel.
                                  #' @param gamma Numeric, angle orientation of the panel.
                                  #' @param by Character, time step. Default is `1 min`.
                                  #' @param tz Character, reference time zone.
                                  solar_angles = function(x, lat, lon, alt, lon_st = 15, beta = 0, gamma = 0, by = "1 min", tz = "Europe/Rome"){
                                    day_date <- as.Date(x[1])
                                    start_date <- as.POSIXct(paste0(day_date, " 00:00:00"), tz = tz)
                                    end_date <- as.POSIXct(paste0(day_date+1, " 00:00:00"), tz = tz)
                                    day_date_seq <- seq.POSIXt(start_date, end_date, by = by)
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    # Solar declination
                                    declination <- self$declination(day_date)
                                    # Solar angle at sunset
                                    omega_max <-  self$sunset_hour_angle(day_date, lat, alt)
                                    # Solar angle at sunrise
                                    omega_min <- -omega_max
                                    # Number of sun hours
                                    sun_hours <- self$sun_hours(day_date, lat, alt)
                                    # Time adjustment in seconds
                                    E <- self$E(day_date)
                                    # Solar time
                                    LST <- self$solar_time(day_date_seq, lon, lon_st, tz = tz)
                                    # Solar angle
                                    omega <- self$hour_angle(LST)
                                    # Incidence angle
                                    theta <- self$incidence_angle(LST, lat, alt, beta, gamma)
                                    dplyr::tibble(date = day_date,
                                                  clocktime = day_date_seq,
                                                  solartime = solartime,
                                                  lat = lat,
                                                  lon = lon,
                                                  omega = omega,
                                                  declination = declination,
                                                  omega_min = omega_min,
                                                  omega_max = omega_max,
                                                  sun_hours = sun_hours,
                                                  theta = theta,
                                                  E = E)
                                  },
                                  #' @description
                                  #' Hottel clearsky
                                  #' @param cosZ solar incidence angle
                                  #' @param G0 solar constant
                                  #' @param alt Numeric, altitude in meters.
                                  #' @param clime clime correction
                                  clearsky = function(cosZ = NULL, G0 = NULL, alt, clime = "No Correction"){
                                    # correction for different climes
                                    clime <- match.arg(clime[1], choices = c("No Correction", "Summer", "Winter", "Subartic Summer", "Tropical"))

                                    # altitude must be converted from metre to km
                                    altitude <- ifelse(missing(alt), 0, alt)
                                    if (altitude > 2.5) {
                                      a0_star <- 0.6*(1-exp(-0.214*(altitude - 1.12)))
                                    } else {
                                      a0_star <- 0.4237 - 0.00821*(6.0 - altitude)^2
                                    }
                                    a1_star <- 0.5055 - 0.00595*(6.5 - altitude)^2
                                    a2_star <- 0.2711 - 0.01858*(2.5 - altitude)^2

                                    # Correction for Clime types
                                    a <- c(a0_star, a1_star, a2_star)
                                    if (clime == "Summer") {
                                      correction_factor <- c(0.97, 0.99, 1.02)
                                      a <- a*correction_factor
                                    } else if (clime == "Winter") {
                                      correction_factor <- c(1.03, 1.01, 1.00)
                                      a <- a*correction_factor
                                    } else if (clime == "Subartic Summer"){
                                      correction_factor <- c(0.99, 0.99, 1.01)
                                      a <- a*correction_factor
                                    } else if (clime == "Tropical") {
                                      correction_factor <- c(0.95, 0.98, 1.02)
                                      a <- a*correction_factor
                                    }
                                    a0 <- a[1]
                                    a1 <- a[2]
                                    a2 <- a[3]
                                    output <- dplyr::tibble(tau_beam = a0 + a1*exp(-a2/cosZ), tau_diffuse = 0.271 - 0.294*tau_beam)
                                    output <- dplyr::mutate(output,
                                                            tau_beam = ifelse(tau_beam > 1 | tau_beam  < 0, 0, tau_beam),
                                                            tau_diffuse = ifelse(tau_diffuse > 1 | tau_diffuse  < 0, 0, tau_diffuse))
                                    skymax <- G0*output$tau_beam + G0*output$tau_diffuse
                                    return(skymax)
                                  }
                                ),
                                private = list(
                                  version = "1.0.0",
                                  ..Gsc = 1367,
                                  method_ = ""
                                ),
                                active = list(
                                  #' @field Gsc solar constant, i,e, `1367`.
                                  Gsc = function(){
                                    private$..Gsc
                                  }
                                ))
