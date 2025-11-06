#deltaTemp <- PVGIS_data$Palermo$DeltaT
#GHI <- PVGIS_data$Palermo$GHI
#date <- PVGIS_data$Palermo$Date
#lat <- PVGIS_data$Palermo$Lat[1]
bristowCampbel <- R6::R6Class("bristowCampbel",
                              public = list(
                                lat = NA_integer_,
                                initialize = function(deltaTemp, GHI, date, lat){
                                  # Complete dataset
                                  data = dplyr::tibble(date = date, Month = lubridate::month(date), GHI = GHI, deltaTemp = deltaTemp)
                                  # Add extraterrestrial radiation
                                  data$H0 = seasonalSolarFunctions$new()$H0(data$date, lat)
                                  # Compute the monthly mean for deltaTemp
                                  seasonal_data = dplyr::group_by(data, Month) %>%
                                    dplyr::summarise(deltaT_avg = mean(deltaTemp, na.rm = TRUE))
                                  # Compute the monthly mean for deltaTemp
                                  data = dplyr::left_join(data, seasonal_data, by = "Month") %>%
                                    dplyr::mutate(
                                      Kt = GHI/H0,
                                      B = 0.036*exp(-0.154 * deltaT_avg)) %>%
                                    dplyr::ungroup()

                                  private$..data <- data
                                  private$..seasonal_data <- seasonal_data
                                  self$lat <- lat
                                  private$..parameters <- list(A = max(data$Kt), C = 1)
                                },
                                B = function(nmonth){
                                  deltaT_avg <- private$..seasonal_data$deltaT_avg
                                  0.036*exp(-0.154 * deltaT_avg[nmonth])
                                },
                                fit = function(){

                                  # Loss function for C parameter
                                  loss_function <- function(par, data){

                                    A = self$parameters$A
                                    C = par
                                    Kt = self$data$Kt
                                    B = self$data$B
                                    deltaTemp = self$data$deltaTemp

                                    errors <- Kt - A*(1- exp(-B*(deltaTemp)^C))

                                    sum(errors^2, na.rm = TRUE)
                                  }
                                  # Optimal parameter C
                                  private$..parameters$C = optimize(loss_function, interval = c(0, 10))$minimum


                                },
                                predict = function(deltaTemp, date){
                                  if (missing(deltaTemp)) {
                                    deltaTemp <- self$data$deltaTemp
                                    date <- self$data$date
                                    H0 <- self$data$H0
                                  } else {
                                    # Add extraterrestrial radiation
                                    H0 = seasonalSolarFunctions$new()$H0(date, self$lat)
                                  }

                                  A = self$parameters$A
                                  C = self$parameters$C
                                  B = self$B(lubridate::month(date))

                                  H0*A*(1- exp(-B*(deltaTemp)^C))
                                }
                              ),
                              private = list(
                                ..data = NA,
                                ..seasonal_data = NA,
                                ..parameters = NA
                              ),
                              active = list(
                                data = function(){
                                  private$..data
                                },
                                parameters = function(){
                                  private$..parameters
                                }
                              ))


#bc <- bristowCampbel$new(deltaTemp, GHI, date, lat)
#bc$fit()
#bc$parameters
#errors <- abs(bc$data$GHI-bc$predict())/bc$data$GHI*100
#errors <- errors[!is.infinite(errors)]
#bc$predict(deltaTemp = 0.1, "2022-06-01")

#ggplot()+
#  geom_line(data = PVGIS_data$Palermo[1:365,], aes(Date, GHI), linewidth = 0.2)+
#  geom_line(data = PVGIS_data$Palermo[1:365,], aes(Date, bc$predict()[1:365]), color = "red", alpha = 0.5)
