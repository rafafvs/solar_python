#' Solar Model transformation functions
#'
#' @examples
#' st <- solarTransform$new()
#'
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarTransform <- R6::R6Class("solarTransform",
                              inherit = boundTransform,
                              public = list(
                                #' @field epsilon Numeric, \eqn{\epsilon} transformation parameter.
                                epsilon = 0,
                                #' @description
                                #' Initialize a `solarTransform` object.
                                #' @param alpha Numeric, \eqn{\alpha} transformation parameter.
                                #' @param beta Numeric, \eqn{\beta} transformation parameter.
                                #' @param link Character, link function.
                                initialize = function(alpha = 0, beta = 1, link = "invgumbel"){
                                  # Control parameters
                                  if (alpha < 0) {
                                    stop("Alpha cannot be lower than zero.")
                                  } else if (alpha + beta > 1) {
                                    stop("`alpha + beta` cannot be greater than one.")
                                  }
                                  # Update the parameters
                                  private$..alpha <- alpha
                                  private$..beta <- beta
                                  # Set the transform
                                  self$set_transform(link)
                                },
                                #' @description
                                #' Map the solar radiation \eqn{R_t} in the risk driver \eqn{X_t}.
                                #' @param Rt Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\text{X}(R_t) = 1 - R_t/C_t}
                                #' @return Numeric, risk driver \eqn{X_t \in (\alpha, \alpha+\beta)}.
                                X = function(Rt, Ct){
                                  1 - Rt / Ct
                                },
                                #' @description
                                #' Map the risk driver \eqn{X_t} in solar radiation \eqn{R_t}
                                #' @param Xt Numeric, risk driver in \eqn{ X_t \in (\alpha, \alpha+\beta)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\text{X}^{-1}(X_t) = C_t(1 - X_t)}
                                #' @return Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                iX = function(Xt, Ct){
                                  Ct * (1 - Xt)
                                },
                                #' @description
                                #' Map the solar radiation \eqn{R_t} in the normalized variable \eqn{X_t^{\prime}}.
                                #' @param Rt Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\eta(R_t) = \frac{1}{\beta}(1 - \alpha - R_t/C_t)}
                                #' @return Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                eta = function(Rt, Ct){
                                  self$X_prime(self$X(Rt, Ct))
                                },
                                #' @description
                                #' Map the normalized variable \eqn{X_t^{\prime}} to the solar radiatio \eqn{R_t}.
                                #' @param Xt_prime Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\eta^{-1}(X_t^{\prime}) = C_t(1 - \alpha - \beta \cdot X_t^{\prime})}
                                #' @return Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                ieta = function(Xt_prime, Ct){
                                  self$iX(self$iX_prime(Xt_prime), Ct)
                                },
                                #' @description
                                #' Convert solar radiation \eqn{R_t} into the transformed variable \eqn{Y_t}.
                                #' @param Rt Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\text{RY}(R_t) = g\left(\frac{1}{\beta}(1 - \alpha- R_t/C_t)\right)}
                                #' @return Transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                RY = function(Rt, Ct){
                                  # Cloudiness index
                                  Xt <- self$X(Rt, Ct)
                                  # Normalized risk driver
                                  Xt_prime <- self$X_prime(Xt)
                                  # Transformed variable
                                  self$Y(Xt_prime)
                                },
                                #' @description
                                #' Convert the transformed variable \eqn{Y_t} into solar radiation \eqn{R_t}.
                                #' @param Yt Numeric, transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{\text{iRY}(Y_t) = C_t(1 - \alpha-\beta g^{-1}(Y_t))}
                                #' @return Numeric, solar radiation \eqn{R_t \in [C_{t}(1-\alpha-\beta), C_{t}(1-\alpha)]}.
                                iRY = function(Yt, Ct){
                                  # Normalized risk driver
                                  Xt_prime <- self$iY(Yt)
                                  # Cloudiness index
                                  Xt <- self$iX_prime(Xt_prime)
                                  # Solar radiation
                                  self$iX(Xt, Ct)
                                }
                              ),
                              private = list(
                                version = "1.0.1"
                              ),
                              active = list()

)
