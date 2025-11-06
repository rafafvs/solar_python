#' @title solarTransform
#' Solar functions
#' @description
#' Solar Model transformation functions
#'
#' @examples
#' st <- solarTransform$new()
#' st$GHI(0.4, 3)
#' st$GHI(st$iGHI(0.4, 3), 3)
#'
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarTransform <- R6::R6Class("solarTransform",
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
                                  private$..link <- match.arg(link, choices = c("invgumbel", "gumbel", "logis", "norm"))
                                },
                                #' @description
                                #' Map the risk driver X in solar radiation
                                #' @param x Numeric values in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{R_t(x) = C_t(1 - x)}
                                #' @return Numeric values in \eqn{C(t)(1-\alpha-\beta, 1-\alpha)}.
                                GHI = function(x, Ct){
                                  Ct*(1 - x)
                                },
                                #' @description
                                #' Map the solar radiation in the risk driver X
                                #' @param x Numeric values in \eqn{[C(t)(1-\alpha-\beta), C(t)(1-\alpha)]}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes the inverse of the `GHI`function.
                                #' \deqn{R_t^{-1}(x) = 1 - \frac{x}{C_t}}
                                #' @return Numeric values in \eqn{[\alpha,\alpha+\beta]}.
                                iGHI = function(x, Ct){
                                  1 - x/Ct
                                },
                                #' @description
                                #' Map the transformed variable Y in solar radiation
                                #' @param y Numeric values in \eqn{(-\infty, \infty)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{R_t(y) = C(t)(1 - \alpha-\beta \exp(-\exp(y)))}
                                #' @return Numeric values in \eqn{[C(t)(1-\alpha-\beta), C(t)(1-\alpha)]}.
                                GHI_y = function(y, Ct){
                                  Ct*(1 - self$iY(y))
                                },
                                #' @description
                                #' Map the risk driver X in the transformed variable Y
                                #' @param x numeric vector in \eqn{[\alpha, \alpha+\beta]}.
                                #' @details The function computes:
                                #' \deqn{Y(x) = \log(\log(\beta) - \log(x - \alpha))}
                                #' @return Numeric values in \eqn{[-\infty, \infty]}.
                                Y = function(x){
                                  x <- ifelse(x >= self$alpha + self$beta, self$alpha + self$beta - self$epsilon, x)
                                  x <- ifelse(x <= self$alpha, self$alpha + self$epsilon, x)
                                  private$..link
                                  self$Q_H(self$ieta(x))
                                },
                                #' @description
                                #' Map the transformed variable Y in the risk driver X.
                                #' @param y numeric vector in \eqn{[-\infty, \infty]}.
                                #' @details The function computes:
                                #' \deqn{Y^{-1}(y) = \alpha + \beta \exp(-\exp(y))}
                                #' @return Numeric values in \eqn{[\alpha, \alpha + \beta]}.
                                iY = function(y){
                                  self$eta(self$F_H(y))
                                },
                                #' @description
                                #' Map the risk driver X in the normalized variable Z.
                                #' Transformation function from X to Y
                                #' @param x numeric vector in \eqn{[\alpha, \alpha+\beta]}.
                                #' @details The function computes:
                                #' \deqn{\eta^{-1}(x) = \frac{x - \alpha}{\beta}}
                                #' @return Numeric values in \eqn{[0, 1]}.
                                ieta =function(x){
                                  (x - self$alpha) / self$beta
                                },
                                #' @description
                                #' Map the normalized variable Z in the risk driver X.
                                #' @param z numeric vector in \eqn{[0, 1]}.
                                #' @details The function computes:
                                #' \deqn{\eta(z) = \alpha + \beta \cdot z}
                                #' @return Numeric values in \eqn{[\alpha, \alpha + \beta]}.
                                eta = function(z){
                                  self$alpha + self$beta * z
                                },
                                #' @description
                                #' Fit the best parameters \eqn{\alpha} and \eqn{\beta} from a given time series
                                #' @param x time series of solar risk drivers in \eqn{(0, 1)}.
                                #' @param epsilon Numeric
                                #' @param min_pos Integer, position of the minimum. For example when `2` the minimum is the second lowest value.
                                #' @param max_pos Integer, position of the maximum. For example when `3` the maximum is the third greatest value.
                                #' @details Return a list that contains:
                                #' \describe{
                                #'  \item{alpha}{Numeric, \eqn{\alpha} transformation parameter.}
                                #'  \item{beta}{Numeric, \eqn{\beta} transformation parameter.}
                                #'  \item{epsilon}{Numeric, threshold used for fitting.}
                                #'  \item{Xt_min}{Numeric, minimum value of the time series.}
                                #'  \item{Xt_min}{Numeric, maximum value of the time series.}
                                #' }
                                #' @return A named list.
                                fit = function(x, epsilon = 0.01, min_pos = 1, max_pos = 1){
                                  # Index for minimum and maximum
                                  idx_xmin <- order(x, decreasing = FALSE)
                                  idx_xmax <- order(x, decreasing = TRUE)
                                  # Upper and lower bounds
                                  xmin <- x[idx_xmin]
                                  xmax <- x[idx_xmax]
                                  # Custom range
                                  range_Xt <- c(xmin[min_pos], xmax[max_pos])
                                  # Approximation parameter
                                  self$epsilon <- range_Xt[1]*epsilon
                                  # Transformation parameters
                                  alpha_ <- range_Xt[1] - self$epsilon
                                  beta_ <- (range_Xt[2] - range_Xt[1]) + 2 * self$epsilon
                                  # Store Transform parameters
                                  # Transform parameters
                                  list(alpha = alpha_, beta = beta_, epsilon = self$epsilon,
                                       Xt_min = range_Xt[1], Xt_max = range_Xt[2],
                                       idx_Xt_min = idx_xmin[min_pos], idx_Xt_max = idx_xmax[max_pos])
                                },
                                #' @description
                                #' Compute the bounds for the transformed variables.
                                #' @param target target variable. Available choices are:
                                #' \describe{
                                #'  \item{`"Xt"`}{Solar risk driver, the bounds returned are \eqn{[\alpha, \alpha + \beta]}.}
                                #'  \item{`"Kt"`}{Clearness index, the bounds returned are \eqn{[1-\alpha-\beta, 1-\alpha]}.}
                                #'  \item{`"Yt"`}{Solar transform, the bounds returned are \eqn{[-\infty, \infty]}.}
                                #'}
                                #' @return A numeric vector where the first element is the lower bound and the second the upper bound.
                                bounds = function(target = "Xt"){
                                  target = match.arg(target, choices = c("Xt", "Yt", "Kt"))
                                  lower <- c(Xt_min = self$alpha)
                                  upper <- c(Xt_max = self$alpha + self$beta)
                                  if (target == "Yt") {
                                    lower <-  c(Yt_min = -Inf)
                                    upper <-  c(Yt_max = Inf)
                                  } else if (target == "Kt") {
                                    lower <-  c(Kt_min = 1-self$alpha-self$beta)
                                    upper <-  c(Kt_max = 1-self$alpha)
                                  }
                                  return(c(lower, upper))
                                },
                                #' @description
                                #' Update the transformation parameters \eqn{\alpha} and \eqn{\beta}.
                                #' @param alpha Numeric, transformation parameter.
                                #' @param beta Numeric, transformation parameter.
                                #' @return Update the slots `$alpha` and `$beta`.
                                update = function(alpha, beta) {
                                  # Old parameters
                                  if (missing(alpha)) {
                                    alpha <- private$..alpha
                                  }
                                  if (missing(beta)) {
                                    beta <- private$..beta
                                  }
                                  # Control alpha > 0
                                  if (alpha < 0) {
                                    warning("Error: alpha is lower than 0")
                                    return(invisible(NULL))
                                  }
                                  # Control alpha + beta <= 1
                                  if (alpha + beta > 1) {
                                    warning("Error: alpha + beta is greater than 1")
                                    return(invisible(NULL))
                                  }
                                  # Update parameters
                                  private$..alpha <- alpha
                                  private$..beta <- beta
                                },
                                #' @description
                                #' Print method for the class `solarTransform`
                                print = function(){
                                  alpha_ <- format(self$alpha, digits = 3, scientific = FALSE)
                                  beta_ <- format(self$beta, digits = 4, scientific = FALSE)
                                  cat("------------------------ \033[1;35m solarTransform \033[0m ------------------------", "\n")
                                  cat(paste0("Link: \033[1;32m ", private$..link, "\033[0m \n"))
                                  cat(paste0("iY(y): \033[1;32m ",  alpha_, "\033[0m + \033[1;32m", beta_, "\033[0m exp(-exp(Y)) \n"))
                                  cat(paste0("Y(x): log(log(\033[1;32m",  beta_, "\033[0m) - log(x- \033[1;32m", alpha_, "\033[0m))"))
                                }
                              ),
                              private = list(
                                version = "1.0.0",
                                ..alpha = NA,
                                ..beta = NA,
                                ..link = NA
                              ),
                              active = list(
                                #' @field alpha Numeric, \eqn{\alpha} transformation parameter.
                                alpha = function(){
                                  private$..alpha
                                },
                                #' @field beta Numeric, \eqn{\beta} transformation parameter.
                                beta = function(){
                                  private$..beta
                                },
                                #' @field F_H Function, distribution function of the transform.
                                F_H = function(){
                                  if (private$..link == "invgumbel") {
                                    function(x) exp(-exp(x))
                                  } else if (private$..link == "gumbel") {
                                    function(x) exp(-exp(-x))
                                  } else if (private$..link == "logis") {
                                    function(x) 1/(1+exp(-x))
                                  } else if (private$..link == "norm") {
                                    function(x) pnorm(x)
                                  }
                                },
                                #' @field Q_H Function, quantile function of the transform.
                                Q_H = function(){
                                  if (private$..link == "invgumbel") {
                                    function(x) log(-log(x))
                                  } else if (private$..link == "gumbel") {
                                    function(x) -log(-log(x))
                                  } else if (private$..link == "logis") {
                                    function(x) log(x/(1-x))
                                  } else if (private$..link == "norm") {
                                    function(x) qnorm(x)
                                  }
                                }
                              )
)


