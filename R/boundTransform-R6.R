#' Bounded transformation functions
#'
#' @examples
#' st <- boundTransform$new()
#'
#' @keywords utils
#' @note Version 1.0.0.
#' @export
boundTransform <- R6::R6Class("boundTransform",
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
                                #' Set the transform function \eqn{g}, \eqn{g^{-1}}, \eqn{g^{\prime}}.
                                #' @param link Character, link function. Valid links are `"invgumbel"`, `"gumbel"`, `"logis"`, `"norm"`.
                                set_transform = function(link){
                                  private$..link <- match.arg(link, choices = c("invgumbel", "gumbel", "logis", "norm"))
                                  # Store the functions: g, g^{-1}, g_prime
                                  if (private$..link == "invgumbel") {
                                    private$..g <- function(x) log(-log(x))
                                    private$..ig <- function(x) exp(-exp(x))
                                    private$..g_prime <- function(x) 1/(x * log(x))
                                    private$..monotonicity <- "decreasing"
                                  } else if (private$..link == "gumbel") {
                                    private$..g <- function(x) -log(-log(x))
                                    private$..ig <- function(x) exp(-exp(-x))
                                    private$..g_prime <- function(x) -1/(x * log(x))
                                    private$..monotonicity <- "increasing"
                                  } else if (private$..link == "logis") {
                                    private$..g <- function(x) log(x/(1-x))
                                    private$..ig <- function(x) 1/(1+exp(-x))
                                    private$..g_prime <- function(x) 1 / (x * (1 - x))
                                    private$..monotonicity <- "increasing"
                                  } else if (private$..link == "norm") {
                                    private$..g <- function(x) qnorm(x)
                                    private$..ig <- function(x) pnorm(x)
                                    private$..g_prime <- function(x) 1 / dnorm(qnorm(x))
                                    private$..monotonicity <- "increasing"
                                  }
                                },
                                #' @description
                                #' Map the risk-driver \eqn{X_t} into the normalized risk-driver \eqn{X_t^{\prime}}.
                                #' @param Xt Numeric, risk driver in \eqn{ X_t \in (\alpha, \alpha+\beta)}.
                                #' @details The function computes:
                                #' \deqn{\text{X}^{\prime}(X_t) = \frac{X_t - \alpha}{\beta}}
                                #' @return Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                X_prime =function(Xt){
                                  Xt <- ifelse(Xt >= self$alpha + self$beta, self$alpha + self$beta - self$epsilon, Xt)
                                  Xt <- ifelse(Xt <= self$alpha, self$alpha + self$epsilon, Xt)
                                  (Xt - self$alpha) / self$beta
                                },
                                #' @description
                                #' Map the normalized variable \eqn{X_t^{\prime}} to the risk driver \eqn{X_t}.
                                #' @param Xt_prime Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                #' @details The function computes:
                                #' \deqn{\text{iX}^{\prime}(X_t^{\prime}) = \alpha + \beta \cdot X_t^{\prime}}
                                #' @return Numeric, risk driver in \eqn{ X_t \in (\alpha, \alpha+\beta)}.
                                iX_prime = function(Xt_prime){
                                  self$alpha + self$beta * Xt_prime
                                },
                                #' @description
                                #' Map the normalized risk driver \eqn{X_t^{\prime}} in the transformed variable \eqn{Y_t}
                                #' @param Xt_prime Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                #' @details The function computes:
                                #' \deqn{\text{Y}(X_t^{\prime}) = g(X_t^{\prime})}
                                #' @return Numeric, transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                Y = function(Xt_prime){
                                  self$g(Xt_prime)
                                },
                                #' @description
                                #' Map the transformed variable \eqn{Y_t} into the normalized risk driver \eqn{X_t^{\prime}}
                                #' @param Yt Numeric, transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                #' @details The function computes:
                                #' \deqn{\text{iY}(Y_t) = g^{-1}(Y_t)}
                                #' @return Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                iY = function(Yt){
                                  self$ig(Yt)
                                },
                                #' @description
                                #' Link function to map \eqn{X_t^{\prime}} to \eqn{Y_t}.
                                #' @param X_prime Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                #' @return Numeric, transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                g = function(X_prime){
                                  private$..g(X_prime)
                                },
                                #' @description
                                #' Inverse of the function to map \eqn{Y_t} to \eqn{X^{\prime}}.
                                #' @param Yt Numeric, transformed variable \eqn{Y_t \in (-\infty, \infty)}.
                                #' @return Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                ig = function(Yt){
                                  private$..ig(Yt)
                                },
                                #' @description
                                #' First derivative of the function \eqn{g}.
                                #' @param X_prime Numeric, normalized risk driver \eqn{X_t^{\prime} \in (0, 1)}.
                                g_prime = function(X_prime){
                                  private$..g_prime(X_prime)
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
                                  cat("------------------------ \033[1;35m boundTransform \033[0m ------------------------", "\n")
                                  cat(paste0("Link: \033[1;32m ", private$..link, "\033[0m \n"))
                                  cat(paste0("iY(y): \033[1;32m ",  alpha_, "\033[0m + \033[1;32m", beta_, "\033[0m exp(-exp(Y)) \n"))
                                  cat(paste0("Y(x): log(log(\033[1;32m",  beta_, "\033[0m) - log(x- \033[1;32m", alpha_, "\033[0m))"))
                                }
                              ),
                              private = list(
                                version = "1.0.1",
                                ..alpha = NA,
                                ..beta = NA,
                                ..link = NA,
                                ..g = NA,
                                ..ig = NA,
                                ..g_prime = NA,
                                ..monotonicity = ""
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
                                #' @field monotonicity Character, type of monotonicity of \eqn{g}.
                                monotonicity = function(){
                                  private$..monotonicity
                                },
                                #' @field link Character, name of the link function \eqn{g}.
                                link = function(){
                                  private$..link
                                }
                              )
)
