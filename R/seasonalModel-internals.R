#' Create a Fourier formula
#'
#' @keywords seasonalModel
#' @noRd
seasonalModel_formula <- function(formula, order = 1, period = 365, sin = TRUE, cos = TRUE, t_idx = "n"){
  if (order == 0) {
    return(formula)
  }
  # Standard coefficients names
  coefs_names <- c()
  # Initial formula
  base_formula <- formula
  # Add sin terms
  if (sin) {
    # Standard name
    new_coef_name <- paste0("sin_", order, "_", period)
    # Check if it was not included
    if (!new_coef_name %in% attr(base_formula, "coef_names")){
      base_formula <- paste0(base_formula, " + ", "I(sin(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
      coefs_names <- c(coefs_names, new_coef_name)
    }
  }
  # Add cos terms
  if (cos) {
    # Standard name
    new_coef_name <- paste0("cos_", order, "_", period)
    # Check if it was not included
    if (!new_coef_name %in% attr(base_formula, "coef_names")){
      base_formula <- paste0(base_formula, " + ", "I(cos(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
      coefs_names <- c(coefs_names, new_coef_name)
    }
  }

  attr(base_formula, "coef_names") <- c(attr(formula, "coef_names"), coefs_names)
  return(base_formula)
}

#' Create a Fourier formula for the differential
#'
#' @keywords seasonalModel
#' @noRd
seasonalModel_formula_dt <- function(formula, order = 1, period = 365, sin = TRUE, cos = TRUE, t_idx = "n"){
  if (order == 0) {
    return(formula)
  }
  # Standard coefficients names
  coefs_names <- c()
  # Initial formula
  base_formula <- formula
  # Add sin terms
  if (sin) {
    # Standard name
    new_coef_name <- paste0("sin_", order, "_", period)
    # Check if it was not included
    if (!new_coef_name %in% attr(base_formula, "coef_names")){
      base_formula <- paste0(base_formula, " + ", "I(", eval(order), " * 2 * base::pi /", eval(period), " * cos(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
      coefs_names <- c(coefs_names, new_coef_name)
    }
  }
  # Add cos terms
  if (cos) {
    # Standard name
    new_coef_name <- paste0("cos_", order, "_", period)
    # Check if it was not included
    if (!new_coef_name %in% attr(base_formula, "coef_names")){
      base_formula <- paste0(base_formula, " + ", "I(", eval(order), " * 2 * base::pi /", eval(period), " * sin(2 * base::pi /", eval(period), " * ", t_idx, " * ", eval(order), "))")
      coefs_names <- c(coefs_names, new_coef_name)
    }
  }

  attr(base_formula, "coef_names") <- c(attr(formula, "coef_names"), coefs_names)
  return(base_formula)
}


#' From constraint to unconstrained parameters
#'
#' @keywords seasonalModel
#' @noRd
#' @export
seasonalModel_params_to_zeta <- function(b){
  # Compute the reparametrized coefficients
  b0_star <- log(b[1])
  b1_star <- atanh(sqrt(b[2]^2 + b[3]^2)/b[1])
  b2_star <- atan2(b[2], b[3])
  # Vector of parameters
  b_star <- c(b0_star, b1_star, b2_star)
  names(b_star) <- paste0(names(b), "_star")
  return(b_star)
}

#' From unconstrained to constraint parameters
#'
#' @keywords seasonalModel
#' @noRd
#' @export
seasonalModel_params_to_phi <- function(b_star){
  # Compute the reparametrized coefficients
  b0 <- exp(b_star[1])
  b1 <- b0 * tanh(b_star[2]) * sin(b_star[3])
  b2 <- b0 * tanh(b_star[2]) * cos(b_star[3])
  # Vector of parameters
  b <- c(b0, b1, b2)
  names(b) <- stringr::str_remove_all(names(b_star), "_star")
  return(b)
}


#' Jacobian from unconstrained to constraint parameters
#'
#' @keywords seasonalModel
#' @noRd
#' @export
seasonalModel_params_to_zeta_jacobian <- function(b_star){
  # Extract the coefficients
  b0_star <- b_star[1]
  b1_star <- b_star[2]
  b2_star <- b_star[3]
  # Initialize
  J <- matrix(0, 3, 3)
  # Derivatives of b0 wrt b0, b1, b2
  J[1,1] <- exp(b0_star)                                  # d_b0_d_b0 = b0
  J[1,2] <- 0                                             # d_b0_d_b1
  J[1,3] <- 0                                             # d_b0_d_b2
  # Derivatives of b1 wrt b0, b1, b2
  J[2,1] <- J[1,1] * tanh(b1_star) * sin(b2_star)         # d_b1_d_b0 = b1
  J[2,2] <- J[1,1] * (1 - tanh(b1_star)^2) * sin(b2_star) # d_b1_d_b1
  J[2,3] <- J[1,1] * tanh(b1_star) * cos(b2_star)         # d_b1_d_b2 = b2
  # Derivatives of b2 wrt b0, b1, b2
  J[3,1] <- J[2,3]                                        # d_b2_d_b0 = b2
  J[3,2] <- J[1,1] * (1 - tanh(b1_star)^2) * cos(b2_star) # d_b2_d_b1
  J[3,3] <- -J[2,1]                                       # d_b2_d_b2 = -b1
  colnames(J) <- c("d_b0_star", "d_b1_star", "d_b2_star")
  rownames(J) <- c("b0", "b1", "b2")
  return(J)
}
