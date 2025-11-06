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

