#' Exponential, Gaussian and Spherical isotropic spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param phi Numeric scalar, parameter.
#' @details Version 1.0.0. Implement the functions:
#' \deqn{\rho_{\exp}(h, \phi^{\exp}) = \exp\left(\frac{h}{\phi^{\exp}}\right)}
#' \deqn{\rho_{\text{gau}}(h, \phi^{\exp}) = \exp\left\{\left(\frac{h}{\phi^{\text{gau}}}\right)^2 \right\}}
#' \deqn{\rho_{\text{sph}}(h, \phi^{\text{sph}}) =  1 - \frac{3}{2}\frac{h}{\phi^{\text{sph}}} + \frac{1}{2}\left(\frac{h}{\phi^{\text{sph}}}\right)^3}
#' @rdname sp_cor_isotr
#' @name sp_cor_isotr
#' @export
sp_cor_isotr_exp <- function(h, phi){
  exp(-h/phi)
}

#' @rdname sp_cor_isotr
#' @name sp_cor_isotr
#' @export
sp_cor_isotr_gau <- function(h, phi){
  exp(-(h/phi)^2)
}

#' @rdname sp_cor_isotr
#' @name sp_cor_isotr
#' @export
sp_cor_isotr_sph <- function(h, phi){
  ifelse(h > phi, 0, 1 - 1.5*(h/phi) + 0.5*(h/phi)^3)
}

#' Exponential and Spherical anisotropic spatial correlation
#'
#' @param h1 Matrix of vector of distances.
#' @param h2 Matrix of vector of distances.
#' @param phi Numeric vector, parameters.
#' @details Version 1.0.0. Implement the functions:
#' \deqn{\rho_{\exp}(h_1, h_2, \phi^{\exp}) = \exp\left(\frac{1}{\phi_1^{\exp}} \sqrt{\phi_2^{\exp} h_1^2 + h_2^2}\right)}
#' \deqn{\rho_{\text{sph}}(h_1, h_2, \phi^{\text{sph}}) =  1 - \frac{3}{2}\frac{\sqrt{\phi_2^{\text{sph}} h_1^2 + h_2}}{\phi_1^{\text{sph}}} + \frac{1}{2}\left(\frac{\sqrt{\phi_2^{\text{sph}} h_1^2 + h_2}}{\phi_1^{\text{sph}}}\right)^3}
#' @rdname sp_cor_aniso
#' @name sp_cor_aniso
#' @export
sp_cor_aniso_exp <- function(h1, h2, phi){
  exp(-1/phi[1] * sqrt(phi[2] * h1^2 + h2^2))
}

#' @rdname sp_cor_aniso
#' @name sp_cor_aniso
#' @export
sp_cor_aniso_sph <- function(h1, h2, phi){
  1 - 1.5 * sqrt(phi[2] * h1^2 + h2^2)/phi[1] + 0.5 * (sqrt(phi[2] * h1^2 + h2^2)/phi[1])^3
}
