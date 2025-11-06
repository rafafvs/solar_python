#' Evaluate VaR test
#'
#' @rdname VaR_test
#' @name VaR_test
#' @export
VaR_test <- function(et, VaR, ci){
  n <- length(et)
  # Test Statistics
  T_stat <- sum((et - VaR) / sqrt(VaR*(1-VaR)))/sqrt(n)
  # P.value
  p.value <- 2*pnorm(-abs(T_stat))
  # Critical values
  rejection_lev <- abs(qnorm(ci/2))
  tibble(
    n = n,
    ci = ci,
    stat = T_stat,
    VaR_hat = mean(et),
    VaR = VaR,
    rejection_lev = rejection_lev,
    p.value = p.value,
    H0 = ifelse(abs(stat) > rejection_lev, "Rejected", "Non-Rejected"),
  )
}
