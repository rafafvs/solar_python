#' Plot a forecast from a solarModel object
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' df_n <- model$moments$conditional[23,]
#' solarModel_predict_plot(solarModel_predict(model, df_n, ci = 0.01))
#' @export
solarModel_predict_plot <-function(object, pdf_components = FALSE, type = "mix"){

  grid <- object$grid
  ci <- object$ci
  emp <- object$df_n
  type <- match.arg(type, choices = c("mix", "up", "dw"))
  vals_colors <- c(realized = "black", expected_value = "orange",seasonal = "red", bounds = "purple", fitted = "magenta")
  vals_labels <- c(realized = "Realized", expected_value = "Expectation", seasonal = "Seasonal",
                   bounds = paste0("CI (", (1-ci)*100, " %)"),  fitted = "Fitted")

  pdf_plot <-  ggplot()+
    # Mixture density
    geom_line(data = grid, aes(x, pdf_Rt_mix))
  if (pdf_components){
    pdf_plot <-  ggplot()+
      # Components densities (weighted by prior probs)
      geom_line(data = grid, aes(x, pdf_Rt_mix_dw), color = "green")+
      geom_line(data = grid, aes(x, pdf_Rt_mix_up), color = "red")
  }

  if (type == "mix") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_lo, y = 0, yend = pdf_ci_mix_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_hi, xend = ci_mix_hi, y = 0, yend = pdf_ci_mix_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_mix_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_mix_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_Rt, 0, color = "expected_value"), size = 3)
  } else if (type == "dw") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_lo, y = 0, yend = pdf_ci_up_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_hi, xend = ci_up_hi, y = 0, yend = pdf_ci_up_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_up_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_up_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_Rt_up, 0, color = "expected_value"), size = 3)
  } else if (type == "up") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_dw_lo, xend = ci_dw_lo, y = 0, yend = pdf_ci_dw_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_dw_hi, xend = ci_dw_hi, y = 0, yend = pdf_ci_dw_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_dw_lo, xend = ci_dw_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_dw_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_dw_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_Rt_dw, 0, color = "expected_value"), size = 3)
  }
  type <- ifelse(type == "mix", "Mixture", ifelse(type == "up", "Sunny", "Cloudy"))
  plot_title <- paste0("Forecast: ", as.character(emp$date), " (", type, ")")
  pdf_plot+
    geom_point(data = emp, aes(Rt, 0, color = "realized"), size = 3)+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL, y = NULL, x = NULL, title = plot_title)+
    scale_color_manual(values = vals_colors, labels = vals_labels)
}

#' Plot scenarios from a solarScenario object
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' scenario <- solarScenario(model, nsim = 70)
#' solarScenario_plot(scenario)
#' @export
solarScenario_plot <- function(x, target = "GHI", nsim = 10, empiric = TRUE, ci = 0.05){
  scenario <- x
  df_emp <- scenario$emp
  df_sim <- scenario$sim
  df_sim <- dplyr::mutate(df_sim, data = purrr::map(data, ~.x[1:nsim,]))
  df_sim <- tidyr::unnest(df_sim, cols = c("data"))
  # Rename target columm
  df_sim$target <- df_sim[[target]]
  df_emp$target <- df_emp[[target]]

  # Plot scenarios
  plot_scenario <- ggplot()+
    geom_line(data = df_sim, aes(date, target, group = nsim, color = "simulated"), alpha = 0.1)

  # Add simulated value at risk
  if (!is.na(ci) & ci > 0) {
    VaR <- tidyr::unnest(scenario$sim, cols = c("data"))
    # Rename target columm
    VaR$target <- VaR[[target]]
    VaR <- VaR %>%
      group_by(date) %>%
      mutate(VaR_up = quantile(target, probs = 1-ci),
             VaR_dw = quantile(target, probs = ci)) %>%
      select(date, VaR_up, VaR_dw)

    plot_scenario <- plot_scenario +
      geom_line(data = VaR, aes(date, VaR_up, color = "VaR_up"), linewidth = 0.5)+
      geom_line(data = VaR, aes(date, VaR_dw, color = "VaR_dw"), linewidth = 0.5)
  }

 # Add realized GHI
 if (empiric) {
   plot_scenario <- plot_scenario +
     geom_point(data = df_emp, aes(date, target), color = "black", size = 1.5)+
     geom_point(data = df_emp, aes(date, target, color = "empiric"), size = 0.8)
 }
  vals_colors <- c(empiric = "orange", VaR_up = "red",VaR_dw = "red", simulated = "black")
  vals_labels <- c(empiric = "Realized", VaR_up = paste0("VaR (", (ci)*100, " %)"),
                   VaR_dw = paste0("VaR (", (1-ci)*100, " %)"),  simulated = paste0(nsim, " simulations"))
  plot_scenario+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL, y = NULL, x = NULL)+
    scale_color_manual(values = vals_colors, labels = vals_labels)

}
