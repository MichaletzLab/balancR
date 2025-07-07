#' Visual comparison of original and balanced data using scatterplots and marginal histograms
#'
#' `plot_balanced_scaling()` juxtaposes the **raw** data with the
#' **bootstrap‑balanced** sample returned by `balanced_scaling()`.
#' Each panel shows a scatterplot overlaid with an SMA line and is wrapped
#' in marginal histograms that reveal how resampling fills sparse bins.
#'
#' @inheritParams balanced_scaling
#' @param original_df A data frame with the original observations.
#' @param balanced_df A data frame returned by `balanced_scaling()`, which
#'   must contain a logical column `Resampled`.
#' @param binwidth Histogram bin width (on the binning scale). Can be called from bin summary (see "Examples").
#' @param tag_letters Two letters used to label the panels.
#' @param x_lab,y_lab Optional axis titles.
#'
#' @return Invisibly returns a **grob** that can be printed or saved.
#'
#' @author
#' Simovic, M. <milos.simovic@botany.ubc.ca>;
#' Michaletz, S.T. <sean.michaletz@ubc.ca>
#'
#' @references
#' Simovic, M., & Michaletz, S.T. (2025). *Harnessing the Full Power of Data
#' to Characterise Biological Scaling Relationships.* **Global Ecology and
#' Biogeography**, 34(2). <https://doi.org/10.1111/geb.70019>
#'
#' @examples
#' library(ggplot2)
#'
#' data(xylem_scaling_simulation_dataset)
#'
#' res<- balanced_scaling(
#'   data        = xylem_scaling_simulation_dataset,
#'   var_x       = L,
#'   var_y       = DAVG,
#'   min_per_bin = 100,
#'   n_boot      = 10,
#'   seed        = 1,
#'   model_type  = "power"
#' )

#' plot_balanced_scaling(
#'   original_df = xylem_scaling_simulation_dataset,
#'   balanced_df = res$first_boot,
#'   var_x       = L,
#'   var_y       = DAVG,
#'   model_type  = "power",
#'   binwidth    = res$bins$summary$bin_width[1],
#'   base        = 10,
#'   tag_letters = c("a", "b"),
#'   x_lab       = NULL,
#'   y_lab       = NULL
#' )
#'
#' @importFrom stats as.formula
#' @importFrom rlang .data ensym as_name
#' @export
plot_balanced_scaling <- function(original_df,
                                  balanced_df,
                                  var_x,
                                  var_y,
                                  model_type  = c("power", "exp", "linear"),
                                  binwidth,
                                  base        = 10,
                                  tag_letters = c("a", "b"),
                                  x_lab       = NULL,
                                  y_lab       = NULL) {

  model_type <- match.arg(model_type)

  var_x_sym <- rlang::ensym(var_x)
  var_y_sym <- rlang::ensym(var_y)
  var_x_chr <- rlang::as_name(var_x_sym)
  var_y_chr <- rlang::as_name(var_y_sym)

  if (is.null(x_lab)) x_lab <- var_x_chr
  if (is.null(y_lab)) y_lab <- var_y_chr

  ##------------------------------------------------------------------##
  ## 1. Back‑transform balanced data so both sets share the same scale
  ##------------------------------------------------------------------##
  orig_df_plot     <- original_df
  balanced_df_plot <- balanced_df

  if (model_type == "power") {
    balanced_df_plot$x_plot <- base ^ balanced_df_plot[[var_x_chr]]
    balanced_df_plot$y_plot <- base ^ balanced_df_plot[[var_y_chr]]
  } else if (model_type == "exp") {
    balanced_df_plot$x_plot <- balanced_df_plot[[var_x_chr]]
    balanced_df_plot$y_plot <- base ^ balanced_df_plot[[var_y_chr]]
  } else {                         # linear
    balanced_df_plot$x_plot <- balanced_df_plot[[var_x_chr]]
    balanced_df_plot$y_plot <- balanced_df_plot[[var_y_chr]]
  }

  orig_df_plot$x_plot <- orig_df_plot[[var_x_chr]]
  orig_df_plot$y_plot <- orig_df_plot[[var_y_chr]]

  ##------------------------------------------------------------------##
  ## 2. Original panel
  ##------------------------------------------------------------------##
  p_orig <- ggplot2::ggplot(
    orig_df_plot,
    ggplot2::aes(x = .data$x_plot, y = .data$y_plot)
  ) +
    ggplot2::geom_point(alpha = 0.10, colour = "navy") +
    ggpmisc::stat_ma_line(method = "SMA", colour = "darkorange", se = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_classic() +
    ggplot2::labs(tag = tag_letters[1]) +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.title   = ggplot2::element_blank(),
      axis.text    = ggplot2::element_text(size = 14),
      plot.tag     = ggplot2::element_text(size = 14, face = "bold")
    )

  p_orig <- ggExtra::ggMarginal(
    p_orig,
    type     = "histogram",
    binwidth = binwidth,
    colour   = "navy",
    fill     = "navy",
    alpha    = 0.10
  )

  ##------------------------------------------------------------------##
  ## 3. Balanced panel
  ##------------------------------------------------------------------##
  p_bal <- ggplot2::ggplot(
    balanced_df_plot,
    ggplot2::aes(x = .data$x_plot, y = .data$y_plot)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$Resampled, alpha = .data$Resampled)
    ) +
    ggplot2::scale_colour_manual(values = c("navy", "cyan2")) +
    ggplot2::scale_alpha_manual(values = c(0.10, 0.03)) +
    ggpmisc::stat_ma_line(method = "SMA", colour = "darkorange", se = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_classic() +
    ggplot2::labs(tag = tag_letters[2]) +
    ggplot2::theme(
      aspect.ratio   = 1,
      axis.title     = ggplot2::element_blank(),
      axis.text      = ggplot2::element_text(size = 14),
      plot.tag       = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "none"
    )

  p_bal <- ggExtra::ggMarginal(
    p_bal,
    type        = "histogram",
    binwidth    = binwidth,
    groupColour = TRUE,
    groupFill   = TRUE,
    alpha       = 0.10
  )

  ##------------------------------------------------------------------##
  ## 4. Assemble & return
  ##------------------------------------------------------------------##
  figure <- gridExtra::grid.arrange(
    p_orig, p_bal, nrow = 1,
    left   = grid::textGrob(y_lab, rot = 90, gp = grid::gpar(fontsize = 14)),
    bottom = grid::textGrob(x_lab,       gp = grid::gpar(fontsize = 14))
  )

  invisible(figure)
}
