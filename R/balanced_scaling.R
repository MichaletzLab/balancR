#' Data balancing across log- or linear-scaled bins and fitting via SMA regression
#'
#' `balanced_scaling()` partitions a data set into equal‑width bins on a
#' log or linear axis, upsamples each bin so they contain the same number
#' of observations, and then fits a **standardised major axis (SMA)**
#' model to every balanced bootstrap replicate.
#'
#' Three model types are supported:
#' \describe{
#'   \item{`"power"`}{Power‑law model *y = a xᵇ* (log10–log10).}
#'   \item{`"exp"`}{Exponential model *y = a exp(b x)* (log–linear).}
#'   \item{`"linear"`}{Ordinary linear model *y = a + b x*.}
#' }
#'
#' @inheritParams create_bins
#' @param var_x,var_y Unquoted column names for the predictor and response.
#' @param n_boot Number of bootstrap iterations.  Default `100`.
#' @param seed  Base seed for reproducibility; iteration *i* uses
#'   `seed + i`.  Set `NULL` for no seeding.
#' @param model_type `"power"` (default), `"exp"`, or `"linear"`.
#'
#' @return A **list** with:
#' \describe{
#'   \item{`stats`}{Regression statistics, including r², p-value, slope, intercept (i.e., elevation).}
#'   \item{`first_boot`}{The first bootstrap-balanced dataset generated via the function. Useful for plotting or statistical comparisons with original, imbalanced data.}
#'   \item{`bins`}{The exact output from [create_bins()] —
#'                 a list with `data` (binned input rows) and
#'                 `summary` (one‑row‑per‑bin metadata including
#'                 `bin_width`, `bin_class`, `bootstrap_count`, etc.).}
#' }
#'
#' @author
#' Simovic, M. <milos.simovic@botany.ubc.ca>;
#' Michaletz, S.T. <sean.michaletz@ubc.ca>
#'
#' @references
#' Simovic, M., & Michaletz, S.T. (2025). *Harnessing the Full Power of Data
#' to Characterise Biological Scaling Relationships.* **Global Ecology and
#' Biogeography**, 34(2). <https://doi.org/10.1111/geb.70019>
#' Warton, D.I., Duursma, R.A., Falster, D.S., & Taskinen, S. (2012).
#' *smatr 3 – an R package for estimation and inference about allometric
#' lines.* **Methods in Ecology and Evolution**, 3(2), 257–259.
#' <https://doi.org/10.1111/j.2041-210X.2011.00153.x>
#'
#' @seealso [create_bins()], [smatr::sma()]
#'
#' @examples
#' if (requireNamespace("smatr", quietly = TRUE)) {
#'   data(xylem_scaling_simulation_dataset)
#'   res <- balanced_scaling(
#'     data        = xylem_scaling_simulation_dataset,
#'     var_x       = L,
#'     var_y       = DAVG,
#'     min_per_bin = 100,
#'     n_boot      = 10,
#'     seed        = 1,
#'     model_type  = "power"
#'   )
#'   head(res$stats)
#' }
#'
#' @importFrom stats as.formula
#' @export

balanced_scaling <- function(data,
                             var_x,
                             var_y,
                             min_per_bin = 100,
                             n_boot      = 100,
                             base        = 10,
                             seed        = 1,
                             model_type  = c("power", "exp", "linear")) {

  if (!requireNamespace("smatr",  quietly = TRUE)) {
    stop("Please install the 'smatr' package.")
  }
  if (!requireNamespace("dplyr",  quietly = TRUE)) {
    stop("Please install the 'dplyr' package.")
  }

  ## -------------------------------------------------------------------------
  ##  Set‑up
  ## -------------------------------------------------------------------------
  var_x_sym <- rlang::ensym(var_x)
  var_y_sym <- rlang::ensym(var_y)
  var_x_chr <- rlang::as_name(var_x_sym)
  var_y_chr <- rlang::as_name(var_y_sym)
  fmla      <- stats::as.formula(paste(var_y_chr, "~", var_x_chr))

  model_type <- match.arg(model_type)
  scale      <- if (model_type == "power") "log" else "linear"

  ## -------------------------------------------------------------------------
  ##  Binning
  ## -------------------------------------------------------------------------
  bins <- create_bins(
    data         = data,
    var          = !!var_x_sym,
    min_per_bin  = min_per_bin,
    scale        = scale,
    base         = base
  )

  df          <- bins$data
  bin_summary <- bins$summary

  ## -------------------------------------------------------------------------
  ##  Transform variables (if needed)
  ## -------------------------------------------------------------------------
  if (model_type == "exp") {
    df[[var_y_chr]] <- log(df[[var_y_chr]])
  } else if (model_type == "power") {
    df[[var_x_chr]] <- log10(df[[var_x_chr]])
    df[[var_y_chr]] <- log10(df[[var_y_chr]])
  }

  ## -------------------------------------------------------------------------
  ##  Bootstrap loop
  ## -------------------------------------------------------------------------
  stats_list <- vector("list", n_boot)
  first_boot <- NULL

  for (i in seq_len(n_boot)) {
    if (!is.null(seed)) set.seed(seed + i)

    df_boot <- dplyr::group_by(df, .data$bin_class) |>
      dplyr::group_modify(~{
        n_add <- bin_summary$bootstrap_count[
          bin_summary$bin_class == .y$bin_class
        ]
        original <- dplyr::mutate(.x, Resampled = FALSE)
        if (n_add > 0) {
          resampled <- dplyr::mutate(
            dplyr::slice_sample(.x, n = n_add, replace = TRUE),
            Resampled = TRUE
          )
          dplyr::bind_rows(original, resampled)
        } else {
          original
        }
      }) |>
      dplyr::ungroup()

    if (i == 1L) first_boot <- df_boot

    fit  <- smatr::sma(fmla, data = df_boot, method = "SMA", log = "")
    gsum <- fit$groupsummary

    stats_list[[i]] <- data.frame(
      iter         = i,
      slope        = gsum$Slope,
      slope_lo     = gsum$Slope_lowCI,
      slope_hi     = gsum$Slope_highCI,
      intercept    = gsum$Int,
      intercept_lo = gsum$Int_lowCI,
      intercept_hi = gsum$Int_highCI,
      r2           = gsum$r2,
      pval         = gsum$pval,
      n            = gsum$n
    )
  }

  ## -------------------------------------------------------------------------
  ##  Return
  ## -------------------------------------------------------------------------
  list(
    stats      = dplyr::bind_rows(stats_list),
    first_boot = first_boot,
    bins       = bins
  )
}
