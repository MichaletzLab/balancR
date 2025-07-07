#' Compare SMA slopes or elevations between original and balanced data
#'
#' `compare_coefficients()` binds an **original** data set to its
#' **balanced/bootstrapped** counterpart, labels each with a `Dataset`
#' factor, and fits a *grouped* standardised major axis (SMA) model to test
#' whether the **slope** or the **elevation** (intercept) differs between
#' groups, and whether that coefficient differs from a user‑supplied
#' reference value.
#'
#' @param original  Data frame with the **unscaled** data.
#' @param balanced  Data frame with the balanced data (e.g.
#'   `balanced_scaling(... )$first_boot`).
#' @param var_x,var_y Unquoted column names for the predictor and response.
#' @param which_coefficient `"slope"` (default) or `"intercept"`.
#' @param model_type  `"power"` (default), `"exp"`, or `"linear"`.
#' @param test_value  Numeric; value to compare the chosen coefficient
#'   against (optional).
#'
#' @return A `"sma.summary"` object (see [smatr::sma()]).
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
#' @examples
#' if (requireNamespace("smatr", quietly = TRUE)) {
#'   data(xylem_scaling_simulation_dataset)
#'
#'   bal <- balanced_scaling(
#'     data        = xylem_scaling_simulation_dataset,
#'     var_x       = L,
#'     var_y       = DAVG,
#'     min_per_bin = 100,
#'     n_boot      = 1,
#'     seed        = 1,
#'     model_type  = "power"
#'   )$first_boot
#'
#'   #Compare slopes between balanced and imbalanced datasets
#'
#'   compare_coefficients(
#'     original  = xylem_scaling_simulation_dataset,
#'     balanced  = bal,
#'     var_x     = L,
#'     var_y     = DAVG,
#'     which_coefficient = "slope",
#'     model_type        = "power"
#'   )
#'
#'   #Test if balancing data results in better agreement with WBE model
#'
#'   compare_coefficients(
#'     original  = xylem_scaling_simulation_dataset,
#'     balanced  = bal,
#'     var_x     = L,
#'     var_y     = DAVG,
#'     which_coefficient = "slope",
#'     model_type        = "power",
#'     test_value        = 0.25
#'   )
#' }
#'
#' @importFrom stats as.formula
#' @export
compare_coefficients <- function(original,
                                 balanced,
                                 var_x,
                                 var_y,
                                 which_coefficient = c("slope", "intercept"),
                                 model_type        = c("power", "exp", "linear"),
                                 test_value        = NA) {

  if (missing(balanced)) {
    stop("`balanced` must be supplied (e.g. balanced_scaling(... )$first_boot).",
         call. = FALSE)
  }

  if (!requireNamespace("smatr", quietly = TRUE)) {
    stop("Please install the 'smatr' package.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package.")
  }

  which_coefficient <- match.arg(which_coefficient)
  model_type        <- match.arg(model_type)

  var_x_sym <- rlang::ensym(var_x)
  var_y_sym <- rlang::ensym(var_y)
  var_x_chr <- rlang::as_name(var_x_sym)
  var_y_chr <- rlang::as_name(var_y_sym)

  original <- dplyr::as_tibble(original)
  balanced <- dplyr::as_tibble(balanced)

  if ("Dataset" %in% c(names(original), names(balanced))) {
    stop("Column 'Dataset' already exists in one of the inputs.", call. = FALSE)
  }

  ## -------------------------------------------------------------------- ##
  ##  Transform ORIGINAL data only
  ## -------------------------------------------------------------------- ##
  if (model_type == "exp") {
    original[[var_y_chr]] <- log(original[[var_y_chr]])
  } else if (model_type == "power") {
    original[[var_x_chr]] <- log10(original[[var_x_chr]])
    original[[var_y_chr]] <- log10(original[[var_y_chr]])
  }

  ## -------------------------------------------------------------------- ##
  ##  Combine & label
  ## -------------------------------------------------------------------- ##
  original <- dplyr::mutate(original, Dataset = "Original")
  balanced <- dplyr::mutate(balanced, Dataset = "Balanced")
  combined <- dplyr::bind_rows(original, balanced)

  ## -------------------------------------------------------------------- ##
  ##  Build SMA formula
  ## -------------------------------------------------------------------- ##
  op   <- if (which_coefficient == "slope") "*" else "+"
  fmla <- stats::as.formula(paste(var_y_chr, "~", var_x_chr, op, "Dataset"))

  slope_test <- if (which_coefficient == "slope")     test_value else NA
  elev_test  <- if (which_coefficient == "intercept") test_value else NA

  fit <- smatr::sma(
    fmla,
    data       = combined,
    method     = "SMA",
    log        = "",
    slope.test = slope_test,
    elev.test  = elev_test
  )

  summary(fit)
}
