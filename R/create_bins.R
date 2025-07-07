#' Flexible equal-width binning with a minimum sample size
#'
#' `create_bins()` finds the *largest* number of equal-width bins
#' (on a **log** or **linear** axis) such that **every** bin still
#' contains at least `min_per_bin` observations.
#' It then reports, for each bin, how many additional rows would need
#' to be bootstrapped so that all bins end up the same size.
#'
#' @param data A data frame or tibble.
#' @param var  Unquoted **numeric** column to bin (tidy-eval).
#' @param min_per_bin Minimum number of observations per bin. Default `100`.
#' @param scale Binning scale: `"log"` (default) or `"linear"`.
#' @param base  Logarithmic base when `scale = "log"`. Default `10`.
#' @param right Logical; should bins be right-closed? Default `TRUE`.
#'
#' @return A **list** with two tibbles:
#' \describe{
#'   \item{`data`}{The original data with an added `bin_class` column.}
#'   \item{`summary`}{Per-bin counts, rows to bootstrap (`bootstrap_count`),
#'                    and the bin width.}
#' }
#'
#' @seealso [balanced_scaling()], [plot_balanced_scaling()]
#'
#' @author
#' Simovic, M. <milos.simovic@botany.ubc.ca>;
#' Michaletz, S.T. <sean.michaletz@ubc.ca>
#'
#' @references
#' Simovic M., & Michaletz S.T. (2025). *Harnessing the Full Power of Data
#' to Characterise Biological Scaling Relationships.* Global Ecology and
#' Biogeography, 34(2). <https://doi.org/10.1111/geb.70019>
#'
#' @examples
#' data(xylem_scaling_simulation_dataset)
#'
#' bins <- create_bins(
#'   data         = xylem_scaling_simulation_dataset,
#'   var          = L,
#'   min_per_bin  = 100,
#'   scale        = "log",
#'   base         = 10
#' )
#'
#' head(bins$data)      # Binned dataset
#' bins$summary         # Counts & bootstrap totals
#'
#' @importFrom rlang :=
#' @export

create_bins <- function(data, var, min_per_bin = 100,
                        scale = "log", base = 10, right = TRUE) {

  var_quo  <- rlang::enquo(var)
  var_name <- rlang::as_name(var_quo)

  x <- dplyr::pull(data, !!var_quo)
  if (!is.numeric(x)) {
    stop("`var` must refer to a numeric column.", call. = FALSE)
  }
  x <- x[!is.na(x)]
  if (length(x) < min_per_bin) {
    stop("Not enough non NA values to satisfy `min_per_bin`.", call. = FALSE)
  }
  if (scale == "log" && any(x <= 0)) {
    stop("Log binning requires all values > 0.", call. = FALSE)
  }

  transformed <- if (scale == "log") log(x, base = base) else x
  min_val <- min(transformed)
  max_val <- max(transformed)

  max_bins <- max(1L, floor(length(x) / min_per_bin))
  num_bins <- max_bins

  repeat {
    bin_size     <- (max_val - min_val) / num_bins
    breaks_trans <- min_val + bin_size * 0:num_bins
    breaks_trans[1]                 <- min_val - 1e-6   # open on left
    breaks_trans[length(breaks_trans)] <- max_val + 1e-6

    bins <- cut(
      transformed,
      breaks         = breaks_trans,
      include.lowest = TRUE,
      right          = right
    )

    counts <- tabulate(bins, nbins = length(breaks_trans) - 1)
    if (min(counts) >= min_per_bin || num_bins == 1L) break
    num_bins <- num_bins - 1L
  }

  # Transform breaks back to original scale for labelling
  breaks <- if (scale == "log") base ^ breaks_trans else breaks_trans
  labels <- paste(breaks[-length(breaks)], "-", breaks[-1])
  class_col <- "bin_class"

  binned_data <- dplyr::mutate(
    data,
    !!class_col := cut(
      !!var_quo,
      breaks         = breaks,
      include.lowest = TRUE,
      right          = right,
      labels         = labels
    )
  )

  bin_widths <- if (scale == "log") diff(breaks_trans) else diff(breaks)

  binned_data_summary <- binned_data |>
    dplyr::group_by(.data[[class_col]]) |>
    dplyr::summarise(original_count = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      bootstrap_count = max(.data$original_count) - .data$original_count,
      bin_width       = bin_widths
    )

  list(
    data    = binned_data,
    summary = binned_data_summary
  )
}
