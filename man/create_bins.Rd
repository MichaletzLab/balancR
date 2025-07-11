% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_bins.R
\name{create_bins}
\alias{create_bins}
\title{Flexible equal-width binning with a minimum sample size}
\usage{
create_bins(
  data,
  var,
  min_per_bin = 100,
  scale = "log",
  base = 10,
  right = TRUE
)
}
\arguments{
\item{data}{A data frame or tibble.}

\item{var}{Unquoted \strong{numeric} column to bin (tidy-eval).}

\item{min_per_bin}{Minimum number of observations per bin. Default \code{100}.}

\item{scale}{Binning scale: \code{"log"} (default) or \code{"linear"}.}

\item{base}{Logarithmic base when \code{scale = "log"}. Default \code{10}.}

\item{right}{Logical; should bins be right-closed? Default \code{TRUE}.}
}
\value{
A \strong{list} with two tibbles:
\describe{
\item{\code{data}}{The original data with an added \code{bin_class} column.}
\item{\code{summary}}{Per-bin counts, rows to bootstrap (\code{bootstrap_count}),
and the bin width.}
}
}
\description{
\code{create_bins()} finds the \emph{largest} number of equal-width bins
(on a \strong{log} or \strong{linear} axis) such that \strong{every} bin still
contains at least \code{min_per_bin} observations.
It then reports, for each bin, how many additional rows would need
to be bootstrapped so that all bins end up the same size.
}
\examples{
data(xylem_scaling_simulation_dataset)

bins <- create_bins(
  data         = xylem_scaling_simulation_dataset,
  var          = L,
  min_per_bin  = 100,
  scale        = "log",
  base         = 10
)

head(bins$data)      # Binned dataset
bins$summary         # Counts & bootstrap totals

}
\references{
Simovic M., & Michaletz S.T. (2025). \emph{Harnessing the Full Power of Data
to Characterise Biological Scaling Relationships.} Global Ecology and
Biogeography, 34(2). \url{https://doi.org/10.1111/geb.70019}
}
\seealso{
\code{\link[=balanced_scaling]{balanced_scaling()}}, \code{\link[=plot_balanced_scaling]{plot_balanced_scaling()}}
}
\author{
Simovic, M. \href{mailto:milos.simovic@botany.ubc.ca}{milos.simovic@botany.ubc.ca};
Michaletz, S.T. \href{mailto:sean.michaletz@ubc.ca}{sean.michaletz@ubc.ca}
}
