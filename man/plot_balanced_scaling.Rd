% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_balanced_scaling.R
\name{plot_balanced_scaling}
\alias{plot_balanced_scaling}
\title{Visual comparison of original and balanced data using scatterplots and marginal histograms}
\usage{
plot_balanced_scaling(
  original_df,
  balanced_df,
  var_x,
  var_y,
  model_type = c("power", "exp", "linear"),
  binwidth,
  base = 10,
  tag_letters = c("a", "b"),
  x_lab = NULL,
  y_lab = NULL
)
}
\arguments{
\item{original_df}{A data frame with the original observations.}

\item{balanced_df}{A data frame returned by \code{balanced_scaling()}, which
must contain a logical column \code{Resampled}.}

\item{var_x, var_y}{Unquoted column names for the predictor and response.}

\item{model_type}{\code{"power"} (default), \code{"exp"}, or \code{"linear"}.}

\item{binwidth}{Histogram bin width (on the binning scale). Can be called from bin summary (see "Examples").}

\item{base}{Logarithmic base when \code{scale = "log"}. Default \code{10}.}

\item{tag_letters}{Two letters used to label the panels.}

\item{x_lab, y_lab}{Optional axis titles.}
}
\value{
Invisibly returns a \strong{grob} that can be printed or saved.
}
\description{
\code{plot_balanced_scaling()} juxtaposes the \strong{raw} data with the
\strong{bootstrap‑balanced} sample returned by \code{balanced_scaling()}.
Each panel shows a scatterplot overlaid with an SMA line and is wrapped
in marginal histograms that reveal how resampling fills sparse bins.
}
\examples{
library(ggplot2)

data(xylem_scaling_simulation_dataset)

res<- balanced_scaling(
  data        = xylem_scaling_simulation_dataset,
  var_x       = L,
  var_y       = DAVG,
  min_per_bin = 100,
  n_boot      = 10,
  seed        = 1,
  model_type  = "power"
)
plot_balanced_scaling(
  original_df = xylem_scaling_simulation_dataset,
  balanced_df = res$first_boot,
  var_x       = L,
  var_y       = DAVG,
  model_type  = "power",
  binwidth    = res$bins$summary$bin_width[1],
  base        = 10,
  tag_letters = c("a", "b"),
  x_lab       = NULL,
  y_lab       = NULL
)

}
\references{
Simovic, M., & Michaletz, S.T. (2025). \emph{Harnessing the Full Power of Data
to Characterise Biological Scaling Relationships.} \strong{Global Ecology and
Biogeography}, 34(2). \url{https://doi.org/10.1111/geb.70019}
}
\author{
Simovic, M. \href{mailto:milos.simovic@botany.ubc.ca}{milos.simovic@botany.ubc.ca};
Michaletz, S.T. \href{mailto:sean.michaletz@ubc.ca}{sean.michaletz@ubc.ca}
}
