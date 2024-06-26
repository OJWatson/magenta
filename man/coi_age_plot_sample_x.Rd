% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{coi_age_plot_sample_x}
\alias{coi_age_plot_sample_x}
\title{Plot COI against age from random sample}
\usage{
coi_age_plot_sample_x(
  sample_coi_out,
  x,
  span = 0.6,
  ylimmax = NULL,
  xlimmax = NULL,
  max_coi = 25
)
}
\arguments{
\item{sample_coi_out}{Output of sample_coi}

\item{x}{Which sample to plot}

\item{span}{Smoothing parameter for loess}

\item{ylimmax}{ylim max}

\item{xlimmax}{xlim max}

\item{max_coi}{max_coi}
}
\description{
\code{sample_coi} samples from a human save according to some age distirbution
}
\keyword{internal}
