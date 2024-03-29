% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainFunctions.R
\name{param_list_simulation_update_create}
\alias{param_list_simulation_update_create}
\title{Parameter List creation for magenta simulation updating}
\usage{
param_list_simulation_update_create(
  years = 1,
  ft = 0.4,
  mu_vec = NULL,
  fv_vec = NULL,
  statePtr,
  spatial_list,
  drug_list,
  barcode_list
)
}
\arguments{
\item{years}{Length of simulation. Default = 1}

\item{ft}{Treatments seeking value}

\item{mu_vec}{Vector of mosquito mortalities for each day with years. Default = NULL, which
will result in rep(0.132,floor(years*365))}

\item{fv_vec}{Vector of mosquito bitings for each day with years. Default = NULL, which
will result in rep(1/3,floor(years*365))}

\item{statePtr}{Pointer for current model state as return by \code{simulation_R}$Ptr}

\item{spatial_list}{Spatial list}

\item{drug_list}{Drug list}

\item{barcode_list}{Barcode parameter list}
}
\description{
\code{param_list_simulation_update_create} creates suitable parameter list for
\code{simulation_R} for continuing a simulation from memory within the active
session.
}
