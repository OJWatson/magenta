#' Equilibrium list creation
#'
#' \code{equilibrium_ss_create} creates a starting state from
#'   the output of \code{equilibrium_init_create}
#'
#' @param eqInit Equilibrium initialisation state as
#'   output from \code{equilibrium_init_create}
#'
#' @keywords internal
equilibrium_ss_create <- function(eqInit) {
  ## create equilibrium state for return
  equilibrium_state <- list(
    "age_brackets" = eqInit$age_bounds,
    "het_brackets" = eqInit$het_bounds,
    "Smat" = eqInit$init_S[, , 1],
    "Dmat" = eqInit$init_D[, , 1],
    "Amat" = eqInit$init_A[, , 1],
    "Umat" = eqInit$init_U[, , 1],
    "Tmat" = eqInit$init_T[, , 1],
    "Pmat" = eqInit$init_P[, , 1],
    "IBmat" = eqInit$init_IB[, , 1],
    "ICAmat" = eqInit$init_ICA[, , 1],
    "ICMmat" = eqInit$init_ICM[, , 1],
    "IDmat" = eqInit$init_ID[, , 1],
    "FOI" = eqInit$FOI,
    "phi" = eqInit$phi_eq,
    "Sv" = eqInit$init_Sv * eqInit$mv0,
    "Ev" = eqInit$init_Ev * eqInit$mv0,
    "Iv" = eqInit$init_Iv * eqInit$mv0,
    "ICM_Init" = eqInit$ICM_init_eq,
    "theta" = seasonal_profile(eqInit$admin, eqInit$country)
  )

  validate_eqss(new_param_list(equilibrium_state, "magenta_eqss"))
}

#' @rdname equilibrium_ss_create
#' @param x List to validate.
validate_eqss <- function(x) {
  what <- "eqSS"
  check_fields(x, cpp_contract$eqss, what)
  na <- length(x$age_brackets)
  nh <- length(x$het_brackets)
  check_number(x$age_brackets, paste0(what, "$age_brackets"), min = 0, n = NULL)
  check_number(x$het_brackets, paste0(what, "$het_brackets"), min = 0, n = NULL)
  mats <- c(
    "Smat", "Dmat", "Amat", "Umat", "Tmat", "Pmat", "IBmat", "ICAmat",
    "ICMmat", "IDmat", "FOI", "phi"
  )
  for (field in mats) {
    m <- x[[field]]
    if (!is.matrix(m) || !is.numeric(m) || anyNA(m) || !identical(dim(m), c(na, nh))) {
      stop(sprintf(
        "%s$%s must be a %d by %d numeric matrix (age by heterogeneity brackets)",
        what, field, na, nh
      ), call. = FALSE)
    }
  }
  for (field in c("Sv", "Ev", "Iv")) {
    check_number(x[[field]], paste0(what, "$", field), min = 0)
  }
  check_number(x$ICM_Init, paste0(what, "$ICM_Init"), min = 0, n = nh)
  check_number(x$theta, paste0(what, "$theta"), n = NULL)
  invisible(x)
}


#' Spatial Parameter List Creation
#'
#' \code{spl_create} creates the user required spatial
#' list to be passed to the initialisation function if needed.
#'
#' @param spatial_type Numeric for spatial type. Default = 0
#' @param human_importation_rate_vector Vector of human importations into this
#'   region. Default = NULL
#' @param mosquito_imporation_rate_vector Vector of mosquito infection
#'   importations into this region. Default = NULL
#' @param cotransmission_freq_vector Vector of size of cotransmissions. This is
#'   preallocated here and regenerated every simulation update length.
#'   Default = rep(1, 10000)
#' @param oocyst_freq_vector Vector of number of ooocysts. This is
#'   preallocated here and regenerated every simulation update length.
#'   Default = rep(1, 10000)
#' @param num_human_infs Number of human infections from last day. Required for
#'   metapopulation. Default = 0
#' @param num_mos_infs Number of mosquito infections from last day. Required for
#'   metapopulation. Default = 0
#' @param plaf Population Level Frequency. Default = rep(0.5, 24)
#' @param island_imports_plaf_linked_flag Boolean. Whether imported barcodes for
#'  island model are drawn dependent on other sites. Currently, if TRUE, imported
#'  barcodes will either be all 0 or will be 1 at all loci where plaf > 0,
#'  if \code{rbinom(1, 1, plaf[plaf>0][1])} is TRUE. I.e. the first loci greater
#'  than 0 is used to determine if the incoming barcode is 1 at all sites where
#'  plaf is greater than 0.
#' @keywords internal
spl_create <- function(spatial_type = 0,
                       human_importation_rate_vector = NULL,
                       mosquito_imporation_rate_vector = NULL,
                       cotransmission_freq_vector = rep(1, 10000),
                       oocyst_freq_vector = rep(1, 10000),
                       num_human_infs = 0,
                       num_mos_infs = 0,
                       plaf = rep(0.5, 24),
                       island_imports_plaf_linked_flag = FALSE) {
  # if non spatial
  if (spatial_type == 0) {
    spatial_list <- list(
      "spatial_type" = spatial_type,
      "cotransmission_freq_vector" = cotransmission_freq_vector,
      "oocyst_freq_vector" = oocyst_freq_vector,
      "imported_cotransmissions_events" = 0,
      "imported_oocyst_events" = 0,
      "plaf" = plaf,
      "island_imports_plaf_linked_flag" = island_imports_plaf_linked_flag
    )

    # if island spatial
  } else if (spatial_type == 1) {
    # first work out how frequent importation is:
    hum_imp_rate <- sum(human_importation_rate_vector, na.rm = TRUE)
    mos_imp_rate <- sum(mosquito_imporation_rate_vector, na.rm = TRUE)

    spatial_list <- list(
      "spatial_type" = spatial_type,
      "cotransmission_freq_vector" = cotransmission_freq_vector,
      "oocyst_freq_vector" = oocyst_freq_vector,
      "imported_cotransmissions_events" = hum_imp_rate,
      "imported_oocyst_events" = mos_imp_rate,
      "plaf" = plaf,
      "island_imports_plaf_linked_flag" = island_imports_plaf_linked_flag
    )

    # metapopulation spatial lists were never finished
  } else {
    stop("metapopulation spatial lists are not implemented")
  }

  validate_spatial_list(new_param_list(spatial_list, "magenta_spatial_list"))
}

#' @rdname spl_create
#' @param x List to validate.
validate_spatial_list <- function(x) {
  what <- "spatial_list"
  check_fields(x, cpp_contract$spatial, what)
  check_number(x$spatial_type, paste0(what, "$spatial_type"), min = 0, max = 1, integer = TRUE)
  check_number(x$cotransmission_freq_vector, paste0(what, "$cotransmission_freq_vector"),
    min = 1, n = NULL, integer = TRUE
  )
  check_number(x$oocyst_freq_vector, paste0(what, "$oocyst_freq_vector"),
    min = 1, n = NULL, integer = TRUE
  )
  check_number(x$imported_cotransmissions_events, paste0(what, "$imported_cotransmissions_events"),
    min = 0
  )
  check_number(x$imported_oocyst_events, paste0(what, "$imported_oocyst_events"), min = 0)
  check_number(x$plaf, paste0(what, "$plaf"), min = 0, max = 1, n = NULL)
  check_flag(x$island_imports_plaf_linked_flag, paste0(what, "$island_imports_plaf_linked_flag"))
  invisible(x)
}
