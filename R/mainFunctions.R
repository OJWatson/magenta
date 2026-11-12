#------------------------------------------------
#' Parameter List creation for magenta simulation initialisation
#'
#' \code{param_list_simulation_init_create} creates suitable parameter list for
#' \code{simulation_R} for the beginning of a simulation. Also takes an argument
#' for feeding in spatial parameters/data.
#'
#' @param N Population size. Default = 1e4
#' @param eqSS Output of \code{equilibrium_steady_state_create}
#' @param barcode_list List of barcode/genetic parameters
#' @param spatial_list Spatial parmeters to come in
#' @param housekeeping_list Housekeeping parameter list
#'   from \code{housekeeping_list_create}
#' @param drug_list Drug parameter list
#'   from \code{drug_list_create}
#' @param nmf_list Non malarial fever parameter list
#'   from \code{nmf_list_create}
#' @param vector_adaptation_list vector adaptation list
#'   from \code{vector_adaptation_list}
#' @param mpl model parameter list from \code{model_param_list_create}
#'
#' @export

param_list_simulation_init_create <- function(N = 1e+04, eqSS, barcode_list,
                                              spatial_list, housekeeping_list,
                                              drug_list, nmf_list,
                                              vector_adaptation_list,
                                              mpl) {
  ## CHECKS ##
  ## ---------------------------------------------
  if (class(eqSS) != "list") stop("eqSS is not of class list")

  if (!identical(
    names(eqSS),
    c(
      "age_brackets", "het_brackets", "Smat",
      "Dmat", "Amat", "Umat", "Tmat", "Pmat",
      "IBmat", "ICAmat", "ICMmat", "IDmat",
      "FOI", "phi",
      "Sv", "Ev", "Iv", "ICM_Init",
      "theta"
    )
  )) {
    stop("Incorrect variable names within equilibrium.steady.state")
  }

  dims.1 <- lapply(eqSS, function(x) {
    return(dim(x)[1])
  })
  dims.2 <- lapply(eqSS, function(x) {
    return(dim(x)[2])
  })
  if (unique(unlist(dims.1[grep("mat", names(eqSS))])) != length(eqSS$age_brackets)) {
    stop("Dimensions 1 error in equilibrium.stead.state matrices")
  }
  if (unique(unlist(dims.2[grep("mat", names(eqSS))])) != length(eqSS$het_brackets)) {
    stop("Dimensions 2 error in equilibrium.stead.state matrices")
  }

  numerics <- which(!unlist(lapply(eqSS, is.numeric)))
  if (length(numerics) != 0) stop(paste(names(numerics), "provided not numeric"))

  ## ---------------------------------------------

  # Create paramlist
  param_list <- list(
    N = N, eqSS = eqSS,
    barcode_list = barcode_list,
    spatial_list = spatial_list,
    housekeeping_list = housekeeping_list,
    drug_list = drug_list,
    nmf_list = nmf_list,
    vector_adaptation_list = vector_adaptation_list,
    core_parameter_list = mpl
  )

  return(param_list)
}

#------------------------------------------------
#' Parameter List creation for magenta simulation updating
#'
#' \code{param_list_simulation_update_create} creates suitable parameter list for
#' \code{simulation_R} for continuing a simulation from memory within the active
#' session.
#'
#' @param years Length of simulation. Default = 1
#' @param ft Treatments seeking value
#' @param mu_vec Vector of mosquito mortalities for each day with years. Default = NULL, which
#' will result in rep(0.132,floor(years*365))
#' @param fv_vec Vector of mosquito bitings for each day with years. Default = NULL, which
#' will result in rep(1/3,floor(years*365))
#' @param statePtr Pointer for current model state as return by \code{simulation_R}$Ptr
#' @param spatial_list Spatial list
#' @param drug_list Drug list
#' @param barcode_list Barcode parameter list
#'
#' @export

param_list_simulation_update_create <- function(years = 1, ft = 0.4,
                                                mu_vec = NULL, fv_vec = NULL,
                                                statePtr,
                                                spatial_list,
                                                drug_list,
                                                barcode_list) {
  ## CHECKS ##
  ## ---------------------------------------------
  if (class(statePtr) != "externalptr") stop("state.ptr is not of class externalptr")
  if (!is.numeric(years)) stop("years provided is not numeric")
  if (!is.numeric(ft)) stop("ft provided is not numeric")
  if (!is.null(mu_vec)) {
    if (!(length(mu_vec) == floor(years * 365))) {
      stop("mu_vec not long enough")
    }
  }
  if (!is.null(fv_vec)) {
    if (!(length(fv_vec) == floor(years * 365))) {
      stop("fv_vec not long enough")
    }
  }

  ## ---------------------------------------------

  # generate mu_vec if needed
  if (is.null(mu_vec)) {
    mu_vec <- rep(0.132, floor(years * 365))
  }
  # generate fv_vec if needed
  if (is.null(fv_vec)) {
    fv_vec <- rep(1 / 3, floor(years * 365))
  }

  # Create paramlist
  param_list <- list(
    years = years, ft = ft, mu_vec = mu_vec,
    fv_vec = fv_vec,
    statePtr = statePtr,
    spatial_list = spatial_list,
    drug_list = drug_list,
    barcode_list = barcode_list
  )

  return(param_list)
}

#------------------------------------------------
#' Parameter List creation for magenta simulation getting (saving to disk)
#'
#' \code{param_list_simulation_get_create} creates suitable parameter list for
#' \code{simulation_R} for continuing a simulation from memory within the active
#' session.
#'
#' @param statePtr Pointer for current model state as return by \code{simulation_R}$Ptr
#'
#' @export

param_list_simulation_get_create <- function(statePtr) {
  ## CHECKS ##
  ## ---------------------------------------------
  if (class(statePtr) != "externalptr") {
    stop("state.ptr is not of class externalptr")
  }

  ## ---------------------------------------------

  # Create paramlist
  param_list <- list(statePtr = statePtr)

  return(param_list)
}

#------------------------------------------------
#' Parameter List creation for magenta simulation finalizer
#'
#' \code{param_list_simulation_finalizer_create} creates suitable parameter list for
#' \code{simulation_R} for free memory used by a simulation
#'
#' @param statePtr Pointer for current model state as return by \code{simulation_R}$Ptr
#'
#' @export

param_list_simulation_finalizer_create <- function(statePtr) {
  ## CHECKS ##
  ## ---------------------------------------------
  if (class(statePtr) != "externalptr") {
    stop("state.ptr is not of class externalptr")
  }

  ## ---------------------------------------------

  # Create paramlist
  param_list <- list(
    statePtr = statePtr,
    finalizer = TRUE
  )

  return(param_list)
}

#------------------------------------------------
#' Parameter List creation for loading saved magenta simulation
#'
#' \code{param_list_simulation_saved_init_create} creates
#' suitable parameter list for \code{simulation_R} for continuing a
#' simulation from memory within the active session.
#'
#' @param savedState Saved state generated by
#'   \code{simulation_R} when provided with a
#'   \code{param_list_simulation_get_create} parameter list
#'
#' @export

param_list_simulation_saved_init_create <- function(savedState) {
  ## CHECKS ##
  ## ---------------------------------------------
  if (class(savedState) != "list") {
    stop("savedState is not of class list")
  }
  required <- c(
    "population_List",
    "populations_event_and_strains_List",
    "scourge_List",
    "parameters_List",
    "parameter_housekeeping_List"
  )

  if (!all(required %in% names(savedState))) {
    stop("Incorrect variable names within savedState")
  }

  # Create paramlist
  param_list <- list(savedState = savedState)

  return(param_list)
}


#------------------------------------------------
#' Deprecated dispatcher for the simulation bindings
#'
#' \code{simulation_R} used to inspect the shape of \code{param_list} to decide
#' which C++ entry point to call. It is kept for one release so that scripts
#' which call it directly keep working, and now forwards to the plainly named
#' bindings \code{sim_init()}, \code{sim_update()}, \code{sim_get_state()},
#' \code{sim_free()} and \code{sim_init_saved()}.
#'
#' @param param_list paramlist passed from
#'   \code{param_list_simulation_init_create} or from
#'   \code{param_list_simulation_update_create}
#' @param seed Seed for the simulation
#' @keywords internal
#' @export
simulation_R <- function(param_list, seed) {
  .Deprecated(
    msg = paste(
      "simulation_R() is deprecated and will be removed in a future release;",
      "pipeline() now calls the sim_init(), sim_update(), sim_get_state(),",
      "sim_free() and sim_init_saved() bindings directly."
    )
  )
  stopifnot(is.list(param_list))
  if (length(param_list) == 1 && is.list(param_list[[1]]) &&
    !is.null(param_list[[1]]$eqSS)) {
    param_list <- param_list[[1]]
  }

  if (!is.null(param_list$eqSS)) {
    return(sim_init(param_list, seed))
  }
  if (!is.null(param_list$savedState)) {
    return(sim_init_saved(param_list$savedState, seed))
  }
  if (!is.null(param_list$finalizer)) {
    return(sim_free(param_list$statePtr, seed))
  }
  if (!is.null(param_list$years) && !is.null(param_list$statePtr)) {
    return(sim_update(param_list, seed))
  }
  if (!is.null(param_list$statePtr)) {
    return(sim_get_state(param_list$statePtr, seed))
  }
  stop("No matching parameter list handler found within simulation_R")
}
