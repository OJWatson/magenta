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
#' @keywords internal

param_list_simulation_init_create <- function(N = 1e+04, eqSS, barcode_list,
                                              spatial_list, housekeeping_list,
                                              drug_list, nmf_list,
                                              vector_adaptation_list,
                                              mpl) {
  check_number(N, "N", min = 1, integer = TRUE)
  validate_eqss(eqSS)
  validate_barcode_list(barcode_list)
  validate_spatial_list(spatial_list)
  validate_housekeeping_list(housekeeping_list)
  validate_drug_list(drug_list)
  validate_nmf_list(nmf_list)
  validate_vector_adaptation_list(vector_adaptation_list)
  validate_model_param_list(mpl)

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
  stopifnot(identical(names(param_list), cpp_contract$init))
  param_list
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
#' @keywords internal

param_list_simulation_update_create <- function(years = 1, ft = 0.4,
                                                mu_vec = NULL, fv_vec = NULL,
                                                statePtr,
                                                spatial_list,
                                                drug_list,
                                                barcode_list) {
  assert_state_ptr(statePtr)
  check_number(years, "years", min = 0)
  check_number(ft, "ft", min = 0)
  n_days <- floor(years * 365)
  if (is.null(mu_vec)) {
    mu_vec <- rep(0.132, n_days)
  }
  if (is.null(fv_vec)) {
    fv_vec <- rep(1 / 3, n_days)
  }
  check_number(mu_vec, "mu_vec", min = 0, n = n_days)
  check_number(fv_vec, "fv_vec", min = 0, n = n_days)
  validate_spatial_list(spatial_list)
  validate_drug_list(drug_list)
  validate_barcode_list(barcode_list)

  param_list <- list(
    years = years, ft = ft, mu_vec = mu_vec,
    fv_vec = fv_vec,
    statePtr = statePtr,
    spatial_list = spatial_list,
    drug_list = drug_list,
    barcode_list = barcode_list
  )
  stopifnot(identical(names(param_list), cpp_contract$update))
  param_list
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
simulation_R <- function(param_list, seed) {
  .Deprecated(
    msg = paste(
      "simulation_R() is deprecated and will be removed in a future release;",
      "pipeline() now calls the sim_init(), sim_update(), sim_get_state(),",
      "sim_free() and sim_init_saved() bindings directly."
    )
  )
  stopifnot(is.list(param_list))
  # the old dispatcher seeded before every call; kept here, and only here
  set.seed(seed)
  if (length(param_list) == 1 && is.list(param_list[[1]]) &&
    !is.null(param_list[[1]]$eqSS)) {
    param_list <- param_list[[1]]
  }

  if (!is.null(param_list$eqSS)) {
    return(sim_init(param_list))
  }
  if (!is.null(param_list$savedState)) {
    return(sim_init_saved(param_list$savedState))
  }
  if (!is.null(param_list$finalizer)) {
    return(sim_free(param_list$statePtr))
  }
  if (!is.null(param_list$years) && !is.null(param_list$statePtr)) {
    return(sim_update(param_list))
  }
  if (!is.null(param_list$statePtr)) {
    return(sim_get_state(param_list$statePtr))
  }
  stop("No matching parameter list handler found within simulation_R")
}

#' Initialise a simulation from an equilibrium parameter list
#'
#' The first of the five thin bindings over the Rcpp entry points
#' (\code{sim_init()}, \code{\link{sim_update}}, \code{\link{sim_get_state}},
#' \code{\link{sim_free}} and \code{\link{sim_init_saved}}). The simulation
#' draws its random numbers from R's generator and the bindings do not seed
#' it: \code{\link{pipeline}} seeds once at the start of a run (or restores
#' the generator state a saved run stopped at) and every call continues from
#' wherever the generator is.
#'
#' @param param_list List created by \code{param_list_simulation_init_create}.
#' @return List with the state pointer (\code{Ptr}) and the \code{Loggers}.
#' @keywords internal
#' @useDynLib magenta
#' @importFrom Rcpp evalCpp
sim_init <- function(param_list) {
  Simulation_Init_cpp(param_list)
}

#' Continue a simulation held in memory
#'
#' @param param_list List created by \code{param_list_simulation_update_create}.
#' @return List with the state pointer (\code{Ptr}) and the \code{Loggers}.
#' @keywords internal
sim_update <- function(param_list) {
  Simulation_Update_cpp(param_list)
}

#' Copy the full model state back to R
#'
#' @param ptr External pointer to the simulation state.
#' @return List of the population, strain, mosquito and parameter lists.
#' @keywords internal
sim_get_state <- function(ptr) {
  assert_state_ptr(ptr)
  Simulation_Get_cpp(list(statePtr = ptr))
}

#' Free the memory held by a simulation
#'
#' @param ptr External pointer to the simulation state.
#' @return Invisibly, the list returned by the finaliser.
#' @keywords internal
sim_free <- function(ptr) {
  assert_state_ptr(ptr)
  invisible(Simulation_Finalizer_cpp(list(statePtr = ptr, finalizer = TRUE)))
}

#' Initialise a simulation from a saved model state
#'
#' @param saved_state List as returned by \code{\link{sim_get_state}}, without
#'   the R-side \code{saved_inputs} element.
#' @return List with the state pointer (\code{Ptr}) and the \code{Loggers}.
#' @keywords internal
sim_init_saved <- function(saved_state) {
  assert_saved_state(saved_state)
  Simulation_Saved_Init_cpp(list(savedState = saved_state))
}

#' Check for a simulation state pointer
#'
#' @param ptr The value.
#' @return \code{ptr}, invisibly, if it is an external pointer.
#' @keywords internal
assert_state_ptr <- function(ptr) {
  if (!inherits(ptr, "externalptr")) {
    stop("statePtr must be an external pointer to a simulation state")
  }
  invisible(ptr)
}

#' Check a saved state before handing it to C++
#'
#' @param saved_state The value.
#' @return \code{saved_state}, invisibly, if it is a list with the five lists
#'   \code{\link{sim_get_state}} returns.
#' @keywords internal
assert_saved_state <- function(saved_state) {
  required <- c(
    "population_List",
    "populations_event_and_strains_List",
    "scourge_List",
    "parameters_List",
    "parameter_housekeeping_List"
  )
  if (!is.list(saved_state)) {
    stop("saved_state must be a list as returned by sim_get_state()")
  }
  missing <- setdiff(required, names(saved_state))
  if (length(missing)) {
    stop(
      "saved_state is missing the elements: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(saved_state)
}
