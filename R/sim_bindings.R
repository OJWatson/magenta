# Thin R bindings over the five Rcpp entry points.
#
# The simulation draws its random numbers from R's generator. The bindings do
# not seed it: pipeline() seeds once at the start of a run (or restores the
# generator state a saved run stopped at) and every call continues from
# wherever the generator is.

#' Initialise a simulation from an equilibrium parameter list
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

assert_state_ptr <- function(ptr) {
  if (!inherits(ptr, "externalptr")) {
    stop("statePtr must be an external pointer to a simulation state")
  }
  invisible(ptr)
}

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
