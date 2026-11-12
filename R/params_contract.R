# The contract between the R parameter lists and the C++ entry points.
#
# `cpp_contract` records, for every list R hands to C++, the field names in
# the order the R constructor produces them (`fields`) and the subset that the
# C++ side reads (`cpp`; see src/main_init.cpp, src/main_update.cpp and
# src/main_saved_init.cpp). The validators below check a list against it, and
# tests/testthat/test-params.R checks that each constructor produces exactly
# the `fields` set, so a change on either side fails a test rather than
# failing inside C++.

#' Field names read by the C++ simulation
#'
#' \code{cpp_contract} lists, for every parameter list passed from R to the
#' C++ simulation, the fields the R constructor produces (\code{fields}, in
#' order) and the subset the C++ side reads (\code{cpp}). The
#' \code{validate_*()} functions check a list against it and the package
#' tests check that every constructor produces exactly the recorded set, so a
#' change on either side is caught by a test rather than inside C++.
#'
#' @format A named list. The entries \code{init}, \code{update},
#'   \code{saved_init}, \code{get}, \code{finalizer} and \code{saved_state}
#'   are character vectors naming the top-level fields of each entry point;
#'   \code{housekeeping}, \code{barcode}, \code{drug}, \code{drug_entry},
#'   \code{vector_adaptation}, \code{nmf}, \code{spatial}, \code{eqss} and
#'   \code{core} are lists with \code{fields} and \code{cpp}.
#' @keywords internal
#' @name cpp_contract
NULL

cpp_contract <- list(
  # top-level lists for each entry point
  init = c(
    "N", "eqSS", "barcode_list", "spatial_list", "housekeeping_list",
    "drug_list", "nmf_list", "vector_adaptation_list", "core_parameter_list"
  ),
  update = c(
    "years", "ft", "mu_vec", "fv_vec", "statePtr",
    "spatial_list", "drug_list", "barcode_list"
  ),
  saved_init = "savedState",
  get = "statePtr",
  finalizer = c("statePtr", "finalizer"),
  saved_state = c(
    "population_List", "populations_event_and_strains_List", "scourge_List",
    "parameters_List", "parameter_housekeeping_List"
  ),
  # the parameter lists themselves
  housekeeping = list(
    fields = c("quiet_print", "quiet_test_print", "cluster", "clear_up"),
    cpp = c("quiet_print", "quiet_test_print")
  ),
  barcode = list(
    fields = c(
      "num_loci", "ibd_length", "starting_ibd", "plaf", "prob_crossover",
      "barcode_type", "mutation_flag", "mutation_rate",
      "mutation_treated_modifier"
    ),
    cpp = c(
      "num_loci", "ibd_length", "starting_ibd", "plaf", "prob_crossover",
      "barcode_type", "mutation_flag", "mutation_rate",
      "mutation_treated_modifier"
    )
  ),
  drug = list(
    fields = c(
      "resistance_flag", "number_of_resistance_loci", "resistance_loci",
      "artemisinin_loci", "cost_of_resistance", "absolute_fitness_cost_flag",
      "drugs", "mft_flag", "temporal_cycling", "next_temporal_cycle",
      "sequential_cycling", "sequential_update", "number_of_drugs",
      "drug_choice", "partner_drug_ratios"
    ),
    cpp = c(
      "resistance_flag", "number_of_resistance_loci", "resistance_loci",
      "artemisinin_loci", "cost_of_resistance", "absolute_fitness_cost_flag",
      "drugs", "mft_flag", "number_of_drugs", "drug_choice",
      "partner_drug_ratios"
    )
  ),
  drug_entry = list(
    fields = c(
      "lpf", "barcode_positions", "prophylactic_positions", "dur_P", "dur_SPC",
      "hill_n", "hill_kA", "hill_res_n", "hill_res_kA"
    ),
    cpp = c(
      "lpf", "barcode_positions", "prophylactic_positions", "dur_P", "dur_SPC",
      "hill_n", "hill_kA", "hill_res_n", "hill_res_kA"
    )
  ),
  vector_adaptation = list(
    fields = c(
      "vector_adaptation_flag", "vector_adaptation_loci",
      "local_oocyst_advantage", "gametocyte_sterilisation_flag",
      "gametocyte_sterilisation", "oocyst_reduction_by_artemisinin"
    ),
    cpp = c(
      "vector_adaptation_flag", "vector_adaptation_loci",
      "local_oocyst_advantage", "gametocyte_sterilisation_flag",
      "gametocyte_sterilisation", "oocyst_reduction_by_artemisinin"
    )
  ),
  nmf = list(
    fields = c(
      "nmf_flag", "mean_nmf_frequency", "nmf_age_brackets",
      "prob_of_testing_nmf"
    ),
    cpp = c(
      "nmf_flag", "mean_nmf_frequency", "nmf_age_brackets",
      "prob_of_testing_nmf"
    )
  ),
  spatial = list(
    fields = c(
      "spatial_type", "cotransmission_freq_vector", "oocyst_freq_vector",
      "imported_cotransmissions_events", "imported_oocyst_events", "plaf",
      "island_imports_plaf_linked_flag"
    ),
    cpp = c(
      "spatial_type", "cotransmission_freq_vector", "oocyst_freq_vector",
      "imported_cotransmissions_events", "imported_oocyst_events", "plaf",
      "island_imports_plaf_linked_flag"
    )
  ),
  eqss = list(
    fields = c(
      "age_brackets", "het_brackets", "Smat", "Dmat", "Amat", "Umat", "Tmat",
      "Pmat", "IBmat", "ICAmat", "ICMmat", "IDmat", "FOI", "phi", "Sv", "Ev",
      "Iv", "ICM_Init", "theta"
    ),
    cpp = c(
      "age_brackets", "het_brackets", "Smat", "Dmat", "Amat", "Umat", "Tmat",
      "Pmat", "IBmat", "ICAmat", "ICMmat", "IDmat", "FOI", "phi", "Sv", "Ev",
      "Iv", "ICM_Init", "theta"
    )
  ),
  # the core (epidemiological) parameter list. C++ reads phi0 only; the rest
  # is consumed by equilibrium_init_create() and the odin model.
  core = list(
    fields = c(
      "DY", "eta", "rho", "a0", "sigma2", "max_age",
      "rA", "rT", "rD", "rU", "rP", "dE", "delayGam",
      "cD", "cT", "cU", "gamma1",
      "d1", "dID", "ID0", "kD", "uD", "aD", "fD0", "gammaD",
      "alphaA", "alphaU",
      "b0", "b1", "dB", "IB0", "kB", "uB",
      "phi0", "phi1", "dCA", "IC0", "kC", "uCA", "PM", "dCM",
      "delayMos", "tau1", "tau2", "mu0", "Q0", "chi", "bites_Bed",
      "bites_Indoors", "fv0", "av0", "Surv0", "p10", "p2",
      "muEL", "muLL", "muPL", "dEL", "dLL", "dPL", "gammaL", "km", "cm",
      "betaL", "eov", "b_lambda", "lambda",
      "num_int", "itn_cov", "irs_cov", "ITN_IRS_on", "pop_split", "cov",
      "d_ITN0", "r_ITN0", "r_ITN1", "r_IRS0", "d_IRS0",
      "irs_half_life", "itn_half_life", "IRS_interval", "ITN_interval",
      "irs_loss", "itn_loss"
    ),
    cpp = "phi0"
  )
)

# ---- generic checks ---------------------------------------------------------
#
# Each check reports the list and field it is looking at, e.g.
# "drug_list$drug_choice must be a single number", and returns the value
# invisibly so the validators read as a list of requirements.

check_fields <- function(x, contract, what) {
  if (!is.list(x)) {
    stop(sprintf("%s must be a list", what), call. = FALSE)
  }
  missing <- setdiff(contract$fields, names(x))
  extra <- setdiff(names(x), contract$fields)
  if (length(missing)) {
    stop(sprintf(
      "%s is missing the field(s): %s", what, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (length(extra)) {
    stop(sprintf(
      "%s has unexpected field(s): %s", what, paste(extra, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(x)
}

check_flag <- function(x, what, n = 1L) {
  if (!is.logical(x) || anyNA(x) || (!is.null(n) && length(x) != n)) {
    stop(sprintf(
      "%s must be %s", what,
      if (identical(n, 1L)) "TRUE or FALSE" else "a logical vector without NAs"
    ), call. = FALSE)
  }
  invisible(x)
}

check_number <- function(x, what, min = -Inf, max = Inf, n = 1L,
                         integer = FALSE) {
  kind <- if (integer) "whole number" else "number"
  if (!(is.numeric(x) || (integer && is.logical(x))) || anyNA(x)) {
    stop(sprintf(
      "%s must be %s", what,
      if (identical(n, 1L)) paste0("a single ", kind) else paste0("a ", kind, " vector")
    ), call. = FALSE)
  }
  if (!is.null(n)) {
    check_length(x, what, n)
  }
  in_range <- all(x >= min) && all(x <= max) && (!integer || all(x == round(x)))
  if (!in_range) {
    range <- if (is.finite(min) && is.finite(max)) {
      sprintf(" between %s and %s", format(min), format(max))
    } else if (is.finite(min)) {
      sprintf(" of at least %s", format(min))
    } else if (is.finite(max)) {
      sprintf(" of at most %s", format(max))
    } else {
      ""
    }
    stop(sprintf(
      "%s must be %s%s", what,
      if (identical(n, 1L)) paste0("a single ", kind) else paste0("a ", kind, " vector"),
      range
    ), call. = FALSE)
  }
  invisible(x)
}

check_length <- function(x, what, n) {
  if (length(x) != n) {
    stop(sprintf("%s must have length %d, not %d", what, n, length(x)),
      call. = FALSE
    )
  }
  invisible(x)
}

check_class <- function(x, class, what) {
  if (!inherits(x, class)) {
    stop(sprintf(
      "%s must be created by %s()", what,
      sub("^magenta_", "", class) |> paste0("_create")
    ), call. = FALSE)
  }
  invisible(x)
}

# Stamp a parameter list with its class, keeping it a plain list underneath so
# that the C++ side and existing code that treats it as a list are unaffected.
new_param_list <- function(x, class) {
  structure(x, class = c(class, "magenta_param_list", "list"))
}

#' @export
print.magenta_param_list <- function(x, ...) {
  cat("<", class(x)[1], "> with ", length(x), " fields\n", sep = "")
  utils::str(unclass(x), max.level = 1, give.attr = FALSE, no.list = TRUE)
  invisible(x)
}
