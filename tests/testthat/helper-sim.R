# Helpers for driving the simulation bindings directly in unit tests.

quiet_housekeeping_test <- function() {
  housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)
}

# A tiny initialised simulation state (two loci, no interventions) built the
# way pipeline() builds one.
tiny_state <- function(seed = 1L, N = 50, years = 1, update_length = 365,
                       drug_list = drug_list_create(), spatial_type = NULL,
                       spatial_incidence = NULL, spatial_mosquitoFOI = NULL,
                       plaf = rep(0.5, 2), ibd_length = 1, saved_state = NULL) {
  schedule <- build_schedule(
    years = years, update_length = update_length, ft = 0.4, plaf = plaf,
    resistance_flag = drug_list$resistance_flag,
    spatial_incidence = spatial_incidence, spatial_mosquitoFOI = spatial_mosquitoFOI
  )
  # pipeline() seeds once before set up; the importation vectors drawn during
  # set up depend on it
  set.seed(seed)
  magenta:::sim_setup(
    N = N, EIR = 1, country = NULL, admin = NULL,
    mpl = model_param_list_create(), schedule = schedule,
    barcode = list(
      num_loci = length(plaf), ibd_length = ibd_length,
      prob_crossover = rep(0.5, length(plaf)),
      starting_ibd = 0, mutation_rate = rep(1e-7, length(plaf)),
      mutation_treated_modifier = 1
    ),
    spatial = list(
      spatial_type = spatial_type, island_imports_plaf_linked_flag = FALSE,
      survival_percentage = 0.2, oocyst_mean = 2.5, oocyst_shape = 1
    ),
    housekeeping_list = quiet_housekeeping_test(), drug_list = drug_list,
    nmf_list = nmf_list_create(),
    vector_adaptation_list = vector_adaptation_list_create(length(plaf)),
    saved_state = saved_state
  )
}

# Constant mosquito death and biting rates for n days.
constant_mu_fv <- function(n) {
  data.frame(mu = rep(0.132, n), fv = rep(0.333, n))
}
