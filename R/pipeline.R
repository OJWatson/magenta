#' Pipeline for cluster submission
#'
#' \code{pipeline} steps through creating the parameter list, the equilibrium
#' initialisation and steady state creation before checking and passing suitable
#' parameters to the simulation. This is then saved. If a path to a
#' savedState is provided then this state is loaded and continued.
#'
#' # Main Params
#'
#' @param N Population Size. Default = 100000
#' @param ft Vector of treatment frequency. Default = 0.4
#' @param years Lenth of simulation. Default = 20
#' @param EIR Numeric for desired annual EIR. Default = 120
#' @param country Character for country within which admin2 is in.
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be
#'   used to match. If not provided then no seasonality is introduced.
#'   Default = NULL
#' @param itn_cov Vector for ITN coverages that change at update_length intervals.
#'   Default = 0
#' @param irs_cov Vector for IRS coverages that change at update_length intervals.
#'   Default = 0
#' @param use_historic_interventions Boolean as to whether to use interventions
#'   on file for the admin and country specified. If TRUE then provide the years
#'   as a year range, e.g. 2000:2015. WARNING - Best to have this as FALSE and
#'   manually specify the itn_cov, irs_cov and ft.
#' @param survival_percentage Mumeric for % of sporozoites surviving.
#'   Default = 0.2
#' @param oocyst_mean Mean for number of oocysts formed from a bite. Default=2.5
#' @param oocyst_shape Shape parameter for oocysts formed. Default=1
#'
#' # Spatial
#'
#' @param spatial_type Default = NULL. If spatial is wanted then provide
#'   a character describing the type of spatial simulation, which must be
#'   one of "island" or "metapop".
# @param redis_host Default = "fi--dideclusthn.dide.ic.ac.uk". This should
#   be the host address where your redis server is running.
# @param spatial_uuid Default = NULL. If spatial is provided, then this will
#   error unless a character string is passed to this argument. This character
#   should be the same for each parallel task within the same job.
#' @param spatial_incidence_matrix Spatial incidence for humans,
#' i.e. importation vector
#' @param spatial_mosquitoFOI_matrix Spatial mosquio FOI, i.e. importation
#'  to mosquitoes vector
#'
#' # Genetic Params
#'
#' @param num_loci Number of loci. Default = 24
#' @param ibd_length If we are simulating IBD dynamics, each loci is now
#'   represented by a bitset of ibd_length. Thus ibd_length needs to be long
#'   enough to ensure that as new identity relationships occur, i.e.
#'   an importation barcode will be a new identity. e.g. If are population is
#'   1000, we may expect at 80% prevalence, with a mean COI of 3 we will need
#'   2400 different identities, i.e. 2^ibd_length > 2400. However, keep in mind
#'   importations as these need to be continually new, i.e. if we are simulating
#'   for 30 years, with 3 importations a day, then we will need at least length
#'   to ensure that 2^ibd_length > 2400 + (30*365*3). This will probably be
#'   automatically calculated in the future. If we are not ding IBD, then this
#'   should be 1, which is the defalt.
#' @param plaf Vector of population level allele frequencies for the barcode.
#'   Default = rep(0.5, num_loci)
#' @param prob_crossover Vector of probabilities for crossover events for the
#'   barcode. Default = rep(0.5, num_loci)
#' @param starting_ibd Starting IBD. Default = 0, which means that each infected
#'   individual at initialisation is given a unique ID for their parasites.
#' @inheritParams barcode_list_create
#' @param mutation_rate Probability of mutation occuring and fixing
#' @param mutation_flag Boolean for simulating mutations
#'
#' # Saving Params
#'
#' @param full_save Boolean detailing whether the entire simulation is saved.
#'   Default = FALSE
#' @param full_update_save Boolean to save entire simualation at each update
#'   save. Default = FALSE
#' @param human_only_full_save Boolean detailing whether just the human
#' component of the simulation is saved within full_save. Default = FALSE
#' @param update_save Boolean detailing whether the logging output is saved
#' each update_length up to years. Default = FALSE
#' @param update_length How long each update is run for in days.
#'   Default = 365
#' @param human_update_save Boolean detailing if the human state is also
#' saved during each update_length. Default = FALSE
#' @param summary_saves_only Boolean if summary tables about COI are saved
#' within human yearly save only. Dataframes of age, clinical status
#' binned COI.
#' @param mean_only Boolean for returning only the mean when summarising the
#'   population COI, COU etc. Default = TRUE
#' @param save_lineages Boolean for whether we save the frequency of each strain
#'   when summarising with \code{genetics_df_without_summarising}=TRUE.
#'   Default = FALSE
#' @param only_allele_freqs Boolean for returning the summarised genetics (allele
#'   frequencies and maybe strain frequencies) or the whole data frame produced
#'   by \code{\link{pop_strains_df}}. Default = TRUE
#' @param saved_state_path Full file path to a saved model state to be
#'   loaded and continued. Default = NULL,
#'   which will trigger initialisation
#' @param genetics_df_without_summarising Boolean for returning just the
#'   genetics data frame without summarising with \code{\link{COI_df_create}}.
#'   Default = FALSE
#' @param update_save_func As opposed to having to provide arguments for the
#'   update behaviour, you can pass in a function. See \code{update_saves}
#'   for the default one.
#' @param set_up_only Boolean for whether to return just the initialised simulation.
#'   Default = FALSE
#' @param sample_size Numeric for number of individuals to be sampled at the end
#'   of each update. Default = Inf, which samples everyone. If you provide a
#'   vector of sample sizes it will sample at each specified sample size.
#' @param sample_states Numeric for which sample infection states are to be
#'   included in sampling. Default = 0:5 (i.e. all states). 1:4 for example
#'   would ensure only infected individuals are included.
#' @param age_breaks What age breaks are used when summarising the population.
#'   Default is `c(-0.001, 5, 15, 100.1)`
#' @param sample_reps Numeric for how many sample reps are done. Default = 1.
#'
#' # Parameter Lists
#'
#' @param housekeeping_list List created by  \code{\link{housekeeping_list_create}}
#' @param drug_list List created by \code{\link{drug_list_create}}
#' @param vector_adaptation_list List created by \code{\link{vector_adaptation_list_create}}
#' @param nmf_list List created by \code{\link{nmf_list_create}}
#' @param seed Random seed. Default is Random
#' @param ... Other parameters to model_param_list_create
#'
#' \code{pipeline}
#'
#' @export
#' @inheritParams spl_create

pipeline <- function(EIR = 120,
                     ft = 0.4,
                     itn_cov = 0,
                     irs_cov = 0,
                     use_historic_interventions = FALSE,
                     survival_percentage = 0.20,
                     oocyst_mean = 2.5,
                     oocyst_shape = 1,
                     N = 100000,
                     years = 20,
                     update_length = 365,
                     country = NULL,
                     admin = NULL,
                     spatial_type = NULL,
                     # redis_host = "fi--dideclusthn.dide.ic.ac.uk",
                     # spatial_uuid = NULL,
                     spatial_incidence_matrix = NULL,
                     spatial_mosquitoFOI_matrix = NULL,
                     island_imports_plaf_linked_flag = FALSE,
                     num_loci = 24,
                     ibd_length = 1,
                     plaf = rep(0.5, num_loci),
                     prob_crossover = rep(0.5, num_loci),
                     starting_ibd = 0.0,
                     mutation_rate = rep(1e-7, num_loci),
                     mutation_flag = FALSE,
                     mutation_treated_modifier = 1,
                     full_save = FALSE,
                     full_update_save = FALSE,
                     human_only_full_save = FALSE,
                     update_save = FALSE,
                     update_save_func = NULL,
                     human_update_save = FALSE,
                     genetics_df_without_summarising = FALSE,
                     summary_saves_only = FALSE,
                     set_up_only = FALSE,
                     mean_only = TRUE,
                     save_lineages = FALSE,
                     saved_state_path = NULL,
                     seed = as.integer(runif(1, 1, 1000000000)),
                     sample_size = Inf,
                     sample_states = 0:5,
                     age_breaks = c(-0.001, 5, 15, 100.1),
                     sample_reps = 1,
                     housekeeping_list = housekeeping_list_create(),
                     drug_list = drug_list_create(),
                     vector_adaptation_list = vector_adaptation_list_create(num_loci),
                     only_allele_freqs = TRUE,
                     nmf_list = nmf_list_create(),
                     ...) {
  # grab the function call for reproducibility
  args <- append(as.list(environment()), list(...))
  call <- match.call(expand.dots = TRUE)
  call$seed <- args$seed

  # PRE-SET UP HOUSEKEEPING ------------------------ ####
  set.seed(seed)
  message(paste0("Seed set to ", seed))
  message("magenta v", utils::packageVersion("magenta"))

  # historic intervention grab if asked
  if (use_historic_interventions && length(years) > 1) {
    ints <- intervention_grab(country = country, admin = admin, year_range = years)
    spl <- spl_grab(country = country, admin = admin, year_range = years)
    spatial_incidence_matrix <- spl$incidence
    spatial_mosquitoFOI_matrix <- spl$mosquitoFOI
    itn_cov <- ints$itn_cov
    irs_cov <- ints$irs_cov
    ft <- ints$ft
  }

  # INITIALISATION --------------------------------- ####
  inputs <- list(
    EIR = EIR, ft = ft, years = years, country = country, admin = admin,
    num_loci = num_loci, ibd_length = ibd_length, plaf = plaf,
    prob_crossover = prob_crossover, starting_ibd = starting_ibd,
    mutation_flag = mutation_flag, mutation_rate = mutation_rate,
    mutation_treated_modifier = mutation_treated_modifier,
    spatial_type = spatial_type,
    spatial_incidence_matrix = spatial_incidence_matrix,
    spatial_mosquitoFOI_matrix = spatial_mosquitoFOI_matrix,
    island_imports_plaf_linked_flag = island_imports_plaf_linked_flag,
    survival_percentage = survival_percentage,
    oocyst_mean = oocyst_mean, oocyst_shape = oocyst_shape
  )

  if (is.null(saved_state_path)) {
    state <- setup_fresh(
      N = N, inputs = inputs, mpl = model_param_list_create(...),
      housekeeping_list = housekeeping_list, drug_list = drug_list,
      nmf_list = nmf_list, vector_adaptation_list = vector_adaptation_list,
      seed = seed
    )
  } else {
    state <- setup_from_saved(saved_state_path, inputs = inputs, seed = seed)
  }
  housekeeping_list <- state$housekeeping_list

  # INITIALISATION ONLY ---------------------------- ####
  if (set_up_only) {
    message("Set Up Grab")
    sim_save <- prepare_state_save(state, seed)
    if (housekeeping_list$clear_up) {
      sim_free(state$sim$Ptr, seed)
      gc()
    }
    return(collect_output(sim_save, human_only = human_only_full_save && !full_save))
  }

  # DETERMINISTIC MOSQUITO PRE-SIMULATION ---------- ####
  mu_fv <- mu_fv_create(
    eqInit = state$eqInit, ft = state$ft_odin, itn_cov = itn_cov,
    irs_cov = irs_cov, years = years
  )

  # SIMULATION ------------------------------------- ####
  if (update_save) {
    if (length(ft) == 1) {
      ft <- rep(ft, ceiling(years))
    }
    saver <- function(res, i, sim_out) {
      if (is.null(update_save_func)) {
        update_saves(
          res = res, i = i, sim.out = sim_out,
          sample_states = sample_states, sample_size = sample_size,
          sample_reps = sample_reps, mean_only = mean_only,
          age_breaks = age_breaks, barcode_list = state$barcode_list,
          num_loci = num_loci, full_update_save = full_update_save,
          genetics_df_without_summarising = genetics_df_without_summarising,
          save_lineages = save_lineages, human_update_save = human_update_save,
          summary_saves_only = summary_saves_only,
          only_allele_freqs = only_allele_freqs, mpl = state$mpl, seed = seed
        )
      } else {
        update_save_func(res, i, sim_out, state$mpl, num_loci)
      }
    }
    run <- run_update_loop(state, mu_fv, years, update_length, ft, saver, seed)
    state <- run$state
    res <- run$res
    if (full_save || human_only_full_save) {
      sim_save <- prepare_state_save(state, seed)
      res[[length(res)]] <- collect_output(sim_save, human_only = human_only_full_save)
    } else {
      res[[length(res)]] <- state$sim
    }
  } else {
    state <- run_single_period(state, mu_fv, years, ft, seed)
    if (full_save || human_only_full_save) {
      sim_save <- prepare_state_save(state, seed)
      res <- collect_output(sim_save,
        human_only = human_only_full_save && !full_save,
        human_vars = c("Infection_States", "Zetas", "Ages")
      )
    } else {
      res <- state$sim
    }
  }

  # FINISH ----------------------------------------- ####
  meta <- list(
    call = call, seed = seed, Ptr = state$sim$Ptr,
    version = utils::packageVersion("magenta")
  )
  if (update_save) {
    meta$times <- run$comp_times
  }
  attr(res, "meta") <- meta

  if (housekeeping_list$clear_up) {
    res <- return_sim_memory(res)
  }
  res
}

# Build the equilibrium, barcode and spatial inputs from the pipeline()
# arguments. This is needed whether the simulation starts fresh or from a saved
# state, because the deterministic mosquito model and the per-period spatial
# lists are always rebuilt from the current arguments.
#
# The order of the random draws matters: equilibrium_init_create() is followed
# by set.seed(seed) and the importation frequency vectors are then drawn by
# spl_create() from that freshly seeded generator.
prepare_inputs <- function(EIR, ft, years, country, admin, mpl, drug_list, seed,
                           num_loci, ibd_length, plaf, prob_crossover,
                           starting_ibd, mutation_flag, mutation_rate,
                           mutation_treated_modifier, spatial_type,
                           spatial_incidence_matrix, spatial_mosquitoFOI_matrix,
                           island_imports_plaf_linked_flag,
                           survival_percentage, oocyst_mean, oocyst_shape) {
  # geometric age brackets
  age_vector <- age_brackets(100, 40, TRUE)

  # change the ft for the equilibrium and odin to reflect 28 day failure rates
  lpfs <- unlist(lapply(drug_list$drugs, function(x) x$lpf[1]))
  ft_odin <- ft * weighted.mean(lpfs, drug_list$partner_drug_ratios)

  # near equilibrium initial condition
  eqInit <- equilibrium_init_create(
    age_vector = age_vector,
    het_brackets = 5,
    ft = ft_odin[1],
    EIR = EIR,
    country = country,
    admin = admin,
    model_param_list = mpl
  )

  # reset seed here as there is some randomness in equilibrium
  set.seed(seed)
  eqSS <- equilibrium_ss_create(eqInit = eqInit)

  # barcode parameters
  plaf_matrix <- plaf_matrix_check(plaf, years)
  barcode_list <- barcode_list_create(
    num_loci = num_loci,
    ibd_length = ibd_length,
    plaf = plaf_matrix[1, ],
    prob_crossover = prob_crossover,
    starting_ibd = starting_ibd,
    mutation_flag = mutation_flag,
    mutation_rate = mutation_rate,
    mutation_treated_modifier = mutation_treated_modifier
  )
  if (length(mutation_flag) == 1) {
    mutation_flag <- rep(mutation_flag, ceiling(years))
  }

  # spatial parameters
  spatial_type <- spatial_type_code(spatial_type)
  spatial_incidence_matrix <- spl_matrix_check(spatial_incidence_matrix, years)
  spatial_mosquitoFOI_matrix <- spl_matrix_check(spatial_mosquitoFOI_matrix, years)
  spatial_list <- spl_create(
    spatial_type = spatial_type,
    human_importation_rate_vector = spatial_incidence_matrix[1, ],
    mosquito_imporation_rate_vector = spatial_mosquitoFOI_matrix[1, ],
    cotransmission_freq_vector = ztrgeomintp(10000, 10, survival_percentage),
    oocyst_freq_vector = ztrnbinom(10000, mean = oocyst_mean, size = oocyst_shape),
    plaf = plaf_matrix[1, ],
    island_imports_plaf_linked_flag = island_imports_plaf_linked_flag
  )

  # drug parameters
  resistance_flags <- drug_list$resistance_flag
  if (length(resistance_flags) == 1) {
    resistance_flags <- rep(resistance_flags, ceiling(years))
  }
  drug_list$resistance_flag <- resistance_flags[1]

  list(
    eqInit = eqInit, eqSS = eqSS, ft_odin = ft_odin,
    plaf_matrix = plaf_matrix, barcode_list = barcode_list,
    mutation_flags = mutation_flag,
    spatial_type = spatial_type,
    spatial_incidence_matrix = spatial_incidence_matrix,
    spatial_mosquitoFOI_matrix = spatial_mosquitoFOI_matrix,
    spatial_list = spatial_list,
    resistance_flags = resistance_flags, drug_list = drug_list,
    survival_percentage = survival_percentage,
    oocyst_mean = oocyst_mean, oocyst_shape = oocyst_shape, mpl = mpl
  )
}

# Map the user-facing spatial_type to the integer code the C++ reads.
spatial_type_code <- function(spatial_type) {
  if (is.null(spatial_type)) {
    return(0)
  }
  if (spatial_type == "metapop") {
    stop("metapopulation simulations are not implemented")
  }
  if (spatial_type == "island") {
    return(1)
  }
  spatial_type
}

# Start a new simulation from the equilibrium solution. Returns the simulation
# "state": the C++ handle plus every R-side list the main loop and the savers
# need.
setup_fresh <- function(N, inputs, mpl, housekeeping_list, drug_list, nmf_list,
                        vector_adaptation_list, seed) {
  state <- do.call(
    prepare_inputs,
    c(inputs, list(mpl = mpl, drug_list = drug_list, seed = seed))
  )
  state$housekeeping_list <- housekeeping_list
  state$nmf_list <- nmf_list
  state$vector_adaptation_list <- vector_adaptation_list

  pl <- param_list_simulation_init_create(
    N = N, eqSS = state$eqSS,
    barcode_list = state$barcode_list,
    spatial_list = state$spatial_list,
    housekeeping_list = housekeeping_list,
    drug_list = state$drug_list,
    nmf_list = nmf_list,
    vector_adaptation_list = vector_adaptation_list,
    mpl = mpl
  )
  state$sim <- sim_init(pl, seed)
  state
}

# Continue a simulation from a saved state on disk. The housekeeping, drug,
# vector adaptation, non-malarial fever and core parameter lists are taken from
# the file; the equilibrium, barcode and spatial inputs are rebuilt from the
# current arguments exactly as for a fresh start.
setup_from_saved <- function(saved_state_path, inputs, seed) {
  saved_state <- readRDS(saved_state_path)
  saved_inputs <- saved_state$saved_inputs
  if (is.null(saved_inputs)) {
    stop("Saved state missing saved_inputs; please re-run the simulation with full_save = TRUE")
  }
  saved_state$saved_inputs <- NULL

  state <- do.call(
    prepare_inputs,
    c(inputs, list(
      mpl = saved_inputs$core_parameter_list,
      drug_list = saved_inputs$drug_list,
      seed = seed
    ))
  )
  state$housekeeping_list <- saved_inputs$housekeeping_list
  state$nmf_list <- saved_inputs$nmf_list
  state$vector_adaptation_list <- saved_inputs$vector_adaptation_list

  state$sim <- sim_init_saved(saved_state, seed)
  rm(saved_state)
  gc()
  state
}

# Rebuild the spatial list for a new simulation year. The importation
# frequency vectors are redrawn each time, from the generator state the
# previous C++ call left behind.
annual_spatial_list <- function(state, year) {
  spl_create(
    spatial_type = state$spatial_type,
    human_importation_rate_vector = state$spatial_incidence_matrix[year, ],
    mosquito_imporation_rate_vector = state$spatial_mosquitoFOI_matrix[year, ],
    cotransmission_freq_vector = ztrgeomintp(10000, 10, state$survival_percentage),
    oocyst_freq_vector = ztrnbinom(10000,
      mean = state$oocyst_mean, size = state$oocyst_shape
    ),
    plaf = state$plaf_matrix[year, ]
  )
}

# Run the simulation forward over one block of days.
run_period <- function(state, days, ft, mu_fv, seed) {
  pl <- param_list_simulation_update_create(
    years = length(days) / 365,
    ft = ft,
    mu_vec = mu_fv$mu[days],
    fv_vec = mu_fv$fv[days],
    spatial_list = state$spatial_list,
    drug_list = state$drug_list,
    barcode_list = state$barcode_list,
    statePtr = state$sim$Ptr
  )
  state$sim <- sim_update(pl, seed)
  state
}

# Run the whole simulation as a single period with no intermediate saves.
run_single_period <- function(state, mu_fv, years, ft, seed) {
  pl <- param_list_simulation_update_create(
    years = years, ft = min(ft[ft > 0]),
    mu_vec = mu_fv$mu,
    fv_vec = mu_fv$fv,
    spatial_list = state$spatial_list,
    drug_list = state$drug_list,
    barcode_list = state$barcode_list,
    statePtr = state$sim$Ptr
  )
  state$sim <- sim_update(pl, seed)
  state
}

# Step through the simulation in blocks of update_length days, calling the
# saver after every block but the last. Annual inputs (ft, resistance and
# mutation flags, importation and plaf rows) are switched when a block starts
# in a new year. The last block is run without a save so that the caller can
# decide what to put in the final slot.
run_update_loop <- function(state, mu_fv, years, update_length, ft, saver, seed) {
  n_periods <- round((years * 365) / update_length)
  res <- vector("list", n_periods)
  update_times <- update_length * seq_len(n_periods)
  comp_times <- rep(0, n_periods - 1)

  year <- 1
  ft_now <- ft[year]

  message("Starting Stochastic Simulation for ", years, " years")
  p <- progress::progress_bar$new(
    format = paste0(" Running: [:bar] :percent eta: :eta"),
    total = n_periods
  )
  p_print <- progress_logging(state$housekeeping_list, res, p, initial = TRUE)

  if (n_periods > 1) {
    for (i in seq_len(n_periods - 1)) {
      p_print <- progress_logging(state$housekeeping_list, res, p, i,
        initial = FALSE, p_print = p_print
      )
      comp_times[i] <- Sys.time()

      # annual updates
      if (floor((update_times[i] - update_length + 1) / 365) == year) {
        year <- year + 1
        ft_now <- ft[year]
        state$drug_list$resistance_flag <- state$resistance_flags[year]
        state$barcode_list$mutation_flag <- state$mutation_flags[year]
        state$spatial_list <- annual_spatial_list(state, year)
      }

      days <- seq_len(update_length) + ((i - 1) * update_length)
      state <- run_period(state, days, ft_now, mu_fv, seed)
      res <- saver(res, i, state$sim)

      # drug resistance updates
      state$drug_list <- drug_list_update(
        state$drug_list, year, recent_treatment_failure(res, i)
      )
    }
  }

  # the last period may not occupy a full update length
  days <- seq_len(update_length) + ((n_periods - 1) * update_length)
  days <- days[!is.na(mu_fv$mu[days])]
  state <- run_period(state, days, ft_now, mu_fv, seed)

  list(state = state, res = res, comp_times = comp_times)
}

# Treatment failure over the last two saved periods, used by the sequential
# cycling rule in drug_list_update(). Only the summary savers record it; for
# the others it is NA.
recent_treatment_failure <- function(res, i) {
  recent <- res[seq(max(1, i - 1), i)]
  tf <- unlist(lapply(recent, function(x) x$overall_treatment_failure))
  if (is.null(tf)) {
    return(NA_real_)
  }
  mean(tf)
}

# Copy the full model state back to R and attach the R-side parameter lists
# that a later call needs to continue the simulation.
prepare_state_save <- function(state, seed) {
  sim_save <- sim_get_state(state$sim$Ptr, seed)
  sim_save$saved_inputs <- list(
    barcode_list = state$barcode_list,
    spatial_list = state$spatial_list,
    housekeeping_list = state$housekeeping_list,
    drug_list = state$drug_list,
    vector_adaptation_list = state$vector_adaptation_list,
    nmf_list = state$nmf_list,
    core_parameter_list = state$mpl
  )
  sim_save
}

# Reduce a full state save to the human and strain components when only those
# are wanted.
collect_output <- function(sim_save, human_only,
                           human_vars = c("Infection_States", "Zetas", "Ages", "ID")) {
  if (!human_only) {
    return(sim_save)
  }
  strain_vars <- c(
    "Strain_infection_state_vectors",
    "Strain_day_of_infection_state_change_vectors",
    "Strain_barcode_vectors"
  )
  strains <- sim_save$populations_event_and_strains_List[strain_vars]
  c(sim_save$population_List[human_vars], strains)
}

#' @noRd
return_sim_memory <- function(sim) {
  meta <- attr(sim, "meta")
  sim_free(meta$Ptr, meta$seed)
  for (i in seq_along(sim)) {
    if ("Ptr" %in% names(sim[[i]])) {
      sim[[i]]$Ptr <- NULL
    }
  }

  meta$Ptr <- NULL
  if ("Ptr" %in% names(sim)) {
    sim$Ptr <- NULL
  }
  attr(sim, "meta") <- meta
  gc()
  return(sim)
}
