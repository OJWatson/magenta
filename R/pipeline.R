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
  if (length(years) > 1) {
    years <- length(years)
  }

  # PARAMETERS ------------------------------------- ####
  # A saved state carries the lists it was created with; those win.
  if (is.null(saved_state_path)) {
    saved_state <- NULL
    mpl <- model_param_list_create(...)
  } else {
    saved <- read_saved_state(saved_state_path)
    saved_state <- saved$state
    mpl <- saved$inputs$core_parameter_list
    drug_list <- saved$inputs$drug_list
    housekeeping_list <- saved$inputs$housekeeping_list
    nmf_list <- saved$inputs$nmf_list
    vector_adaptation_list <- saved$inputs$vector_adaptation_list
  }

  # SCHEDULE --------------------------------------- ####
  schedule <- build_schedule(
    years = years,
    update_length = if (update_save) update_length else round(years * 365),
    ft = ft, itn_cov = itn_cov, irs_cov = irs_cov,
    resistance_flag = drug_list$resistance_flag,
    mutation_flag = mutation_flag, plaf = plaf,
    spatial_incidence = spatial_incidence_matrix,
    spatial_mosquitoFOI = spatial_mosquitoFOI_matrix,
    drug_choice = drug_list$drug_choice
  )
  if (!update_save) {
    # a single period has always run at the smallest non-zero ft; kept until
    # the behaviour-changing fixes in Phase 4 of R_CODE_REVIEW.md
    schedule$ft <- min(ft[ft > 0])
  }

  # INITIALISATION --------------------------------- ####
  state <- sim_setup(
    N = N, EIR = EIR, country = country, admin = admin, mpl = mpl,
    schedule = schedule,
    barcode = list(
      num_loci = num_loci, ibd_length = ibd_length,
      prob_crossover = prob_crossover, starting_ibd = starting_ibd,
      mutation_rate = mutation_rate,
      mutation_treated_modifier = mutation_treated_modifier
    ),
    spatial = list(
      spatial_type = spatial_type,
      island_imports_plaf_linked_flag = island_imports_plaf_linked_flag,
      survival_percentage = survival_percentage,
      oocyst_mean = oocyst_mean, oocyst_shape = oocyst_shape
    ),
    housekeeping_list = housekeeping_list, drug_list = drug_list,
    nmf_list = nmf_list, vector_adaptation_list = vector_adaptation_list,
    saved_state = saved_state, seed = seed
  )

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
  annual <- attr(schedule, "annual")
  mu_fv <- mu_fv_create(
    eqInit = state$eqInit, ft = annual$ft * state$lpf_weight,
    itn_cov = annual$itn_cov, irs_cov = annual$irs_cov, years = years
  )

  # SIMULATION ------------------------------------- ####
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
  run <- run_schedule(state, schedule, mu_fv, saver, years, seed)
  state <- run$state
  res <- run$res
  schedule <- run$schedule

  if (full_save || human_only_full_save) {
    sim_save <- prepare_state_save(state, seed)
    final <- collect_output(sim_save,
      human_only = human_only_full_save && (update_save || !full_save),
      human_vars = if (update_save) {
        c("Infection_States", "Zetas", "Ages", "ID")
      } else {
        c("Infection_States", "Zetas", "Ages")
      }
    )
  } else {
    final <- state$sim
  }
  if (update_save) {
    res[[length(res)]] <- final
  } else {
    res <- final
  }

  # FINISH ----------------------------------------- ####
  meta <- list(
    call = call, seed = seed, Ptr = state$sim$Ptr,
    version = utils::packageVersion("magenta"), schedule = schedule
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

# Read a saved state from disk and split it into the C++ state and the R-side
# parameter lists it was created with.
read_saved_state <- function(saved_state_path) {
  saved_state <- readRDS(saved_state_path)
  saved_inputs <- saved_state$saved_inputs
  if (is.null(saved_inputs)) {
    stop("Saved state missing saved_inputs; please re-run the simulation with full_save = TRUE")
  }
  saved_state$saved_inputs <- NULL
  list(state = saved_state, inputs = saved_inputs)
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

# Build the equilibrium, barcode and spatial inputs for the first period of
# the schedule and start the simulation, fresh or from a saved state. Returns
# the simulation "state": the C++ handle plus every R-side list the main loop
# and the savers need.
#
# The order of the random draws matters: equilibrium_init_create() is followed
# by set.seed(seed), the importation frequency vectors are then drawn by
# spl_create() from that freshly seeded generator, and the binding seeds
# again before the C++ call.
sim_setup <- function(N, EIR, country, admin, mpl, schedule, barcode, spatial,
                      housekeeping_list, drug_list, nmf_list,
                      vector_adaptation_list, saved_state = NULL, seed) {
  validate_model_param_list(mpl)
  validate_housekeeping_list(housekeeping_list)
  validate_drug_list(drug_list)
  validate_nmf_list(nmf_list)
  validate_vector_adaptation_list(vector_adaptation_list)
  annual <- attr(schedule, "annual")

  # change the ft for the equilibrium and odin to reflect 28 day failure rates
  lpfs <- unlist(lapply(drug_list$drugs, function(x) x$lpf[1]))
  lpf_weight <- weighted.mean(lpfs, drug_list$partner_drug_ratios)

  # near equilibrium initial condition
  eqInit <- equilibrium_init_create(
    age_vector = age_brackets(100, 40, TRUE),
    het_brackets = 5,
    ft = annual$ft[1] * lpf_weight,
    EIR = EIR,
    country = country,
    admin = admin,
    model_param_list = mpl
  )

  # reset seed here as there is some randomness in equilibrium
  set.seed(seed)
  eqSS <- equilibrium_ss_create(eqInit = eqInit)

  state <- list(
    N = N, eqInit = eqInit, mpl = mpl, lpf_weight = lpf_weight,
    housekeeping_list = housekeeping_list, nmf_list = nmf_list,
    vector_adaptation_list = vector_adaptation_list,
    spatial_type = spatial_type_code(spatial$spatial_type),
    island_imports_plaf_linked_flag = spatial$island_imports_plaf_linked_flag,
    survival_percentage = spatial$survival_percentage,
    oocyst_mean = spatial$oocyst_mean, oocyst_shape = spatial$oocyst_shape
  )

  state$barcode_list <- barcode_list_create(
    num_loci = barcode$num_loci,
    ibd_length = barcode$ibd_length,
    plaf = schedule$plaf[[1]],
    prob_crossover = barcode$prob_crossover,
    starting_ibd = barcode$starting_ibd,
    mutation_flag = schedule$mutation_flag[1],
    mutation_rate = barcode$mutation_rate,
    mutation_treated_modifier = barcode$mutation_treated_modifier
  )
  state$spatial_list <- period_spatial_list(state, schedule, 1)
  drug_list$resistance_flag <- schedule$resistance_flag[1]
  state$drug_list <- drug_list

  if (is.null(saved_state)) {
    pl <- param_list_simulation_init_create(
      N = N, eqSS = eqSS,
      barcode_list = state$barcode_list,
      spatial_list = state$spatial_list,
      housekeeping_list = housekeeping_list,
      drug_list = state$drug_list,
      nmf_list = nmf_list,
      vector_adaptation_list = vector_adaptation_list,
      mpl = mpl
    )
    state$sim <- sim_init(pl, seed)
  } else {
    state$sim <- sim_init_saved(saved_state, seed)
    rm(saved_state)
    gc()
  }
  state
}

# The spatial list for period i of the schedule. The importation frequency
# vectors are redrawn each time, from the generator state the previous C++
# call left behind.
period_spatial_list <- function(state, schedule, i) {
  spl_create(
    spatial_type = state$spatial_type,
    human_importation_rate_vector = schedule$importation_human[[i]],
    mosquito_imporation_rate_vector = schedule$importation_mosquito[[i]],
    cotransmission_freq_vector = ztrgeomintp(10000, 10, state$survival_percentage),
    oocyst_freq_vector = ztrnbinom(10000,
      mean = state$oocyst_mean, size = state$oocyst_shape
    ),
    plaf = schedule$plaf[[i]],
    island_imports_plaf_linked_flag = state$island_imports_plaf_linked_flag
  )
}

# Switch the yearly inputs when period i starts a new year.
apply_period_year <- function(state, schedule, i) {
  state$drug_list$resistance_flag <- schedule$resistance_flag[i]
  state$barcode_list$mutation_flag <- schedule$mutation_flag[i]
  state$spatial_list <- period_spatial_list(state, schedule, i)
  state
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

# Step through the schedule one period at a time, calling the saver after
# every period but the last. The last period is run without a save so that
# the caller can decide what to put in the final slot. The schedule is
# returned with the drug used in each period filled in.
run_schedule <- function(state, schedule, mu_fv, saver, years, seed) {
  n_periods <- nrow(schedule)
  res <- vector("list", n_periods)
  comp_times <- rep(0, max(n_periods - 1, 0))

  if (n_periods > 1) {
    message("Starting Stochastic Simulation for ", years, " years")
    p <- progress::progress_bar$new(
      format = paste0(" Running: [:bar] :percent eta: :eta"),
      total = n_periods
    )
    p_print <- progress_logging(state$housekeeping_list, res, p, initial = TRUE)
  }

  for (i in seq_len(n_periods)) {
    if (i > 1 && schedule$year[i] != schedule$year[i - 1]) {
      state <- apply_period_year(state, schedule, i)
    }
    schedule$drug_choice[i] <- state$drug_list$drug_choice
    last <- i == n_periods

    if (!last) {
      p_print <- progress_logging(state$housekeeping_list, res, p, i,
        initial = FALSE, p_print = p_print
      )
      comp_times[i] <- Sys.time()
    }

    days <- seq(schedule$day_start[i], schedule$day_end[i])
    state <- run_period(state, days, schedule$ft[i], mu_fv, seed)

    if (!last) {
      res <- saver(res, i, state$sim)
      state$drug_list <- drug_list_update(
        state$drug_list, schedule$year[i], recent_treatment_failure(res, i)
      )
    }
  }

  list(state = state, res = res, comp_times = comp_times, schedule = schedule)
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
