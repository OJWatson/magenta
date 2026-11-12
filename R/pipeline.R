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
#' # Output
#'
#' @param output What to keep after each update period and at the end: a
#'   saver from \code{\link{save_loggers}}, \code{\link{save_summary}},
#'   \code{\link{save_humans}} or \code{\link{save_full}}, or a
#'   \code{function(state, period, params)} of your own (see
#'   \code{\link{savers}}). Default = NULL, which means \code{save_loggers()}.
#' @param update_save Boolean detailing whether the simulation stops every
#'   \code{update_length} days to extract the saver's output. Default = FALSE,
#'   which runs year by year, so that the yearly inputs (\code{ft},
#'   \code{plaf}, the flags, the importation matrices) still change when they
#'   should, and returns the final slot alone.
#' @param update_length How long each update is run for in days when
#'   \code{update_save = TRUE}. Default = 365
#' @param set_up_only Boolean for whether to return just the initialised
#'   simulation state, as a full save, without running it. Default = FALSE
#' @param saved_state_path Full file path to a saved model state to be
#'   loaded and continued (as written by \code{\link{sim_save}} or kept by a
#'   saver with \code{final = "full"}), or the object itself. Default = NULL,
#'   which will trigger initialisation. A resumed run continues with the
#'   parameter lists the saved run was using (core, drug, housekeeping, nmf,
#'   vector adaptation, barcode and spatial settings): passing a different
#'   one is an error. The yearly inputs (\code{ft}, \code{itn_cov},
#'   \code{irs_cov}, \code{plaf}, \code{mutation_flag} and the importation
#'   matrices) are for the new years and may be given; \code{plaf},
#'   \code{mutation_flag} and the importation vectors default to the saved
#'   run's last values.
#'
#' # Deprecated output flags
#'
#' These arguments are kept for one release and map onto a saver with a
#' warning saying which: pass \code{output = } instead.
#'
#' @param full_save Deprecated: \code{output = save_loggers(final = "full")}.
#' @param full_update_save Deprecated: \code{output = save_full()}.
#' @param human_only_full_save Deprecated: \code{final = "humans"}.
#' @param human_update_save,summary_saves_only,genetics_df_without_summarising,only_allele_freqs
#'   Deprecated: \code{output = save_humans()} or
#'   \code{output = save_summary(genetics = )}.
#' @param mean_only,save_lineages,sample_size,sample_states,age_breaks,sample_reps
#'   Deprecated: the arguments of \code{\link{save_summary}}.
#' @param update_save_func Deprecated: pass the function as \code{output}.
#'
#' # Parameter Lists
#'
#' @param housekeeping_list List created by  \code{\link{housekeeping_list_create}}
#' @param drug_list List created by \code{\link{drug_list_create}}
#' @param vector_adaptation_list List created by \code{\link{vector_adaptation_list_create}}
#' @param nmf_list List created by \code{\link{nmf_list_create}}
#' @param seed Random seed, set once at the start of the run. Default is
#'   random. A resumed run restores the generator state the saved run stopped
#'   at (kept in \code{meta$rng_state} and in the saved state) unless
#'   \code{seed} is given.
#' @param ... Other parameters to model_param_list_create
#'
#' \code{pipeline}
#'
#' @export
#' @inheritParams spl_create
#' @examples
#' \donttest{
#' # a small, fast run: 100 people for one year, loggers only
#' out <- pipeline(
#'   EIR = 2, years = 1, N = 100, seed = 1L,
#'   housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
#' )
#' out
#' out$output$Loggers$S
#'
#' # monthly saves of the loggers, as a data frame
#' out_m <- pipeline(
#'   EIR = 2, years = 1, N = 100, seed = 1L, update_save = TRUE, update_length = 30,
#'   housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
#' )
#' head(as.data.frame(out_m))
#' }

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
                     output = NULL,
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
  supplied <- setdiff(names(call), "")
  call$seed <- args$seed
  started <- Sys.time()

  # the deprecated output flags, and which of them the caller set
  flags <- list(
    full_save = full_save, full_update_save = full_update_save,
    human_only_full_save = human_only_full_save,
    human_update_save = human_update_save,
    genetics_df_without_summarising = genetics_df_without_summarising,
    summary_saves_only = summary_saves_only, only_allele_freqs = only_allele_freqs,
    mean_only = mean_only, save_lineages = save_lineages,
    sample_size = sample_size, sample_states = sample_states,
    age_breaks = age_breaks, sample_reps = sample_reps,
    update_save_func = update_save_func
  )
  used <- c(
    full_save = !missing(full_save), full_update_save = !missing(full_update_save),
    human_only_full_save = !missing(human_only_full_save),
    human_update_save = !missing(human_update_save),
    genetics_df_without_summarising = !missing(genetics_df_without_summarising),
    summary_saves_only = !missing(summary_saves_only),
    only_allele_freqs = !missing(only_allele_freqs),
    mean_only = !missing(mean_only), save_lineages = !missing(save_lineages),
    sample_size = !missing(sample_size), sample_states = !missing(sample_states),
    age_breaks = !missing(age_breaks), sample_reps = !missing(sample_reps),
    update_save_func = !missing(update_save_func)
  )
  saver <- resolve_saver(output, flags, used)

  # PRE-SET UP HOUSEKEEPING ------------------------ ####
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
  # The random number generator is seeded once, here. A saved state carries
  # the lists it was created with and the generator state it stopped at;
  # those win over the arguments.
  if (is.null(saved_state_path)) {
    saved_state <- NULL
    check_dots_are_model_params(...)
    mpl <- model_param_list_create(...)
    set.seed(seed)
    message("Seed set to ", seed)
  } else {
    saved <- read_saved_state(saved_state_path)
    saved_state <- saved$state
    if (use_historic_interventions) {
      supplied <- c(supplied, "spatial_incidence_matrix", "spatial_mosquitoFOI_matrix")
    }
    if (length(list(...))) {
      check_dots_are_model_params(...)
      check_resume_arg(
        "...", model_param_list_create(...), saved$inputs$core_parameter_list, supplied = "..."
      )
    }
    check_resume_arg("drug_list", drug_list, saved$inputs$drug_list, supplied,
      ignore = c("drug_choice", "next_temporal_cycle", "resistance_flag")
    )
    check_resume_arg(
      "housekeeping_list", housekeeping_list, saved$inputs$housekeeping_list, supplied
    )
    check_resume_arg("nmf_list", nmf_list, saved$inputs$nmf_list, supplied)
    check_resume_arg(
      "vector_adaptation_list", vector_adaptation_list,
      saved$inputs$vector_adaptation_list, supplied
    )
    mpl <- saved$inputs$core_parameter_list
    drug_list <- saved$inputs$drug_list
    housekeeping_list <- saved$inputs$housekeeping_list
    nmf_list <- saved$inputs$nmf_list
    vector_adaptation_list <- saved$inputs$vector_adaptation_list

    bl <- saved$inputs$barcode_list
    for (nm in c(
      "num_loci", "ibd_length", "prob_crossover", "starting_ibd", "mutation_rate",
      "mutation_treated_modifier"
    )) {
      check_resume_arg(nm, get(nm), bl[[nm]], supplied)
      assign(nm, bl[[nm]])
    }
    if (!("plaf" %in% supplied)) {
      plaf <- bl$plaf
    }
    if (!("mutation_flag" %in% supplied)) {
      mutation_flag <- bl$mutation_flag
    }

    ss <- saved$inputs$spatial_settings
    if (!is.null(ss)) {
      for (nm in names(ss)) {
        check_resume_arg(nm, get(nm), ss[[nm]], supplied)
        assign(nm, ss[[nm]])
      }
    }
    sl <- saved$inputs$spatial_list
    if (!("spatial_incidence_matrix" %in% supplied) && !is.null(sl$human_importation_rate_vector)) {
      spatial_incidence_matrix <- matrix(sl$human_importation_rate_vector, nrow = 1)
    }
    if (!("spatial_mosquitoFOI_matrix" %in% supplied) &&
      !is.null(sl$mosquito_imporation_rate_vector)) {
      spatial_mosquitoFOI_matrix <- matrix(sl$mosquito_imporation_rate_vector, nrow = 1)
    }

    if ("seed" %in% supplied || is.null(saved$rng_state)) {
      set.seed(seed)
      message("Seed set to ", seed)
    } else {
      set_rng_state(saved$rng_state)
      if (!is.null(saved$seed)) {
        seed <- saved$seed
      }
      # replaying the call must restore the generator again, not reseed
      call$seed <- NULL
      message("Random number generator restored from the saved state")
    }
  }

  # SCHEDULE --------------------------------------- ####
  # One period per update when saving periodically, otherwise one per year so
  # that the yearly inputs still take effect.
  schedule <- build_schedule(
    years = years,
    update_length = if (update_save) update_length else 365,
    ft = ft, itn_cov = itn_cov, irs_cov = irs_cov,
    resistance_flag = drug_list$resistance_flag,
    mutation_flag = mutation_flag, plaf = plaf,
    spatial_incidence = spatial_incidence_matrix,
    spatial_mosquitoFOI = spatial_mosquitoFOI_matrix,
    drug_choice = drug_list$drug_choice
  )
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
    saved_state = saved_state
  )
  params <- new_magenta_params(list(
    N = N, EIR = EIR, years = years,
    update_length = if (update_save) update_length else NA,
    country = country, admin = admin, seed = seed, num_loci = num_loci,
    core = mpl, housekeeping = housekeeping_list, drug = drug_list,
    nmf = nmf_list, vector_adaptation = vector_adaptation_list,
    barcode = state$barcode_list, spatial = state$spatial_list
  ))

  # INITIALISATION ONLY ---------------------------- ####
  if (set_up_only) {
    message("Set Up Grab")
    sim_save <- prepare_state_save(state)
    if (housekeeping_list$clear_up) {
      sim_free(state$sim$Ptr)
      gc()
    }
    return(collect_output(sim_save, human_only = saver$final == "humans"))
  }

  # DETERMINISTIC MOSQUITO PRE-SIMULATION ---------- ####
  annual <- attr(schedule, "annual")
  mu_fv <- mu_fv_create(
    eqInit = state$eqInit, ft = annual$ft * state$lpf_weight,
    itn_cov = annual$itn_cov, irs_cov = annual$irs_cov, years = years
  )

  # SIMULATION ------------------------------------- ####
  run <- run_schedule(state, schedule, mu_fv, saver, params, years, save_periods = update_save)
  state <- run$state
  res <- run$res
  schedule <- run$schedule

  final <- switch(saver$final,
    loggers = state$sim,
    full = prepare_state_save(state),
    humans = collect_output(prepare_state_save(state), human_only = TRUE)
  )
  if (update_save) {
    res[[length(res)]] <- final
  } else {
    res <- final
  }

  # FINISH ----------------------------------------- ####
  meta <- list(
    call = call, seed = seed, rng_state = get_rng_state(), Ptr = state$sim$Ptr,
    version = utils::packageVersion("magenta"),
    saver = saver$name, final = saver$final, update_save = update_save,
    times = run$comp_times, elapsed = Sys.time() - started
  )
  sim <- new_magenta_sim(res, params, schedule, meta)

  if (housekeeping_list$clear_up) {
    sim <- return_sim_memory(sim)
  }
  sim
}

# The dots of pipeline() go to model_param_list_create(); anything else is
# almost always a misspelt pipeline() argument, so say so.
check_dots_are_model_params <- function(...) {
  dots <- names(list(...))
  unknown <- setdiff(dots, names(formals(model_param_list_create)))
  if (length(dots) && (is.null(dots) || any(!nzchar(dots)) || length(unknown))) {
    stop(
      "unknown argument(s) to pipeline(): ",
      paste(if (length(unknown)) unknown else "<unnamed>", collapse = ", "),
      ". Arguments in ... are passed to model_param_list_create()",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Read a saved state from disk and split it into the C++ state, the R-side
# parameter lists it was created with, and the generator state and seed of
# the run that saved it (absent from saves made before 1.4.0).
read_saved_state <- function(saved_state_path) {
  if (is.character(saved_state_path)) {
    saved <- sim_restore(saved_state_path)
  } else {
    saved <- as_saved_state(saved_state_path)
  }
  list(
    state = unclass(saved)[cpp_contract$saved_state],
    inputs = saved$saved_inputs,
    rng_state = saved$rng_state,
    seed = saved$meta$seed
  )
}

# A resumed run uses the parameter lists of the saved run. An argument the
# caller supplied that differs from the saved one is an error rather than a
# silent override: say which, and that the saved run's value is the one used.
check_resume_arg <- function(name, current, saved, supplied, ignore = character()) {
  if (!(name %in% supplied)) {
    return(invisible(NULL))
  }
  strip <- function(x) {
    if (is.list(x)) {
      x <- unclass(x)
      x[setdiff(names(x), ignore)]
    } else {
      x
    }
  }
  if (!isTRUE(all.equal(strip(current), strip(saved), check.attributes = FALSE))) {
    stop(
      "a resumed run continues with the saved run's parameter lists, but `",
      name, "` differs from the saved value; drop it from the call",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# The state of R's random number generator, so that a run can be continued
# from exactly where it stopped.
get_rng_state <- function() {
  get0(".Random.seed", envir = globalenv(), inherits = FALSE)
}

set_rng_state <- function(state) {
  assign(".Random.seed", state, envir = globalenv())
  invisible(state)
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
# The generator was seeded (or restored) by pipeline() before this is called;
# the importation frequency vectors drawn by spl_create() and the C++
# initialisation consume it in that order.
sim_setup <- function(N, EIR, country, admin, mpl, schedule, barcode, spatial,
                      housekeeping_list, drug_list, nmf_list,
                      vector_adaptation_list, saved_state = NULL) {
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

  eqSS <- equilibrium_ss_create(eqInit = eqInit)

  state <- list(
    N = N, eqInit = eqInit, mpl = mpl, lpf_weight = lpf_weight,
    housekeeping_list = housekeeping_list, nmf_list = nmf_list,
    vector_adaptation_list = vector_adaptation_list,
    spatial_settings = spatial,
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
    state$sim <- sim_init(pl)
  } else {
    state$sim <- sim_init_saved(saved_state)
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
run_period <- function(state, days, ft, mu_fv) {
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
  state$sim <- sim_update(pl)
  state
}

# Step through the schedule one period at a time, calling the saver after
# every period but the last when save_periods is TRUE. The last period is run
# without a save so that the caller can decide what to put in the final slot.
# The drug cycling rules run between periods on the treatment failure the
# loggers report, whatever the saver. The schedule is returned with the drug
# used in each period filled in.
run_schedule <- function(state, schedule, mu_fv, saver, params, years, save_periods) {
  n_periods <- nrow(schedule)
  res <- vector("list", n_periods)
  comp_times <- rep(0, max(n_periods - 1, 0))
  treatment_failure <- rep(NA_real_, n_periods)

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
    state <- run_period(state, days, schedule$ft[i], mu_fv)
    treatment_failure[i] <- period_treatment_failure(state$sim)

    if (!last) {
      if (save_periods) {
        res[i] <- list(saver$fun(state, schedule_period(schedule, i), params))
      }
      # sequential cycling looks at the failure over the last two periods
      state$drug_list <- drug_list_update(
        state$drug_list, schedule$year[i],
        mean(treatment_failure[seq(max(1, i - 1), i)])
      )
    }
  }

  list(state = state, res = res, comp_times = comp_times, schedule = schedule)
}

# Period i of the schedule as a list, with its index and the number of
# periods, for the savers.
schedule_period <- function(schedule, i) {
  row <- lapply(unclass(schedule), function(col) col[[i]])
  c(list(i = i, n = nrow(schedule)), row)
}

# Copy the full model state back to R and attach the R-side parameter lists
# and the generator state that a later call needs to continue the simulation.
prepare_state_save <- function(state) {
  sim_save <- sim_get_state(state$sim$Ptr)
  sim_save$saved_inputs <- list(
    barcode_list = state$barcode_list,
    spatial_list = state$spatial_list,
    housekeeping_list = state$housekeeping_list,
    drug_list = state$drug_list,
    vector_adaptation_list = state$vector_adaptation_list,
    nmf_list = state$nmf_list,
    core_parameter_list = state$mpl,
    spatial_settings = state$spatial_settings
  )
  sim_save$rng_state <- get_rng_state()
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
  sim_free(sim$meta$Ptr)
  out <- sim$output
  for (i in seq_along(out)) {
    if (is.list(out[[i]]) && "Ptr" %in% names(out[[i]])) {
      out[[i]]$Ptr <- NULL
    }
  }
  if ("Ptr" %in% names(out)) {
    out$Ptr <- NULL
  }
  sim$output <- out
  sim$meta$Ptr <- NULL
  gc()
  sim
}
