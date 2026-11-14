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
#' @return What the saver kept, as a \code{\link{magenta_sim}}: the final
#'   slot alone, or with \code{update_save = TRUE} a list with one element
#'   per period, carrying the run's \code{params}, \code{schedule} and
#'   \code{meta} as attributes. With \code{set_up_only = TRUE}, the
#'   initialised model state as a plain list.
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
#' out$Loggers$S
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
  # grab the function call for reproducibility, with the seed filled in so
  # that replaying it gives the same run. Which arguments were given
  # explicitly matters when resuming a saved run, so note them before the
  # seed is added to the call
  args <- append(as.list(environment()), list(...))
  call <- match.call(expand.dots = TRUE)
  supplied <- setdiff(names(call), "")
  call$seed <- args$seed
  started <- Sys.time()

  # the deprecated output flags, and which of them the caller set: the output
  # argument, or the flags, decide what is kept after each period
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

  # historic intervention grab if asked: the net, spray and treatment
  # coverage and the importation matrices on file for the admin unit, one
  # row per calendar year in the range
  if (use_historic_interventions && length(years) > 1) {
    ints <- intervention_grab(country = country, admin = admin, year_range = years)
    spl <- spl_grab(country = country, admin = admin, year_range = years)
    spatial_incidence_matrix <- spl$incidence
    spatial_mosquitoFOI_matrix <- spl$mosquitoFOI
    itn_cov <- ints$itn_cov
    irs_cov <- ints$irs_cov
    ft <- ints$ft
  }
  # a calendar year range means its length from here on
  if (length(years) > 1) {
    years <- length(years)
  }

  # PARAMETERS ------------------------------------- ####
  # The random number generator is seeded once, here. A saved state carries
  # the lists it was created with and the generator state it stopped at;
  # those win over the arguments.
  if (is.null(saved_state_path)) {
    # a fresh run: the model parameters come from the dots, and the seed is
    # set once for the whole run
    saved_state <- NULL
    check_dots_are_model_params(...)
    mpl <- model_param_list_create(...)
    set.seed(seed)
    message("Seed set to ", seed)
  } else {
    # a resumed run: the C++ state and the lists it was created with come
    # from the save. Anything the caller gave explicitly has to agree with
    # the saved value, or it is an error, so that a run cannot silently mix
    # saved and current inputs
    saved <- read_saved_state(saved_state_path)
    saved_state <- saved$state
    # the historic grab above counts as giving the matrices explicitly
    if (use_historic_interventions) {
      supplied <- c(supplied, "spatial_incidence_matrix", "spatial_mosquitoFOI_matrix")
    }
    # the model parameters, if any were given, must match the saved list
    if (length(list(...))) {
      check_dots_are_model_params(...)
      check_resume_arg(
        "...", model_param_list_create(...), saved$inputs$core_parameter_list,
        supplied = "..."
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
    # from here on the saved lists are the ones in use
    mpl <- saved$inputs$core_parameter_list
    drug_list <- saved$inputs$drug_list
    housekeeping_list <- saved$inputs$housekeeping_list
    nmf_list <- saved$inputs$nmf_list
    vector_adaptation_list <- saved$inputs$vector_adaptation_list

    # the barcode list was built from separate pipeline() arguments, so each
    # is checked and taken from the save in turn; plaf and the mutation flag
    # are yearly inputs and may be given for the new years
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

    # the spatial settings likewise (saves from before 1.4.0 have none)
    ss <- saved$inputs$spatial_settings
    if (!is.null(ss)) {
      for (nm in names(ss)) {
        check_resume_arg(nm, get(nm), ss[[nm]], supplied)
        assign(nm, ss[[nm]])
      }
    }
    # no importation matrices given: the saved run's last vectors carry on
    sl <- saved$inputs$spatial_list
    if (!("spatial_incidence_matrix" %in% supplied) && !is.null(sl$human_importation_rate_vector)) {
      spatial_incidence_matrix <- matrix(sl$human_importation_rate_vector, nrow = 1)
    }
    if (!("spatial_mosquitoFOI_matrix" %in% supplied) &&
      !is.null(sl$mosquito_imporation_rate_vector)) {
      spatial_mosquitoFOI_matrix <- matrix(sl$mosquito_imporation_rate_vector, nrow = 1)
    }

    # the generator: a seed given explicitly reseeds, otherwise the run
    # carries on from exactly where the saved run stopped (saves from before
    # 1.4.0 have no generator state and reseed too)
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
    mutation_flag = mutation_flag,
    drug_choice = drug_list$drug_choice
  )
  # the allele frequencies and importation rates are a vector per year rather
  # than a value, so they sit beside the schedule as one matrix each, padded by
  # the same rule, and are read by the period's year
  annual <- annual_inputs(
    years = years, plaf = plaf,
    spatial_incidence = spatial_incidence_matrix,
    spatial_mosquitoFOI = spatial_mosquitoFOI_matrix
  )
  # INITIALISATION --------------------------------- ####
  # build the equilibrium and the first period's lists and start the C++
  # simulation, fresh or from the saved state
  state <- sim_setup(
    N = N, EIR = EIR, country = country, admin = admin, mpl = mpl,
    schedule = schedule, annual = annual,
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
  # the record of the inputs the run was created with; it comes back with the
  # result and goes to the savers
  params <- new_magenta_params(list(
    N = N, EIR = EIR, years = years,
    update_length = if (update_save) update_length else NA,
    country = country, admin = admin, seed = seed, num_loci = num_loci,
    plaf = annual$plaf, spatial_incidence = annual$spatial_incidence,
    spatial_mosquitoFOI = annual$spatial_mosquitoFOI,
    core = mpl, housekeeping = housekeeping_list, drug = drug_list,
    nmf = nmf_list, vector_adaptation = vector_adaptation_list,
    barcode = state$barcode_list, spatial = state$spatial_list
  ))

  # INITIALISATION ONLY ---------------------------- ####
  if (set_up_only) {
    message("Set Up Grab")
    # copy the full state back to R, free the C++ memory, and return the
    # state as the saver's final asks for it
    sim_save <- prepare_state_save(state)
    if (housekeeping_list$clear_up) {
      sim_free(state$sim$Ptr)
      gc()
    }
    return(collect_output(sim_save, human_only = saver$final == "humans"))
  }

  # DETERMINISTIC MOSQUITO PRE-SIMULATION ---------- ####
  # the yearly coverage is the schedule's value in the first period of each
  # year; a year with no period is not simulated (build_schedule() checks)
  yearly <- schedule[!duplicated(schedule$year), ]
  mu_fv <- mu_fv_create(
    eqInit = state$eqInit, ft = yearly$ft * state$lpf_weight,
    itn_cov = yearly$itn_cov, irs_cov = yearly$irs_cov, years = years
  )

  # SIMULATION ------------------------------------- ####
  # run every period; the saver fills res after each one, and the schedule
  # comes back with the drug used in each period
  run <- run_schedule(state, schedule, mu_fv, saver, params, years, save_periods = update_save)
  state <- run$state
  res <- run$res
  schedule <- run$schedule

  # the last slot holds what the saver's final asks for: the loggers, the
  # full model state so that the run can be continued, or the humans
  final <- switch(saver$final,
    loggers = state$sim,
    full = prepare_state_save(state),
    humans = collect_output(prepare_state_save(state), human_only = TRUE)
  )
  # with no periodic saving the final slot is the whole result
  if (update_save) {
    res[[length(res)]] <- final
  } else {
    res <- final
  }

  # FINISH ----------------------------------------- ####
  # the record of how the run was produced: the call, the seed and where the
  # generator ended up, the version, the saver, and the timings
  meta <- list(
    call = call, seed = seed, rng_state = get_rng_state(), Ptr = state$sim$Ptr,
    version = utils::packageVersion("magenta"),
    saver = saver$name, final = saver$final, update_save = update_save,
    times = run$comp_times, elapsed = Sys.time() - started
  )
  sim <- new_magenta_sim(res, params, schedule, meta)

  # if we want all memory to be freed
  if (housekeeping_list$clear_up) {
    sim <- return_sim_memory(sim)
  }
  sim
}

#' Check that the dots of pipeline() are model parameters
#'
#' The \code{...} of \code{\link{pipeline}} go to
#' \code{\link{model_param_list_create}}; anything else is almost always a
#' misspelt \code{pipeline()} argument, so say so rather than let it pass
#' silently.
#'
#' @param ... The dots as given to \code{pipeline()}.
#' @return \code{NULL}, invisibly; an error names the unknown arguments.
#' @keywords internal
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

#' Read a saved model state for a resumed run
#'
#' Reads a saved state from disk, or takes the object itself, and splits it
#' into the C++ state, the R-side parameter lists it was created with, and the
#' generator state and seed of the run that saved it (absent from saves made
#' before 1.4.0).
#'
#' @param saved_state_path Path to an \code{.rds} file written by
#'   \code{\link{sim_save}} or kept by a saver with \code{final = "full"},
#'   or the saved state itself.
#' @return List with \code{state} (the five C++ lists), \code{inputs} (the
#'   \code{saved_inputs}), \code{rng_state} and \code{seed}.
#' @keywords internal
read_saved_state <- function(saved_state_path) {
  # a path is read and checked by sim_restore(); an object is checked only
  if (is.character(saved_state_path)) {
    saved <- sim_restore(saved_state_path)
  } else {
    saved <- as_saved_state(saved_state_path)
  }
  # the five C++ lists, then what R needs to carry on from them
  list(
    state = unclass(saved)[cpp_contract$saved_state],
    inputs = saved$saved_inputs,
    rng_state = saved$rng_state,
    seed = saved$meta$seed
  )
}

#' Check an argument against the saved run's value
#'
#' A resumed run uses the parameter lists of the saved run. An argument the
#' caller supplied that differs from the saved one is an error rather than a
#' silent override: the message says which, and that the saved value is the
#' one used. Arguments the caller did not supply are not checked.
#'
#' @param name Name of the argument, as in the call.
#' @param current The value in this call.
#' @param saved The value the saved run used.
#' @param supplied Names of the arguments the caller gave explicitly.
#' @param ignore Fields of a list argument that may differ, such as the drug
#'   cycling state the saved run moved on.
#' @return \code{NULL}, invisibly.
#' @keywords internal
check_resume_arg <- function(name, current, saved, supplied, ignore = character()) {
  # an argument the caller did not give cannot disagree
  if (!(name %in% supplied)) {
    return(invisible(NULL))
  }
  # lists are compared without the fields allowed to differ (the drug
  # cycling state the saved run moved on) and without their class
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

#' Read the state of R's random number generator
#'
#' \code{.Random.seed} as it stands, so that a run can be continued from
#' exactly where it stopped. Kept in \code{meta$rng_state} and in every full
#' save.
#'
#' @return The integer vector \code{.Random.seed}, or \code{NULL} before the
#'   generator has been used.
#' @keywords internal
get_rng_state <- function() {
  get0(".Random.seed", envir = globalenv(), inherits = FALSE)
}

#' Restore the state of R's random number generator
#'
#' @param state A vector as returned by \code{\link{get_rng_state}}.
#' @return \code{state}, invisibly.
#' @keywords internal
set_rng_state <- function(state) {
  assign(".Random.seed", state, envir = globalenv())
  invisible(state)
}

#' Spatial type as the integer code the C++ reads
#'
#' @param spatial_type \code{NULL} for no spatial structure, \code{"island"}
#'   or \code{"metapop"} (not implemented), as given to \code{pipeline()}.
#' @return \code{0} for none, \code{1} for island; a numeric input is
#'   returned as it is.
#' @keywords internal
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

#' Initialise the simulation state for a run
#'
#' Builds the equilibrium, barcode and spatial inputs for the first period of
#' the schedule and starts the simulation, fresh or from a saved state. The
#' generator was seeded (or restored) by \code{\link{pipeline}} before this
#' is called; the importation frequency vectors drawn by
#' \code{\link{spl_create}} and the C++ initialisation consume it in that
#' order.
#'
#' @param N Population size.
#' @param EIR Annual EIR the equilibrium is set up for.
#' @param country,admin Country and admin unit for seasonality, or
#'   \code{NULL}.
#' @param mpl Model parameter list from \code{\link{model_param_list_create}}.
#' @param schedule The schedule from \code{\link{build_schedule}}.
#' @param annual The per-year matrices from \code{\link{annual_inputs}}.
#' @param barcode List of the barcode arguments to \code{pipeline()}:
#'   \code{num_loci}, \code{ibd_length}, \code{prob_crossover},
#'   \code{starting_ibd}, \code{mutation_rate} and
#'   \code{mutation_treated_modifier}.
#' @param spatial List of the spatial arguments to \code{pipeline()}:
#'   \code{spatial_type}, \code{island_imports_plaf_linked_flag},
#'   \code{survival_percentage}, \code{oocyst_mean} and \code{oocyst_shape}.
#' @param housekeeping_list,drug_list,nmf_list,vector_adaptation_list The
#'   parameter lists.
#' @param saved_state The five C++ lists of a saved run to continue from, or
#'   \code{NULL} to initialise from the equilibrium.
#' @return The simulation state: a list with the C++ handle and loggers
#'   (\code{sim}), the equilibrium (\code{eqInit}), the parameter lists as
#'   they stand (\code{barcode_list}, \code{spatial_list},
#'   \code{drug_list} and the rest), the per-year matrices (\code{annual})
#'   and the spatial settings, which is everything the main loop and the
#'   savers need.
#' @keywords internal
sim_setup <- function(N, EIR, country, admin, mpl, schedule, annual, barcode, spatial,
                      housekeeping_list, drug_list, nmf_list,
                      vector_adaptation_list, saved_state = NULL) {
  # check every list before it goes anywhere near C++
  validate_model_param_list(mpl)
  validate_housekeeping_list(housekeeping_list)
  validate_drug_list(drug_list)
  validate_nmf_list(nmf_list)
  validate_vector_adaptation_list(vector_adaptation_list)

  # change the ft for the equilibrium and odin to reflect 28 day failure rates
  lpfs <- unlist(lapply(drug_list$drugs, function(x) x$lpf[1]))
  lpf_weight <- weighted.mean(lpfs, drug_list$partner_drug_ratios)

  # near equilibrium initial condition for the first year's treatment rate
  eqInit <- equilibrium_init_create(
    age_vector = age_brackets(100, 40, TRUE),
    het_brackets = 5,
    ft = schedule$ft[1] * lpf_weight,
    EIR = EIR,
    country = country,
    admin = admin,
    model_param_list = mpl
  )

  # and the starting state of the population from it
  eqSS <- equilibrium_ss_create(eqInit = eqInit)

  # everything the main loop and the savers need, in one list
  state <- list(
    N = N, eqInit = eqInit, mpl = mpl, lpf_weight = lpf_weight, annual = annual,
    housekeeping_list = housekeeping_list, nmf_list = nmf_list,
    vector_adaptation_list = vector_adaptation_list,
    spatial_settings = spatial,
    spatial_type = spatial_type_code(spatial$spatial_type),
    island_imports_plaf_linked_flag = spatial$island_imports_plaf_linked_flag,
    survival_percentage = spatial$survival_percentage,
    oocyst_mean = spatial$oocyst_mean, oocyst_shape = spatial$oocyst_shape
  )

  # the barcode list for year one; the spatial list draws the importation
  # and cotransmission frequency vectors, so it is built after the seed is set
  state$barcode_list <- barcode_list_create(
    num_loci = barcode$num_loci,
    ibd_length = barcode$ibd_length,
    plaf = annual$plaf[1, ],
    prob_crossover = barcode$prob_crossover,
    starting_ibd = barcode$starting_ibd,
    mutation_flag = schedule$mutation_flag[1],
    mutation_rate = barcode$mutation_rate,
    mutation_treated_modifier = barcode$mutation_treated_modifier
  )
  state$spatial_list <- period_spatial_list(state, schedule, 1)
  # the resistance flag is a yearly input, so the drug list carries the
  # year-one value from here on
  drug_list$resistance_flag <- schedule$resistance_flag[1]
  state$drug_list <- drug_list

  # start the C++ simulation, fresh from the equilibrium or from the saved
  # state, which can be large so is dropped once C++ has copied it
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

#' The spatial list for a period of the schedule
#'
#' The importation rates and allele frequencies of the period's year, with
#' the co-transmission and oocyst frequency vectors redrawn from the generator
#' state the previous C++ call left behind, as the old loop did at every year
#' change.
#'
#' @param state The simulation state from \code{\link{sim_setup}}.
#' @param schedule The schedule.
#' @param i The period.
#' @return A spatial list from \code{\link{spl_create}}.
#' @keywords internal
period_spatial_list <- function(state, schedule, i) {
  year <- schedule$year[i]
  spl_create(
    spatial_type = state$spatial_type,
    human_importation_rate_vector = annual_row(state$annual$spatial_incidence, year),
    mosquito_imporation_rate_vector = annual_row(state$annual$spatial_mosquitoFOI, year),
    cotransmission_freq_vector = ztrgeomintp(10000, 10, state$survival_percentage),
    oocyst_freq_vector = ztrnbinom(10000,
      mean = state$oocyst_mean, size = state$oocyst_shape
    ),
    plaf = state$annual$plaf[year, ],
    island_imports_plaf_linked_flag = state$island_imports_plaf_linked_flag
  )
}

#' Switch the yearly inputs when a period starts a new year
#'
#' @param state The simulation state.
#' @param schedule The schedule.
#' @param i The period that starts the new year.
#' @return The state with the resistance flag, mutation flag and spatial list
#'   of the new year.
#' @keywords internal
apply_period_year <- function(state, schedule, i) {
  state$drug_list$resistance_flag <- schedule$resistance_flag[i]
  state$barcode_list$mutation_flag <- schedule$mutation_flag[i]
  state$spatial_list <- period_spatial_list(state, schedule, i)
  state
}

#' Run the simulation forward over one block of days
#'
#' @param state The simulation state.
#' @param days The days to simulate, as a contiguous integer vector indexing
#'   the rows of \code{mu_fv}.
#' @param ft Treatment seeking rate for the block.
#' @param mu_fv Daily mosquito death and feeding rates for the whole run, from
#'   \code{\link{mu_fv_create}}.
#' @return The state with its \code{sim} element replaced by the C++ handle
#'   and loggers after the update.
#' @keywords internal
run_period <- function(state, days, ft, mu_fv) {
  # the update list is checked against the C++ contract before the call; the
  # mosquito rates are the days' slice of the whole run's daily series
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

#' Step through the schedule one period at a time
#'
#' The main loop of \code{\link{pipeline}}. Every period is run with
#' \code{\link{run_period}}; the saver is called after every period but the
#' last when \code{save_periods} is \code{TRUE}, the last being left for the
#' caller to fill with the saver's \code{final} choice. The drug cycling
#' rules run between periods on the treatment failure the loggers report,
#' whatever the saver.
#'
#' @param state The simulation state from \code{\link{sim_setup}}.
#' @param schedule The schedule.
#' @param mu_fv Daily mosquito death and feeding rates for the whole run.
#' @param saver The \code{magenta_saver} in use.
#' @param params The \code{magenta_params} recorded with the run, passed to
#'   the saver.
#' @param years Length of the run in years, for the progress message.
#' @param save_periods Whether to call the saver after each period.
#' @return List with the final \code{state}, the per-period results
#'   \code{res} (empty slots when \code{save_periods} is \code{FALSE}), the
#'   \code{comp_times} at which each period started and the
#'   \code{schedule} with the drug used in each period filled in.
#' @keywords internal
run_schedule <- function(state, schedule, mu_fv, saver, params, years, save_periods) {
  # one slot per period; the last is filled by the caller
  n_periods <- nrow(schedule)
  res <- vector("list", n_periods)
  comp_times <- rep(0, max(n_periods - 1, 0))
  treatment_failure <- rep(NA_real_, n_periods)

  # messaging: the progress bar is only worth having with more than one period
  if (n_periods > 1) {
    message("Starting Stochastic Simulation for ", years, " years")
    p <- progress::progress_bar$new(
      format = paste0(" Running: [:bar] :percent eta: :eta"),
      total = n_periods
    )
    p_print <- progress_logging(state$housekeeping_list, res, p, initial = TRUE)
  }

  # START MAIN SIMULATION LOOP
  for (i in seq_len(n_periods)) {
    # annual updates: a period that starts a new year switches the yearly
    # inputs before it runs
    if (i > 1 && schedule$year[i] != schedule$year[i - 1]) {
      state <- apply_period_year(state, schedule, i)
    }
    # record the drug this period uses
    schedule$drug_choice[i] <- state$drug_list$drug_choice
    last <- i == n_periods

    # messaging and timing; the last period has no save after it
    if (!last) {
      p_print <- progress_logging(state$housekeeping_list, res, p, i,
        initial = FALSE, p_print = p_print
      )
      comp_times[i] <- Sys.time()
    }

    # carry out the period and read the treatment failure from its loggers
    days <- seq(schedule$day_start[i], schedule$day_end[i])
    state <- run_period(state, days, schedule$ft[i], mu_fv)
    treatment_failure[i] <- period_treatment_failure(state$sim)

    if (!last) {
      # save what we want to save
      if (save_periods) {
        res[i] <- list(saver$fun(state, schedule_period(schedule, i), params))
      }
      # drug resistance updates: sequential cycling looks at the failure
      # over the last two periods
      state$drug_list <- drug_list_update(
        state$drug_list, schedule$year[i],
        mean(treatment_failure[seq(max(1, i - 1), i)])
      )
    }
  }

  list(state = state, res = res, comp_times = comp_times, schedule = schedule)
}

#' One period of the schedule as a list, for the savers
#'
#' @param schedule The schedule.
#' @param i The period.
#' @return The row as a named list, preceded by its index \code{i} and the
#'   number of periods \code{n}.
#' @keywords internal
schedule_period <- function(schedule, i) {
  # column by column, so a saver reads period$ft rather than a data frame row
  row <- lapply(unclass(schedule), function(col) col[[i]])
  c(list(i = i, n = nrow(schedule)), row)
}

#' The full model state, ready to be saved
#'
#' Copies the full model state back from C++ and attaches the R-side parameter
#' lists and the generator state that a later call needs to continue the
#' simulation.
#'
#' @param state The simulation state.
#' @return The five C++ lists plus \code{saved_inputs} (the parameter lists
#'   and spatial settings) and \code{rng_state}.
#' @keywords internal
prepare_state_save <- function(state) {
  # the five C++ lists, then the R-side lists as they stand at this point
  # of the run, which is what a resumed run continues with
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
  # and the generator, so that the resumed run carries on from here
  sim_save$rng_state <- get_rng_state()
  sim_save
}

#' Reduce a full state save to the human population
#'
#' @param sim_save A full state save from \code{\link{prepare_state_save}}.
#' @param human_only Whether to keep only the human and strain components.
#' @param human_vars The elements of the population list kept.
#' @return \code{sim_save} as given, or the human and strain vectors alone.
#' @keywords internal
collect_output <- function(sim_save, human_only,
                           human_vars = c("Infection_States", "Zetas", "Ages", "ID")) {
  if (!human_only) {
    return(sim_save)
  }
  # just the humans: their states, immunity and ages, and their strains
  strain_vars <- c(
    "Strain_infection_state_vectors",
    "Strain_day_of_infection_state_change_vectors",
    "Strain_barcode_vectors"
  )
  strains <- sim_save$populations_event_and_strains_List[strain_vars]
  c(sim_save$population_List[human_vars], strains)
}

#' Free the simulation's memory and drop the pointers from the result
#'
#' @param sim The \code{magenta_sim} with its C++ pointer still live.
#' @return \code{sim} with every \code{Ptr} removed, including from
#'   \code{meta}.
#' @keywords internal
return_sim_memory <- function(sim) {
  meta <- sim_meta(sim)
  sim_free(meta$Ptr)
  # the pointer is dead now, so drop it from wherever a saver kept it
  for (i in seq_along(sim)) {
    if (is.list(sim[[i]]) && "Ptr" %in% names(sim[[i]])) {
      sim[[i]]$Ptr <- NULL
    }
  }
  if ("Ptr" %in% names(sim)) {
    sim$Ptr <- NULL
  }
  meta$Ptr <- NULL
  attr(sim, "meta") <- meta
  gc()
  sim
}
