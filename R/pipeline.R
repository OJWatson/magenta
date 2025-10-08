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
                     full_update_save=FALSE,
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
  # if no seed is specified then save the seed
  set.seed(seed)
  message(paste0("Seed set to ", seed))
  message("magenta v", utils::packageVersion("magenta"))
  
  # simulation save variables
  strain_vars <- c(
    "Strain_infection_state_vectors",
    "Strain_day_of_infection_state_change_vectors",
    "Strain_barcode_vectors"
  )
  human_vars <- c("Infection_States", "Zetas", "Ages", "ID")
  
  # historic intervetnion grab if asked
  if (use_historic_interventions & length(years) > 1) {
    
    ints <- intervention_grab(country = country, admin = admin, year_range = years)
    spl <- spl_grab(country = country, admin = admin, year_range = years)
    spatial_incidence_matrix <- spl$incidence
    spatial_mosquitoFOI_matrix <- spl$mosquitoFOI
    itn_cov <- ints$itn_cov
    irs_cov <- ints$irs_cov
    ft <- ints$ft
  }
  
  # INITIALISATION --------------------------------- ####
  # If we don't have a saved state then we initialise first
  if (is.null(saved_state_path)) {
    
    # Create parameter list, changine any key parameters
    mpl <- model_param_list_create(...)
    
    # Create geometric age brackets
    age_vector <- age_brackets(100, 40, TRUE)
    
    # check to change the ft for the initial and odin to reflect 28 day failure rates
    lpfs <- unlist(lapply(drug_list$drugs, function(x) {x$lpf[1]}))
    ft_odin <- ft * weighted.mean(lpfs, drug_list$partner_drug_ratios)
    
    # Create a near equilibirum initial condition
    eqInit <- equilibrium_init_create(
      age_vector = age_vector,
      het_brackets = 5,
      ft = ft_odin[1],
      EIR = EIR,
      country = country,
      admin = admin,
      model_param_list = mpl
    )
    
    # reset seed here as there is some randomness in equlibirum (Need?)
    set.seed(seed)
    
    # Next create the starting state
    eqSS <- equilibrium_ss_create(eqInit = eqInit)
    
    # and the barcode parms list
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
    
    # handle mutations parms
    if (length(mutation_flag) == 1) {
      mutation_flag <- rep(mutation_flag, ceiling(years))
    }
    
    # spatial checks and formatting
    if (!is.null(spatial_type)) {
      if (spatial_type == "metapop") {
        # spatial_type <- 2
        # if (is.null(spatial_uuid)) {
        #   stop("Spatial_uuid is required if running spatial simulations")
        # }
        # redis_id <- paste0("magenta_", spatial_uuid)
      } else if (spatial_type == "island") {
        spatial_type <- 1
      }
    } else {
      spatial_type <- 0
    }
    
    # make spatial list
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
    
    
    # handle drug parms
    resistance_flags <- drug_list$resistance_flag
    if (length(resistance_flags) == 1) {
      resistance_flags <- rep(resistance_flags, ceiling(years))
    }
    drug_list$resistance_flag <- resistance_flags[1]
    
    # Now check and create the parameter list for use in the Rcpp simulation
    pl <- param_list_simulation_init_create(
      N = N, eqSS = eqSS,
      barcode_list = barcode_list,
      spatial_list = spatial_list,
      housekeeping_list = housekeeping_list,
      drug_list = drug_list,
      nmf_list = nmf_list,
      vector_adaptation_list = vector_adaptation_list,
      mpl = mpl
    )
  } else {

    # If we have provided the saved state then load this and then delete as can be large
    saved_state <- readRDS(saved_state_path)
    saved_inputs <- saved_state$saved_inputs
    if (is.null(saved_inputs)) {
      stop("Saved state missing saved_inputs; please re-run the simulation with full_save = TRUE")
    }
    barcode_list <- saved_inputs$barcode_list
    spatial_list <- saved_inputs$spatial_list
    housekeeping_list <- saved_inputs$housekeeping_list
    drug_list <- saved_inputs$drug_list
    vector_adaptation_list <- saved_inputs$vector_adaptation_list
    nmf_list <- saved_inputs$nmf_list
    mpl <- saved_inputs$core_parameter_list
    
    # Create geometric age brackets
    age_vector <- age_brackets(100, 40, TRUE)
    
    # check to change the ft for the initial and odin to reflect 28 day failure rates
    lpfs <- unlist(lapply(drug_list$drugs, function(x) {x$lpf[1]}))
    ft_odin <- ft * weighted.mean(lpfs, drug_list$partner_drug_ratios)
    
    # Create a near equilibirum initial condition
    eqInit <- equilibrium_init_create(
      age_vector = age_vector,
      het_brackets = 5,
      ft = ft_odin[1],
      EIR = EIR,
      country = country,
      admin = admin,
      model_param_list = mpl
    )
    
    # reset seed here as there is some randomness in equlibirum (Need?)
    set.seed(seed)
    
    # Next create the starting state
    eqSS <- equilibrium_ss_create(eqInit = eqInit)
    
    # and the barcode parms list
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
    
    # handle mutations parms
    if (length(mutation_flag) == 1) {
      mutation_flag <- rep(mutation_flag, ceiling(years))
    }
    
    # spatial checks and formatting
    if (!is.null(spatial_type)) {
      if (spatial_type == "metapop") {
        # spatial_type <- 2
        # if (is.null(spatial_uuid)) {
        #   stop("Spatial_uuid is required if running spatial simulations")
        # }
        # redis_id <- paste0("magenta_", spatial_uuid)
      } else if (spatial_type == "island") {
        spatial_type <- 1
      }
    } else {
      spatial_type <- 0
    }
    
    # make spatial list
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
    
    
    # handle drug parms
    resistance_flags <- drug_list$resistance_flag
    if (length(resistance_flags) == 1) {
      resistance_flags <- rep(resistance_flags, ceiling(years))
    }
    drug_list$resistance_flag <- resistance_flags[1]

    # clear the saved state
    saved_state$saved_inputs <- NULL
    
    # set up
    pl <- param_list_simulation_saved_init_create(savedState = saved_state)
    rm(saved_state)
    gc()
  }
  
  # Create model simulation state
  sim.out <- simulation_R(param_list = pl, seed = seed)
  
  # If we just want the set up
  # INITIALISATIION ONLY --------------------------- ####
  if (set_up_only) {
    message("Set Up Grab")
    # Now let's save the simulation in full
    pl2 <- param_list_simulation_get_create(statePtr = sim.out$Ptr)
    sim_save <- simulation_R(pl2, seed = seed)
    sim_save$saved_inputs <- list(
      barcode_list = barcode_list,
      spatial_list = spatial_list,
      housekeeping_list = housekeeping_list,
      drug_list = drug_list,
      vector_adaptation_list = vector_adaptation_list,
      nmf_list = nmf_list,
      core_parameter_list = mpl
    )
    
    if (housekeeping_list$clear_up) {
      pl5 <- param_list_simulation_finalizer_create(sim.out$Ptr)
      sim.out <- simulation_R(pl5, seed = seed)
      gc()
    }
    
    # If we want just the humans then get the keybits and save that instead
    if (full_save) {
      return(sim_save)
    }
    if (human_only_full_save) {
      Strains <- sim_save$populations_event_and_strains_List[strain_vars]
      Humans <- c(sim_save$population_List[human_vars], Strains)
      res <- Humans
      return(res)
    } else {
      res <- sim_save
      return(res)
    }
  }
  
  # if it's metapopulation sim set up the necessary redis lists
  # set this up here anyway and then less if loops later
  # METAPOP SETUP (NOT FULLY IMPLEMENTED YET) ------ #####
  imported_barcodes <- NULL
  if (spatial_type == 2) {
    stop("metapopulation simulation not implemented yet")
    
    # redis <- redux::hiredis(host = "fi--dideclusthn.dide.ic.ac.uk")
    # 
    # # create barcode as binary string
    # barcodes <- lapply(sim.out$Exported_Barcodes, as.numeric) %>%
    #   lapply(paste0, collapse = "") %>%
    #   unlist()
    # 
    # # work out here how many barcodes are going to which other simulation
    # export_proportions <- round((spatial / sum(spatial, na.rm = TRUE)) * 
    #                               round(sum(spatial * N, na.rm = TRUE)))
    # export_positions <- list()
    # import_positions <- list()
    # check_barcodes <- list()
    # epc <- 1
    # metapopulation_number <- which(is.na(spatial))
    # other_metapopulations <- seq_len(length(spatial))[-metapopulation_number]
    # 
    # # Create the necessary redis lists and push the barcodes
    # for (i in 1:length(export_proportions)) {
    #   if (!is.na(export_proportions[i])) {
    #     export_positions[[i]] <- epc:(epc + export_proportions[i] - 1)
    #     import_positions[[i]] <- paste0(i, metapopulation_number)
    #     check_barcodes[[i]] <- paste0(i, "DONE")
    #     epc <- epc + export_proportions[i]
    #     redis$SET(
    #       paste0(
    #         redis_id, "_", metapopulation_number, i
    #       ),
    #       redux::object_to_bin(barcodes[export_positions[[i]]])
    #     )
    #   } else {
    #     export_positions[[i]] <- export_proportions[i]
    #     import_positions[[i]] <- NULL
    #     check_barcodes[[i]] <- NULL
    #   }
    # }
    # 
    # # clear the blank ones
    # import_positions[[metapopulation_number]] <- NULL
    # check_barcodes[[metapopulation_number]] <- NULL
    # 
    # # Say that this metapopulation is finished
    # redis$SET(paste0(metapopulation_number, "DONE"), 1)
    # 
    # # create the redis object for grabbing other barcodes
    # redis <- redux::redis
    # get_barcodes_cmds <- lapply(import_positions, function(x) {
    #   (redis$GET(x))
    # })
    # check_done <- lapply(check_barcodes, function(x) {
    #   (redis$GET(x))
    # })
    # 
    # # sit on a while loop here untill all populations are ready
    # cannot_import <- TRUE
    # while (cannot_import) {
    #   barcode_checks <- unlist(redis$pipeline(.commands = check_done))
    #   if (sum(is.null(barcode_checks)) == 0) {
    #     cannot_import <- FALSE
    #   }
    # }
    # 
    # # once they are all ready set the done to 0
    # redis$SET(paste0(redis_id, "_", metapopulation_number, "DONE"), 0)
    # 
    # # fetch and format the barcodes for import
    # imported_barcodes <- redis$pipeline(.commands = get_barcodes_cmds) %>%
    #   lapply(redux::bin_to_object) %>%
    #   unlist()
    # imported_barcodes <- imported_barcodes %>%
    #   lapply(function(x) {
    #     strsplit(x, "") %>% unlist() %>% as.numeric() %>% as.logical()
    #   })
  }
  
  
  # DETERMINISTIC MOSQUITO PRE-SIMULATION ---------- #####
  out <- mu_fv_create(
    eqInit = eqInit, ft = ft_odin, itn_cov = itn_cov,
    irs_cov = irs_cov, years = years
  )
  
  if (length(ft) == 1 && update_save) {
    ft <- rep(ft, ceiling(years))
  }
  
  
  # SIMULATION WITH LOGGING ------------------------ #####
  # If we have specified a yearly save we iterate through
  # the total time in chunks saving the loggers at each stage
  if (update_save) {
    
    # set up results list
    res <- list()
    length(res) <- round((years * 365) / update_length)
    
    # Update times in days
    update_times <- (update_length * ((1:(round((years * 365) / update_length)))))
    
    # Vector for stroing computation time
    comp_times <- rep(0, length(res) - 1)
    
    # and set up annual checks for variables that change discretely
    year <- 1
    next_drug_cycle <- drug_list$temporal_cycling
    ft_now <- ft[year]
    
    # messaging
    message("Starting Stochastic Simulation for ", years, " years")
    p <- progress::progress_bar$new(
      format = paste0(" Running: [:bar] :percent eta: :eta"),
      total = length(res)
    )
    p_print <- progress_logging(housekeeping_list, res, p, initial = TRUE)
    
    # START MAIN SIMULATION LOOP
    if (length(res) > 1) {
      for (i in 1:(length(res) - 1)) {
        
        # messaging
        p_print <- progress_logging(housekeeping_list, res, p, i,
                                    initial = FALSE, p_print = p_print
        )
        comp_times[i] <- Sys.time()
        
        # annual updates
        if (floor((update_times[i] - update_length + 1)/365) == year) {
          
          # update the year, ft and resistance flag
          year <- year + 1
          ft_now <- ft[year]
          drug_list$resistance_flag <- resistance_flags[year]
          barcode_list$mutation_flag <- mutation_flag[year]
          
          # update the spatial list
          spatial_list <- spl_create(
            spatial_type = spatial_type,
            human_importation_rate_vector = spatial_incidence_matrix[year, ],
            mosquito_imporation_rate_vector = spatial_mosquitoFOI_matrix[year, ],
            cotransmission_freq_vector = ztrgeomintp(10000, 10, survival_percentage),
            oocyst_freq_vector = ztrnbinom(10000, mean = oocyst_mean, size = oocyst_shape),
            plaf = plaf_matrix[year, ]
          )
        }
        
        # prepare simulation for update
        pl2 <- param_list_simulation_update_create(
          years = update_length / 365,
          ft = ft_now,
          mu_vec = out$mu[1:update_length + ((i - 1) * update_length)],
          fv_vec = out$fv[1:update_length + ((i - 1) * update_length)],
          spatial_list = spatial_list,
          drug_list = drug_list,
          barcode_list = barcode_list,
          statePtr = sim.out$Ptr
        )
        
        # carry out simulation
        sim.out <- simulation_R(pl2, seed = seed)
        
        # save what we want to save
        if (is.null(update_save_func)) {
          res <- update_saves(
            res = res, i = i, sim.out = sim.out, 
            sample_states = sample_states,
            sample_size = sample_size, sample_reps = sample_reps, 
            mean_only = mean_only, barcode_list = barcode_list, 
            num_loci = num_loci, 
            genetics_df_without_summarising = genetics_df_without_summarising, 
            save_lineages = save_lineages,
            human_update_save = human_update_save, 
            summary_saves_only = summary_saves_only,
            only_allele_freqs = only_allele_freqs, 
            mpl = mpl, seed = seed, 
            full_update_save = full_update_save
          )
        } else {
          res <- update_save_func(res, i, sim.out, mpl, num_loci)
        }
        
        # spatial export
        # Metapopulation not fully implemeted yet so commenting out
        if (spatial_type == 2) {
          
          # # Push barcodes to redis
          # for (i in 1:length(export_proportions)) {
          #   if (!is.na(export_proportions[i])) {
          #     redis$SET(
          #       paste0(redis_id, "_", metapopulation_number, i),
          #       redux::object_to_bin(barcodes[export_positions[[i]]])
          #     )
          #   }
          # }
          # 
          # # once pushed set the done to 1
          # redis$SET(paste0(redis_id, "_", metapopulation_number, "DONE"), 1)
          # 
          # # sit on a while loop here untill all populations are ready
          # cannot_import <- TRUE
          # while (cannot_import) {
          #   barcode_checks <- redis$pipeline(.commands = check_done) %>% unlist()
          #   if (sum(barcode_checks == 0) == 0) {
          #     cannot_import <- FALSE
          #   }
          # }
          # 
          # # once they are all ready set the done to 0
          # redis$SET(paste0(redis_id, "_", metapopulation_number, "DONE"), 0)
          # 
          # # fetch and format the barcodes for import
          # imported_barcodes <- redis$pipeline(.commands = get_barcodes_cmds) %>%
          #   lapply(redux::bin_to_object) %>%
          #   unlist()
          # imported_barcodes <- imported_barcodes %>%
          #   lapply(function(x) {
          #     strsplit(x, "") %>% unlist() %>% as.numeric() %>% as.logical()
          #   })
        }
        
        # drug resistance updates
        if (i == 1) {
          drug_list <- drug_list_update(
            drug_list, year, res[[i]]$overall_treatment_failure
          )
        } else {
          drug_list <- drug_list_update(
            drug_list, year,
            mean(c(res[[i - 1]]$overall_treatment_failure, res[[i]]$overall_treatment_failure))
          )
        }
      }
    }
    
    
    # last rep may not occupy a full update length
    temp_mu <- out$mu[1:update_length + ((i) * update_length)]
    abridged_length <- sum(!is.na(temp_mu))
    
    # final run
    pl2 <- param_list_simulation_update_create(
      years = abridged_length / 365, ft = ft_now,
      mu_vec = out$mu[1:abridged_length + ((i) * update_length)],
      fv_vec = out$fv[1:abridged_length + ((i) * update_length)],
      spatial_list = spatial_list,
      drug_list = drug_list,
      barcode_list = barcode_list,
      statePtr = sim.out$Ptr
    )
    sim.out <- simulation_R(pl2, seed = seed)
    
    # If we have specified a full save then we grab that and save
    # it or just the human bits of interest
    if (full_save || human_only_full_save) {
      
      # Now let's save the simulation in full
      pl2 <- param_list_simulation_get_create(statePtr = sim.out$Ptr)
      sim_save <- simulation_R(pl2, seed = seed)
      sim_save$saved_inputs <- list(
        barcode_list = barcode_list,
        spatial_list = spatial_list,
        housekeeping_list = housekeeping_list,
        drug_list = drug_list,
        vector_adaptation_list = vector_adaptation_list,
        nmf_list = nmf_list,
        core_parameter_list = mpl
      )
      
      # If we want just the humans then get the keybits and save that instead
      if (human_only_full_save) {
        Strains <- sim_save$populations_event_and_strains_List[strain_vars]
        Humans <- c(sim_save$population_List[human_vars], Strains)
        res[[length(res)]] <- Humans
      }
      else {
        res[[length(res)]] <- sim_save
      }
    }
    else {
      # If we don't want a full save then just save the Loggers as usual
      res[[length(res)]] <- sim.out
    }
    
    # OR SIMULATION NO LOGGING ----------------------- #####
  } else {
    
    # Set up update for years long
    pl2 <- param_list_simulation_update_create(
      years = years, ft = min(ft[ft > 0]),
      mu_vec = out$mu,
      fv_vec = out$fv,
      spatial_list = spatial_list,
      drug_list = drug_list,
      barcode_list = barcode_list,
      statePtr = sim.out$Ptr
    )
    
    sim.out <- simulation_R(param_list = pl2, seed = seed)
    
    # If we have specified a full save or human save then we grab that
    # and save it or just the human bits of interest
      if (full_save || human_only_full_save) {

        # Now let's save the simulation in full
        pl2 <- param_list_simulation_get_create(statePtr = sim.out$Ptr)
        sim_save <- simulation_R(pl2, seed = seed)
        sim_save$saved_inputs <- list(
          barcode_list = barcode_list,
          spatial_list = spatial_list,
          housekeeping_list = housekeeping_list,
          drug_list = drug_list,
          vector_adaptation_list = vector_adaptation_list,
          nmf_list = nmf_list,
          core_parameter_list = mpl
        )

        # If we want just the humans then get the keybits and save that instead
        if (human_only_full_save & !full_save) {
          human_vars <- c("Infection_States", "Zetas", "Ages")
          Strains <- sim_save$populations_event_and_strains_List[strain_vars]
        Humans <- c(sim_save$population_List[human_vars], Strains)
        res <- Humans
      }
      else {
        res <- sim_save
      }
    }
    else {
      # If we don't want a full save then just save the Loggers as usual
      res <- sim.out
    }
  }
  
  # FINISH ----------------------------------------- #####
  
  # Save the seed and function call
  meta <- list()
  meta$call <- call
  meta$seed <- seed
  meta$Ptr <- sim.out$Ptr
  meta$version <- utils::packageVersion("magenta")
  
  # append times if an update simulation was done
  if (update_save) {
    meta$times <- comp_times
  }
  
  # add this as an attribute
  attr(res, "meta") <- meta
  
  # if we want all memory to be freed
  if (housekeeping_list$clear_up) {
    res <- return_sim_memory(res)
    meta$Ptr <- NULL
  }
  return(res)
}


#' @noRd
return_sim_memory <- function(sim) {
  
  meta <- attr(sim, "meta")
  pl <- param_list_simulation_finalizer_create(meta$Ptr)
  sim_fin <- simulation_R(pl, seed = meta$seed)
  for(i in 1:length(sim)){
    if("Ptr" %in% names(sim[[i]])) {
      sim[[i]]$Ptr <- NULL 
    }
  }
  
  meta$Ptr <- NULL
  if("Ptr" %in% names(sim)) {
    sim$Ptr <- NULL
  }
  attr(sim, "meta") <- meta
  gc()
  return(sim)
}
  