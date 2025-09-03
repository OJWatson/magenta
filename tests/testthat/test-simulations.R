context("simulations")

test_that("dependent plaf draws with 100% importation fixes plaf", {
  
  nl <- 6
  hd <- 1
  tl <- 5
  m <- 0.00084
  sp <- 0.22
  start_p <- 0.5
  plaf <- matrix(c(rep(c(start_p,start_p,0,0,0,0),hd),rep(0,nl*tl)),ncol=nl,byrow=TRUE)
  
  
  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
      mft_flag = FALSE,
      artemisinin_loci = 5,
      absolute_fitness_cost_flag = FALSE,
      partner_drug_ratios = 1,
      drugs = list(magenta:::drug_create_asaq()),
      cost_of_resistance = rep(1 - (0.0005*1), 6),
      sequential_cycling = -1,
      number_of_drugs = 1,
      sequential_update = -1,
      number_of_resistance_loci = 6
    )
  
  x <- list(
    EIR = 1,
    N = 1000,
    years = hd+tl,
    itn_cov = 0,
    save_lineages = TRUE,
    irs_cov = 0,
    ft = 0.2,
    num_loci = nl,
    sample_reps = 1,
    mutation_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
    mutation_treated_modifier = 100,
    mutation_rate = rep(m, nl),
    survival_percentage = sp,
    genetics_df_without_summarising = TRUE,
    spatial_incidence_matrix = c(rep(1, hd), rep(0, tl)),
    spatial_mosquitoFOI_matrix = c(rep(1, hd), rep(0, tl)),
    human_only_full_save = FALSE,
    spatial_type = "island",
    use_historic_interventions = TRUE,
    seed = as.integer(runif(1, 1, 1000000000)),
    sample_size = c(100, 1000),
    sample_states = c(1, 2, 4),
    genetics_df_without_summarising = TRUE,
    ibd_length = 1,
    update_length = 30,
    plaf = plaf,
    update_save = TRUE,
    human_update_save = TRUE,
    summary_saves_only = TRUE,
    housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE, cluster = TRUE),
    drug_list = drug_list,
    nmf_list = magenta:::nmf_list_create(nmf_flag = FALSE, prob_of_testing_nmf = 1), 
    island_imports_plaf_linked_flag = TRUE
  )
  
  out <- magenta::pipeline(EIR=x$EIR,
                           seed=x$seed,
                           save_lineages = x$save_lineages,
                           N=x$N,
                           mutation_rate = x$mutation_rate, 
                           mutation_flag = x$mutation_flag,
                           sample_size = x$sample_size,
                           sample_reps = x$sample_reps,
                           years=x$years,
                           survival_percentage=x$survival_percentage,
                           itn_cov=x$itn_cov,
                           irs_cov=x$irs_cov,
                           ft=x$ft,
                           genetics_df_without_summarising = x$genetics_df_without_summarising,
                           spatial_incidence_matrix=x$spatial_incidence_matrix,
                           spatial_mosquitoFOI_matrix=x$spatial_mosquitoFOI_matrix,
                           spatial_type=x$spatial_type,
                           use_historic_interventions=x$use_historic_interventions,
                           human_only_full_save=x$human_only_full_save,
                           ibd_length=x$ibd_length,
                           num_loci=x$num_loci,
                           sample_states = x$sample_states,
                           update_length=x$update_length,
                           update_save=x$update_save,
                           human_update_save=x$human_update_save,
                           summary_saves_only=x$summary_saves_only,
                           housekeeping_list=x$housekeeping_list,
                           island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag,
                           drug_list = x$drug_list, 
                           nmf_list = x$nmf_list, plaf = x$plaf)
  
  
  expect_equal(as.numeric(which(out[[1]]$lineage>0)), c(1,4))
  expect_false(identical(as.numeric(which(out[[length(out)-1]]$lineage>0)), c(1,4)))
  
  
})



test_that("Asymptomatic Prophylaxis Test", {
  
  nl <- 6
  hd <- 1
  tl <- 5
  m <- 0.00084
  sp <- 0.22
  start_p <- 0.5
  plaf <- matrix(c(rep(c(start_p,start_p,0,0,start_p,0),hd),rep(0,nl*tl)),ncol=nl,byrow=TRUE)
  
  
  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
      mft_flag = FALSE,
      artemisinin_loci = 5,
      absolute_fitness_cost_flag = FALSE,
      partner_drug_ratios = 1,
      drugs = list(magenta:::drug_create_asaq()),
      cost_of_resistance = rep(1 - (0.0005*1), 6),
      sequential_cycling = -1,
      number_of_drugs = 1,
      sequential_update = -1,
      number_of_resistance_loci = 6
    )
  
  # Force the drug to always be LPF
  drug_list$drugs[[1]]$lpf <- rep(0, length(drug_list$drugs[[1]]$lpf))
  
  x <- list(
    EIR = 1,
    N = 1000,
    years = hd+tl,
    itn_cov = 0,
    save_lineages = TRUE,
    irs_cov = 0,
    ft = 1,
    num_loci = nl,
    sample_reps = 1,
    mutation_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
    mutation_treated_modifier = 100,
    mutation_rate = rep(m, nl),
    survival_percentage = sp,
    genetics_df_without_summarising = TRUE,
    spatial_incidence_matrix = c(rep(1, hd), rep(0, tl)),
    spatial_mosquitoFOI_matrix = c(rep(1, hd), rep(0, tl)),
    human_only_full_save = FALSE,
    spatial_type = "island",
    use_historic_interventions = TRUE,
    seed = as.integer(runif(1, 1, 1000000000)),
    sample_size = c(100, 1000),
    sample_states = c(1, 2, 4),
    genetics_df_without_summarising = TRUE,
    ibd_length = 1,
    update_length = 30,
    plaf = plaf,
    update_save = TRUE,
    human_update_save = TRUE,
    summary_saves_only = TRUE,
    housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE, cluster = TRUE, quiet_test = FALSE),
    drug_list = drug_list,
    nmf_list = magenta:::nmf_list_create(nmf_flag = FALSE, prob_of_testing_nmf = 1), 
    island_imports_plaf_linked_flag = TRUE
  )
  
  expect_output(out <- magenta::pipeline(EIR=x$EIR,
                                         seed=x$seed,
                                         save_lineages = x$save_lineages,
                                         N=x$N,
                                         mutation_rate = x$mutation_rate, 
                                         mutation_flag = x$mutation_flag,
                                         sample_size = x$sample_size,
                                         sample_reps = x$sample_reps,
                                         years=x$years,
                                         survival_percentage=x$survival_percentage,
                                         itn_cov=x$itn_cov,
                                         irs_cov=x$irs_cov,
                                         ft=x$ft,
                                         genetics_df_without_summarising = x$genetics_df_without_summarising,
                                         spatial_incidence_matrix=x$spatial_incidence_matrix,
                                         spatial_mosquitoFOI_matrix=x$spatial_mosquitoFOI_matrix,
                                         spatial_type=x$spatial_type,
                                         use_historic_interventions=x$use_historic_interventions,
                                         human_only_full_save=x$human_only_full_save,
                                         ibd_length=x$ibd_length,
                                         num_loci=x$num_loci,
                                         sample_states = x$sample_states,
                                         update_length=x$update_length,
                                         update_save=x$update_save,
                                         human_update_save=x$human_update_save,
                                         summary_saves_only=x$summary_saves_only,
                                         housekeeping_list=x$housekeeping_list,
                                         island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag,
                                         drug_list = x$drug_list, 
                                         nmf_list = x$nmf_list, plaf = x$plaf), 
                regexp = "Asymptomatic LPF"
  )
  
  
})



test_that("Fitness Test", {
  
  nl <- 6
  hd <- 1
  tl <- 10
  m <- 0.00084
  sp <- 0.22
  start_p <- 0.5
  plaf <- matrix(c(rep(c(start_p,start_p,0,0,start_p,0),hd),rep(0,nl*tl)),ncol=nl,byrow=TRUE)
  
  
  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
      mft_flag = FALSE,
      artemisinin_loci = 5,
      absolute_fitness_cost_flag = TRUE,
      partner_drug_ratios = 1,
      drugs = list(magenta:::drug_create_asaq()),
      cost_of_resistance = rep(1 - (0.0005*1), 6),
      sequential_cycling = -1,
      number_of_drugs = 1,
      sequential_update = -1,
      number_of_resistance_loci = 6
    )
  
  # Force the drug to always be LPF
  drug_list$drugs[[1]]$lpf <- rep(0, length(drug_list$drugs[[1]]$lpf))
  
  x <- list(
    EIR = 1,
    N = 1000,
    years = hd+tl,
    itn_cov = 0,
    save_lineages = TRUE,
    irs_cov = 0,
    ft = 1,
    num_loci = nl,
    sample_reps = 1,
    mutation_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
    mutation_treated_modifier = 100,
    mutation_rate = rep(m, nl),
    survival_percentage = sp,
    genetics_df_without_summarising = TRUE,
    spatial_incidence_matrix = c(rep(1, hd), rep(0, tl)),
    spatial_mosquitoFOI_matrix = c(rep(1, hd), rep(0, tl)),
    human_only_full_save = FALSE,
    spatial_type = "island",
    use_historic_interventions = TRUE,
    seed = as.integer(runif(1, 1, 1000000000)),
    sample_size = c(100, 1000),
    sample_states = c(1, 2, 4),
    genetics_df_without_summarising = TRUE,
    ibd_length = 1,
    update_length = 30,
    plaf = plaf,
    update_save = TRUE,
    human_update_save = TRUE,
    summary_saves_only = TRUE,
    housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE, cluster = TRUE, quiet_test = TRUE),
    drug_list = drug_list,
    nmf_list = magenta:::nmf_list_create(nmf_flag = FALSE, prob_of_testing_nmf = 1), 
    island_imports_plaf_linked_flag = TRUE
  )
  
  out <- magenta::pipeline(EIR=x$EIR,
                           seed=x$seed,
                           save_lineages = x$save_lineages,
                           N=x$N,
                           mutation_rate = x$mutation_rate, 
                           mutation_flag = x$mutation_flag,
                           sample_size = x$sample_size,
                           sample_reps = x$sample_reps,
                           years=x$years,
                           survival_percentage=x$survival_percentage,
                           itn_cov=x$itn_cov,
                           irs_cov=x$irs_cov,
                           ft=x$ft,
                           genetics_df_without_summarising = x$genetics_df_without_summarising,
                           spatial_incidence_matrix=x$spatial_incidence_matrix,
                           spatial_mosquitoFOI_matrix=x$spatial_mosquitoFOI_matrix,
                           spatial_type=x$spatial_type,
                           use_historic_interventions=x$use_historic_interventions,
                           human_only_full_save=x$human_only_full_save,
                           ibd_length=x$ibd_length,
                           num_loci=x$num_loci,
                           sample_states = x$sample_states,
                           update_length=x$update_length,
                           update_save=x$update_save,
                           human_update_save=x$human_update_save,
                           summary_saves_only=x$summary_saves_only,
                           housekeeping_list=x$housekeeping_list,
                           island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag,
                           drug_list = x$drug_list, 
                           nmf_list = x$nmf_list, plaf = x$plaf)
  
  x$drug_list$cost_of_resistance <- x$drug_list$cost_of_resistance * 0.5
  x$drug_list$cost_of_resistance[1] <- 1
  
  out2 <- magenta::pipeline(EIR=x$EIR,
                            seed=x$seed,
                            save_lineages = x$save_lineages,
                            N=x$N,
                            mutation_rate = x$mutation_rate, 
                            mutation_flag = x$mutation_flag,
                            sample_size = x$sample_size,
                            sample_reps = x$sample_reps,
                            years=x$years,
                            survival_percentage=x$survival_percentage,
                            itn_cov=x$itn_cov,
                            irs_cov=x$irs_cov,
                            ft=x$ft,
                            genetics_df_without_summarising = x$genetics_df_without_summarising,
                            spatial_incidence_matrix=x$spatial_incidence_matrix,
                            spatial_mosquitoFOI_matrix=x$spatial_mosquitoFOI_matrix,
                            spatial_type=x$spatial_type,
                            use_historic_interventions=x$use_historic_interventions,
                            human_only_full_save=x$human_only_full_save,
                            ibd_length=x$ibd_length,
                            num_loci=x$num_loci,
                            sample_states = x$sample_states,
                            update_length=x$update_length,
                            update_save=x$update_save,
                            human_update_save=x$human_update_save,
                            summary_saves_only=x$summary_saves_only,
                            housekeeping_list=x$housekeeping_list,
                            island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag,
                            drug_list = x$drug_list, 
                            nmf_list = x$nmf_list, plaf = x$plaf)
  
  expect_true(sum(out[[100]]$af) > sum(out2[[100]]$af))
  
  
})


test_that("Seasonality Test", {
  
  # set up the years
  year_range <- 2010:2015
  
  # grab the interventions for region of interest
  country <- "Democratic Republic of the Congo"
  admin <- "Kinshasa"
  ints <- intervention_grab(country = country, 
                            admin = admin,
                            year_range = year_range)
  ## Requested: Kinshasa, Democratic Republic of the Congo
  ## Returned: Kinshasa, Democratic Republic of the Congo
  # grab the importation from the migration model
  spl <- spl_grab(country, admin, year_range)
  ## Requested: Kinshasa, Democratic Republic of the Congo
  ## Returned: Kinshasa_City, Democratic Republic of the Congo
  # get the intervention parmeters
  ft <- ints$ft
  itn_cov <- ints$itn_cov
  irs_cov <- ints$irs_cov
  
  # load the Malaria Atlas Project Prevalence data
  MAP_file <- system.file("extdata/MAP_PfPR_population_weighted.csv", 
                          package = "magenta")
  MAP <- read.csv(MAP_file)
  
  # work out the EIR for the prevalence in year 2000 (i.e. pre interventions)
  pfpr_micro <- MAP$X2000[which(MAP$Name==admin)]
  EIR <- magenta:::pfpr_to_eir_heuristic(ft = 0, PfPR_micro = pfpr_micro)
  ## 0.00152858753907675
  # run the simulation
  out12 <- pipeline(
    EIR = EIR,
    N = 1000,
    years = length(itn_cov),
    country = country,
    admin = admin,
    itn_cov = itn_cov,
    irs_cov = irs_cov,
    ft = ft,
    spatial_incidence_matrix = spl$incidence,
    spatial_mosquitoFOI_matrix = spl$mosquitoFOI,
    spatial_type = "island",
    human_only_full_save = TRUE,
    update_length = 30,
    update_save = TRUE,
    human_update_save = FALSE,
    summary_saves_only = FALSE
  )
  
  
  out_no_seas <- pipeline(
    EIR = EIR,
    N = 1000,
    years = length(itn_cov),
    itn_cov = itn_cov,
    irs_cov = irs_cov,
    ft = ft,
    spatial_incidence_matrix = spl$incidence,
    spatial_mosquitoFOI_matrix = spl$mosquitoFOI,
    spatial_type = "island",
    human_only_full_save = TRUE,
    update_length = 30,
    update_save = TRUE,
    human_update_save = FALSE,
    summary_saves_only = FALSE
  )
  
  inc_ses <- unlist(lapply(out12, function(x) {sum(x$Incidence %>% unlist())}))
  inc <- unlist(lapply(out_no_seas, function(x) {sum(x$Incidence %>% unlist())}))
  ranges <- magenta:::ranges(12, length(inc))
  varses <- varno <- c()
  for(i in ranges) {
    varses <- c(varses, inc_ses[i] %>% var)
    varno <- c(varno, inc[i] %>% var)
  }
  
  expect_true(mean(varses) > mean(varno))
  
})