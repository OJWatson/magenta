# Longer end-to-end runs of the resistance, fitness and seasonality machinery.
# These take a few minutes together and are skipped on CRAN.

resistance_drug_list <- function(hd, tl, absolute = FALSE, always_lpf = FALSE) {
  dl <- drug_list_create(
    resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
    mft_flag = FALSE,
    artemisinin_loci = 5,
    absolute_fitness_cost_flag = absolute,
    partner_drug_ratios = 1,
    drugs = list(drug_create_asaq()),
    cost_of_resistance = rep(1 - (0.0005 * 1), 6),
    sequential_cycling = -1,
    number_of_drugs = 1,
    sequential_update = -1,
    number_of_resistance_loci = 6
  )
  if (always_lpf) {
    dl$drugs[[1]]$lpf <- rep(0, length(dl$drugs[[1]]$lpf))
  }
  dl
}

# An island simulation with one year of importation at the given plaf and
# then tl years of resistance and mutation, saving allele and lineage
# frequencies every 30 days.
run_resistance_sim <- function(hd, tl, plaf_row, drug_list, ft, seed,
                               housekeeping = housekeeping_list_create(
                                 quiet = TRUE, cluster = TRUE
                               )) {
  nl <- 6
  pipeline(
    EIR = 1, N = 1000, years = hd + tl, seed = seed,
    itn_cov = 0, irs_cov = 0, ft = ft,
    num_loci = nl, ibd_length = 1,
    plaf = matrix(c(rep(plaf_row, hd), rep(0, nl * tl)), ncol = nl, byrow = TRUE),
    mutation_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
    mutation_treated_modifier = 100,
    mutation_rate = rep(0.00084, nl),
    survival_percentage = 0.22,
    spatial_type = "island",
    spatial_incidence_matrix = c(rep(1, hd), rep(0, tl)),
    spatial_mosquitoFOI_matrix = c(rep(1, hd), rep(0, tl)),
    island_imports_plaf_linked_flag = TRUE,
    update_save = TRUE, update_length = 30,
    output = save_summary(
      genetics = "allele_freqs", lineages = TRUE,
      sample_size = c(100, 1000), sample_states = c(1, 2, 4), sample_reps = 1
    ),
    housekeeping_list = housekeeping,
    drug_list = drug_list,
    nmf_list = nmf_list_create(nmf_flag = FALSE, prob_of_testing_nmf = 1)
  )
}

test_that("dependent plaf draws with 100% importation fixes plaf", {
  skip_on_cran()
  out <- run_resistance_sim(
    hd = 1, tl = 5, plaf_row = c(0.5, 0.5, 0, 0, 0, 0),
    drug_list = resistance_drug_list(1, 5), ft = 0.2,
    seed = as.integer(runif(1, 1, 1000000000))
  )
  res <- out$output
  expect_equal(as.numeric(which(res[[1]]$lineage > 0)), c(1, 4))
  expect_false(identical(as.numeric(which(res[[length(res) - 1]]$lineage > 0)), c(1, 4)))
})

test_that("Asymptomatic Prophylaxis Test", {
  skip_on_cran()
  expect_output(
    run_resistance_sim(
      hd = 1, tl = 5, plaf_row = c(0.5, 0.5, 0, 0, 0.5, 0),
      drug_list = resistance_drug_list(1, 5, always_lpf = TRUE), ft = 1,
      seed = as.integer(runif(1, 1, 1000000000)),
      housekeeping = housekeeping_list_create(quiet = TRUE, cluster = TRUE, quiet_test = FALSE)
    ),
    regexp = "Asymptomatic LPF"
  )
})

test_that("Fitness Test", {
  skip_on_cran()
  seed <- as.integer(runif(1, 1, 1000000000))
  dl <- resistance_drug_list(1, 10, absolute = TRUE, always_lpf = TRUE)
  out <- run_resistance_sim(
    hd = 1, tl = 10, plaf_row = c(0.5, 0.5, 0, 0, 0.5, 0),
    drug_list = dl, ft = 1, seed = seed
  )
  dl$cost_of_resistance <- dl$cost_of_resistance * 0.5
  dl$cost_of_resistance[1] <- 1
  out2 <- run_resistance_sim(
    hd = 1, tl = 10, plaf_row = c(0.5, 0.5, 0, 0, 0.5, 0),
    drug_list = dl, ft = 1, seed = seed
  )
  expect_true(sum(out$output[[100]]$af) > sum(out2$output[[100]]$af))
})

test_that("Seasonality Test", {
  skip_on_cran()
  year_range <- 2010:2015
  country <- "Democratic Republic of the Congo"
  admin <- "Kinshasa"
  ints <- intervention_grab(country = country, admin = admin, year_range = year_range)
  spl <- spl_grab(country, admin, year_range)

  MAP_file <- system.file("extdata/MAP_PfPR_population_weighted.csv", package = "magenta")
  MAP <- read.csv(MAP_file)
  pfpr_micro <- MAP$X2000[which(MAP$Name == admin)]
  EIR <- magenta:::pfpr_to_eir_heuristic(ft = 0, PfPR_micro = pfpr_micro)

  seasonal_run <- function(seasonal) {
    pipeline(
      EIR = EIR, N = 1000, years = length(ints$itn_cov),
      country = if (seasonal) country else NULL,
      admin = if (seasonal) admin else NULL,
      itn_cov = ints$itn_cov, irs_cov = ints$irs_cov, ft = ints$ft,
      spatial_incidence_matrix = spl$incidence,
      spatial_mosquitoFOI_matrix = spl$mosquitoFOI,
      spatial_type = "island",
      update_save = TRUE, update_length = 30,
      output = save_loggers(final = "humans")
    )
  }
  out12 <- seasonal_run(TRUE)
  out_no_seas <- seasonal_run(FALSE)

  incidence <- function(sim) {
    unlist(lapply(sim$output, function(x) sum(unlist(x$Incidence))))
  }
  inc_ses <- incidence(out12)
  inc <- incidence(out_no_seas)
  ranges <- magenta:::ranges(12, length(inc))
  varses <- varno <- c()
  for (i in ranges) {
    varses <- c(varses, var(inc_ses[i]))
    varno <- c(varno, var(inc[i]))
  }

  expect_true(mean(varses) > mean(varno))
})

test_that("Gametocyte sterilisation runs and leaves the default untouched", {

  # NOTE ON WHAT IS NOT TESTED HERE.
  #
  # An earlier version of this test asserted that switching sterilisation on
  # raised the final artemisinin resistant allele frequency. That assertion is
  # not supported: over 14 seeds the paired difference is +0.045 (p = 0.30),
  # and the same test also passes on the code with the sterilisation condition
  # the other way round. It passed only because seeds 1 to 5 happened to agree.
  #
  # The reason the effect is unmeasurable is structural. The sterilisation
  # multiplier is applied to m_active_strain_contribution, which sample1()
  # uses as relative weights among the strains that one host carries. A host
  # whose strains are all sensitive has every weight scaled by the same factor,
  # so the scaling cancels, and a host with a single strain returns early
  # before the weights are used at all. Sterilisation therefore only reweights
  # which strain a mixed host transmits; it never reduces how infectious that
  # host is. Until that is addressed there is no population level signal to
  # assert on, so this test only checks that the option runs and that it is
  # genuinely off by default.

  nl <- 6
  art_locus <- 4          # kelch13 580Y, zero indexed

  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = TRUE,
      number_of_resistance_loci = nl,
      artemisinin_loci = art_locus,
      cost_of_resistance = rep(1, nl),
      drugs = list(magenta:::drug_create_al()),
      number_of_drugs = 1,
      partner_drug_ratios = 1
    )

  run <- function(vector_adaptation_list) {
    magenta::pipeline(
      EIR = 20, N = 500, years = 2, ft = 0.8,
      num_loci = nl, plaf = rep(0.15, nl),
      drug_list = drug_list,
      vector_adaptation_list = vector_adaptation_list,
      seed = 1L
    )
  }

  off <- run(magenta:::vector_adaptation_list_create(vector_adaptation_loci = art_locus))
  on <- run(magenta:::vector_adaptation_list_create(
    vector_adaptation_loci = art_locus,
    gametocyte_sterilisation_flag = TRUE,
    gametocyte_sterilisation = 0.01,
    oocyst_reduction_by_artemisinin = 1
  ))

  # both complete and report a sane prevalence
  expect_true(all(vapply(list(off, on), function(x) {
    p <- sum(unlist(x$Loggers[c("D", "A", "U", "T")]))
    p > 0 && p < 1
  }, logical(1))))

  # the flag is off by default, so the default list must not switch it on
  expect_false(
    magenta:::vector_adaptation_list_create(
      vector_adaptation_loci = art_locus
    )$gametocyte_sterilisation_flag
  )

})
