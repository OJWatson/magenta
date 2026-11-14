# Fixed-seed characterisation tests for pipeline().
#
# Each test runs a small, fast configuration with a fixed seed and pins a
# compact summary of the result. See helper-characterisation.R for what is
# pinned and why. Keep these configurations small: together they should run in
# well under a minute.

test_that("characterisation: default single-run loggers", {
  skip_characterisation_if_unsupported()

  out <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 200, seed = 1L,
    housekeeping_list = quiet_housekeeping()
  ))

  expect_s3_class(out, "magenta_sim")
  expect_named(out, "Loggers")
  expect_type(sim_meta(out), "list")
  expect_equal(sim_meta(out)$seed, 1L)
  expect_snapshot_value(snapshot_summary(sim_output(out)), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: update_save loggers over monthly periods", {
  skip_characterisation_if_unsupported()

  out <- suppressMessages(pipeline(
    EIR = 2, years = 1, N = 200, seed = 2L,
    update_save = TRUE, update_length = 30,
    ft = 0.5, itn_cov = c(0.2), irs_cov = 0,
    housekeeping_list = quiet_housekeeping()
  ))

  expect_length(out, round(365 / 30))
  expect_equal(nrow(sim_schedule(out)), round(365 / 30))
  expect_snapshot_value(snapshot_summary(sim_output(out)), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: summary saves with COI summaries", {
  skip_characterisation_if_unsupported()

  out <- suppressMessages(pipeline(
    EIR = 5, years = 1, N = 300, seed = 3L, num_loci = 6,
    update_save = TRUE, update_length = 60,
    output = save_summary(sample_size = Inf, sample_reps = 1, mean_only = TRUE),
    housekeeping_list = quiet_housekeeping()
  ))

  first <- out[[1]]
  expect_true(all(c(
    "summary_coi", "summary_coi_detected", "summary_polygenom",
    "summary_unique", "summary_cou", "clonality", "barcode_freq",
    "overall_treatment_failure", "pcr_prev", "micro_2_10"
  ) %in% names(first)))
  expect_snapshot_value(snapshot_summary(sim_output(out)), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: resistance, importation, mutation and lineages", {
  skip_characterisation_if_unsupported()

  nl <- 4
  plaf <- matrix(c(0.5, 0.5, 0, 0, rep(0, nl)), ncol = nl, byrow = TRUE)
  drug_list <- drug_list_create(
    resistance_flag = c(FALSE, TRUE),
    number_of_resistance_loci = nl,
    artemisinin_loci = 3,
    cost_of_resistance = rep(0.99, nl),
    drugs = list(drug_create(
      prob_of_lpf = seq(1, 0.5, length.out = 2^nl),
      barcode_res_pos = seq_len(nl) - 1,
      prophylactic_pos = nl - 1
    )),
    partner_drug_ratios = 1
  )

  out <- suppressMessages(pipeline(
    EIR = 2, years = 2, N = 200, seed = 4L, num_loci = nl, ft = 0.5,
    plaf = plaf, mutation_flag = c(FALSE, TRUE), mutation_rate = rep(1e-3, nl),
    mutation_treated_modifier = 10,
    spatial_type = "island", island_imports_plaf_linked_flag = TRUE,
    spatial_incidence_matrix = c(1, 0), spatial_mosquitoFOI_matrix = c(1, 0),
    update_save = TRUE, update_length = 60,
    output = save_summary(
      genetics = "allele_freqs", lineages = TRUE,
      sample_size = c(50, 100), sample_states = c(1, 2, 4)
    ),
    drug_list = drug_list,
    nmf_list = nmf_list_create(nmf_flag = TRUE, prob_of_testing_nmf = 0.5),
    housekeeping_list = quiet_housekeeping()
  ))

  expect_true(all(c("af", "lineage", "mutations") %in% names(out[[1]])))
  expect_length(out[[1]]$lineage, 2^nl)
  expect_snapshot_value(snapshot_summary(sim_output(out)), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: seasonal admin unit with historic interventions", {
  skip_characterisation_if_unsupported()

  ints <- suppressMessages(intervention_grab(
    country = "Democratic Republic of the Congo", admin = "Kinshasa",
    year_range = 2012:2013
  ))
  spl <- suppressMessages(spl_grab(
    "Democratic Republic of the Congo", "Kinshasa", 2012:2013
  ))
  expect_snapshot_value(
    list(ints = ints, spl = snapshot_summary(spl)),
    style = "deparse", tolerance = 1e-6
  )

  out <- suppressMessages(pipeline(
    EIR = 1, N = 200, years = 2, seed = 5L,
    country = "Democratic Republic of the Congo", admin = "Kinshasa",
    itn_cov = ints$itn_cov, irs_cov = ints$irs_cov, ft = ints$ft,
    spatial_incidence_matrix = spl$incidence,
    spatial_mosquitoFOI_matrix = spl$mosquitoFOI,
    spatial_type = "island",
    update_save = TRUE, update_length = 90,
    housekeeping_list = quiet_housekeeping()
  ))

  expect_snapshot_value(snapshot_summary(sim_output(out)), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: set_up_only state and full save inputs", {
  skip_characterisation_if_unsupported()

  state <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 6L, num_loci = 3,
    set_up_only = TRUE, output = save_loggers(final = "full"),
    housekeeping_list = quiet_housekeeping()
  ))

  expect_true(all(c(
    "population_List", "populations_event_and_strains_List", "scourge_List",
    "parameters_List", "parameter_housekeeping_List", "saved_inputs"
  ) %in% names(state)))
  expect_snapshot_value(snapshot_summary(state), style = "deparse", tolerance = 1e-6)

  full <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 6L, num_loci = 3,
    update_save = TRUE, update_length = 180,
    output = save_loggers(final = "full"),
    housekeeping_list = quiet_housekeeping()
  ))
  final <- full[[length(full)]]
  expect_named(
    final$saved_inputs,
    c(
      "barcode_list", "spatial_list", "housekeeping_list", "drug_list",
      "vector_adaptation_list", "nmf_list", "core_parameter_list", "spatial_settings"
    )
  )
  expect_snapshot_value(snapshot_summary(final), style = "deparse", tolerance = 1e-6)
})

test_that("characterisation: equilibrium initialisation is stable", {
  skip_characterisation_if_unsupported()

  mpl <- model_param_list_create()
  age_vector <- magenta:::age_brackets(100, 40, TRUE)
  eq <- suppressMessages(equilibrium_init_create(
    age_vector = age_vector, het_brackets = 5, ft = 0.4, EIR = 10,
    model_param_list = mpl, quiet = TRUE
  ))
  ss <- magenta:::equilibrium_ss_create(eq)

  expect_snapshot_value(
    list(
      mv0 = eq$mv0, FOIv_eq = eq$FOIv_eq, omega = eq$omega,
      state_totals = vapply(ss[grep("mat$", names(ss))], sum, numeric(1)),
      Sv = ss$Sv, Ev = ss$Ev, Iv = ss$Iv
    ),
    style = "deparse", tolerance = 1e-8
  )
})
