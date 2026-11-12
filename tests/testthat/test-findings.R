# The behaviour-changing fixes from R_CODE_REVIEW.md section 4.6, applied
# after review (see the review's section 7): one seed per run, yearly periods
# without periodic saving, resuming from the saved run's lists, and the
# spl_grab() year range.

test_that("the generator is seeded once, so periods do not replay the same stream", {
  skip_on_cran()
  # two identical states updated for the same ten days differ once the
  # generator has been moved on between them
  a <- tiny_state(seed = 3L)
  b <- tiny_state(seed = 3L)
  mu_fv <- constant_mu_fv(365)
  a <- magenta:::run_period(a, 1:10, 0.4, mu_fv)
  runif(1000)
  b <- magenta:::run_period(b, 1:10, 0.4, mu_fv)
  expect_false(identical(a$sim$Loggers, b$sim$Loggers))
  magenta:::sim_free(a$sim$Ptr)
  magenta:::sim_free(b$sim$Ptr)
})

test_that("the same seed reproduces a run and different seeds differ", {
  skip_on_cran()
  run <- function(seed) {
    suppressMessages(pipeline(
      EIR = 2, years = 1, N = 80, seed = seed, num_loci = 2,
      update_save = TRUE, update_length = 90,
      housekeeping_list = quiet_housekeeping_test()
    ))
  }
  one <- run(7L)
  again <- run(7L)
  other <- run(8L)
  expect_identical(as.data.frame(one), as.data.frame(again))
  expect_false(identical(as.data.frame(one), as.data.frame(other)))
  expect_identical(one$meta$rng_state, again$meta$rng_state)
  expect_false(identical(one$meta$rng_state, other$meta$rng_state))
})

test_that("a resumed run restores the generator state and so is reproducible", {
  skip_on_cran()
  first <- suppressMessages(pipeline(
    EIR = 2, years = 1, N = 80, seed = 9L, num_loci = 2,
    output = save_loggers(final = "full"),
    housekeeping_list = quiet_housekeeping_test()
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  sim_save(first, path)
  expect_identical(first$output$rng_state, first$meta$rng_state)
  expect_identical(sim_restore(path)$rng_state, first$meta$rng_state)

  resume <- function(...) {
    suppressMessages(pipeline(
      EIR = 2, years = 1, N = 80, num_loci = 2, saved_state_path = path,
      update_save = TRUE, update_length = 90,
      housekeeping_list = quiet_housekeeping_test(), ...
    ))
  }
  r1 <- resume()
  r2 <- resume()
  expect_identical(as.data.frame(r1), as.data.frame(r2))
  expect_equal(r1$meta$seed, 9L)
  expect_null(r1$meta$call$seed)
  # an explicit seed reseeds instead, and changes the result
  r3 <- resume(seed = 10L)
  expect_equal(r3$meta$seed, 10L)
  expect_false(identical(as.data.frame(r1), as.data.frame(r3)))
  # a run continued in one go gives the same loggers as the saved run's end
  expect_equal(r1$output[[1]]$Time, 366 + 90)
})

test_that("yearly inputs apply year by year without periodic saving", {
  skip_on_cran()
  ramp <- seq(0.1, 0.28, by = 0.02)
  out <- suppressMessages(pipeline(
    EIR = 1, years = 10, N = 60, seed = 4L, ft = ramp,
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(nrow(out$schedule), 10)
  expect_equal(out$schedule$ft, ramp)
  expect_equal(out$schedule$year, 1:10)
  expect_equal(attr(out$schedule, "annual")$ft, ramp)
  expect_false(out$meta$update_save)
  expect_true(is.na(out$params$update_length))
  # only the final loggers come back
  expect_named(out$output, "Loggers")
  expect_equal(out$output$Loggers$Time, 10 * 365 + 1)

  # an all-zero ft is simply zero
  out0 <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 4L, ft = 0,
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(out0$schedule$ft, 0)
  expect_equal(nrow(out0$schedule), 1)
})

test_that("temporal drug cycling happens without periodic saving", {
  skip_on_cran()
  dl <- drug_list_create(
    resistance_flag = TRUE, number_of_drugs = 2, temporal_cycling = 2,
    drugs = list(drug_create_al(), drug_create_dhappq()), partner_drug_ratios = c(0.5, 0.5)
  )
  out <- suppressMessages(pipeline(
    EIR = 1, years = 5, N = 60, seed = 4L, ft = 0.4, drug_list = dl,
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(out$schedule$drug_choice, c(0, 0, 1, 1, 0))
  # sequential cycling needs the summary saver's treatment failure; without
  # it the rule sees NA and does not switch
  dl_seq <- drug_list_create(
    resistance_flag = TRUE, number_of_drugs = 2, sequential_cycling = 0.1,
    drugs = list(drug_create_al(), drug_create_dhappq()), partner_drug_ratios = c(0.5, 0.5)
  )
  expect_equal(magenta:::drug_list_update(dl_seq, 2, NA_real_)$drug_choice, 0)
})

test_that("a resumed run continues with the saved run's parameter lists", {
  skip_on_cran()
  first <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 5L, num_loci = 3, plaf = rep(0.1, 3),
    output = save_loggers(final = "full"),
    drug_list = drug_list_create(
      drug_choice = 0, number_of_drugs = 1, drugs = list(drug_create_al()),
      number_of_resistance_loci = 3, cost_of_resistance = rep(1, 3)
    ),
    housekeeping_list = quiet_housekeeping_test()
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  sim_save(first, path)
  expect_named(first$output$saved_inputs, c(
    "barcode_list", "spatial_list", "housekeeping_list", "drug_list",
    "vector_adaptation_list", "nmf_list", "core_parameter_list", "spatial_settings"
  ))

  # nothing supplied: everything, including num_loci and plaf, comes from the save
  resumed <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, saved_state_path = path,
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(resumed$params$num_loci, 3)
  expect_equal(resumed$params$barcode$plaf, rep(0.1, 3))
  expect_equal(resumed$params$drug, first$params$drug)
  expect_identical(resumed$params$drug$drugs[[1]], drug_create_al())
  expect_equal(resumed$params$vector_adaptation, first$params$vector_adaptation)
  expect_equal(resumed$params$core, first$params$core)

  # the same list again is fine; a different one is an error
  same <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, saved_state_path = path, num_loci = 3,
    drug_list = first$params$drug, nmf_list = nmf_list_create(),
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(same$params$num_loci, 3)
  expect_error(
    suppressMessages(pipeline(
      EIR = 1, years = 1, N = 60, saved_state_path = path,
      drug_list = drug_list_create(), housekeeping_list = quiet_housekeeping_test()
    )),
    "`drug_list` differs from the saved value"
  )
  expect_error(
    suppressMessages(pipeline(
      EIR = 1, years = 1, N = 60, saved_state_path = path, num_loci = 4,
      housekeeping_list = quiet_housekeeping_test()
    )),
    "`num_loci` differs from the saved value"
  )
  expect_error(
    suppressMessages(pipeline(
      EIR = 1, years = 1, N = 60, saved_state_path = path, eta = 0.001,
      housekeeping_list = quiet_housekeeping_test()
    )),
    "`...` differs from the saved value"
  )
  # yearly inputs are for the new years and may change
  new_plaf <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, saved_state_path = path, plaf = rep(0.9, 3),
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(new_plaf$params$barcode$plaf, rep(0.9, 3))
})

test_that("spl_grab keeps the year 2000 in a range that ends there", {
  cc <- "Democratic Republic of the Congo"
  six <- suppressMessages(spl_grab(cc, "Kinshasa", 1995:2000))
  expect_equal(nrow(six$incidence), 6)
  expect_true(all(six$incidence[1:5, ] == 0))
  one <- suppressMessages(spl_grab(cc, "Kinshasa", 2000:2000))
  expect_equal(nrow(one$incidence), 1)
  expect_equal(six$incidence[6, ], one$incidence[1, ])
  expect_equal(six$mosquitoFOI[6, ], one$mosquitoFOI[1, ])
})

test_that("misspelt pipeline() arguments are reported as such", {
  expect_error(
    suppressMessages(pipeline(EIR = 1, years = 1, N = 50, update_lenght = 30)),
    "unknown argument\\(s\\) to pipeline\\(\\): update_lenght"
  )
})

test_that("a resumed island run keeps its importation rates and spatial settings", {
  skip_on_cran()
  cc <- "Democratic Republic of the Congo"
  spl <- suppressMessages(spl_grab(cc, "Kinshasa", 2012:2013))
  first <- suppressMessages(pipeline(
    EIR = 1, years = 2, N = 60, seed = 36L, num_loci = 2,
    country = cc, admin = "Kinshasa", spatial_type = "island",
    spatial_incidence_matrix = spl$incidence,
    spatial_mosquitoFOI_matrix = spl$mosquitoFOI,
    output = save_loggers(final = "full"),
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(first$output$saved_inputs$spatial_settings$spatial_type, "island")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  sim_save(first, path)

  # no matrices given: the saved run's last importation vectors carry on
  resumed <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, country = cc, admin = "Kinshasa",
    saved_state_path = path, housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(resumed$params$spatial$spatial_type, 1)
  expect_equal(
    resumed$schedule$importation_human[[1]],
    first$output$saved_inputs$spatial_list$human_importation_rate_vector
  )
  expect_equal(
    resumed$schedule$importation_mosquito[[1]],
    first$output$saved_inputs$spatial_list$mosquito_imporation_rate_vector
  )
  expect_error(
    suppressMessages(pipeline(
      EIR = 1, years = 1, N = 60, saved_state_path = path, spatial_type = NULL,
      housekeeping_list = quiet_housekeeping_test()
    )),
    "`spatial_type` differs from the saved value"
  )
  # historic interventions on resume supply the matrices for the new years
  historic <- suppressMessages(pipeline(
    EIR = 1, years = 2014:2015, N = 60, country = cc, admin = "Kinshasa",
    use_historic_interventions = TRUE, saved_state_path = path,
    housekeeping_list = quiet_housekeeping_test()
  ))
  spl2 <- suppressMessages(spl_grab(cc, "Kinshasa", 2014:2015))
  expect_equal(historic$schedule$importation_human[[1]], spl2$incidence[1, ])
})

test_that("spl_grab gives all-zero rows for a range that ends before 2000", {
  early <- suppressMessages(spl_grab("Democratic Republic of the Congo", "Kinshasa", 1996:1999))
  expect_equal(dim(early$incidence), c(4, 593))
  expect_true(all(early$incidence == 0))
})
