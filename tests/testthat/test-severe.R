test_that("severe probability matches a hand calculation", {

  theta <- function(age, iva, ivm) {
    fv <- 1 - (1 - 0.141195) / (1 + (age / 2493.41)^2.91282)
    0.0749886 * (0.0001191 + (1 - 0.0001191) /
                   (1 + fv * ((iva + ivm) / 1.09629)^2.00048))
  }

  for (x in list(c(0, 0, 0), c(365, 0, 0.5), c(3650, 5, 0), c(20 * 365, 30, 0))) {
    expect_equal(test_severe_probability(x[1], x[2], x[3]),
                 theta(x[1], x[2], x[3]))
  }

  # no immunity gives theta0 and immunity reduces it
  expect_equal(test_severe_probability(1000, 0, 0), 0.0749886)
  expect_lt(test_severe_probability(1000, 10, 0), test_severe_probability(1000, 1, 0))

})

test_that("severe outcomes follow the postie formula", {

  out <- severe_outcomes(c(0.01, 0.02), ft = 0.5)
  sev <- c(0.01, 0.02) / 0.8 * (0.5 * 0.58 + 0.5) / (0.1 * 0.58 + 0.9)
  expect_equal(out$severe, sev)
  expect_equal(out$severe_hospital, sev * 0.8)
  expect_equal(out$mortality, sev * 0.8 * 0.065 + sev * 0.2 * 0.6)
  expect_equal(out$severe_model, c(0.01, 0.02))
  expect_error(severe_outcomes(0.01, ft = 0.5, treatment_scaler = 2))

})

test_that("equilibrium IVA uses the ICA recursion", {

  # with the ICA boosting and decay, IVA is identical to ICA
  mpl <- model_param_list_create(uVA = 6.06349, dVA = 10950)
  eq <- equilibrium_init_create(age_vector = age_brackets(100, 40, TRUE),
                                het_brackets = 5, EIR = 20, ft = 0.4,
                                model_param_list = mpl, country = NULL,
                                admin = NULL, quiet = TRUE)
  expect_equal(eq$init_IVA, eq$init_ICA)
  expect_equal(eq$IVM_init_eq, eq$ICM_init_eq * 0.195768 / mpl$PM)

})

test_that("severe outputs are produced and sensible", {
  skip_on_cran()

  out <- pipeline(EIR = 20, N = 2000, years = 1, num_loci = 2, seed = 3L,
                  housekeeping_list = housekeeping_list_create(quiet = TRUE))
  lg <- out$Loggers

  expect_true(all(c("Severe", "Severe_05", "IVA", "IVM") %in% names(lg)))
  expect_gt(lg$Severe, 0)
  expect_gte(lg$Severe_05, 0)
  expect_lte(lg$Severe_05, lg$Severe)
  expect_length(lg$IVA, 2000)
  expect_true(all(lg$IVA >= 0) && all(lg$IVM >= 0))

  # severe immunity builds with age and maternal immunity wanes
  expect_gt(mean(lg$IVA[lg$Ages > 20 * 365]), mean(lg$IVA[lg$Ages < 5 * 365]))
  expect_lt(mean(lg$IVM[lg$Ages > 365]), mean(lg$IVM[lg$Ages < 90]))

  # boosting and decay hold severe immunity near its equilibrium
  eq <- equilibrium_init_create(age_vector = age_brackets(100, 40, TRUE),
                                het_brackets = 5, EIR = 20, ft = 0.4,
                                model_param_list = model_param_list_create(),
                                country = NULL, admin = NULL, quiet = TRUE)
  eq_mean <- sum(sweep(eq$init_IVA[, , 1], 2, eq$het_wt, "*") * eq$den)
  expect_equal(mean(lg$IVA), eq_mean, tolerance = 0.25)

  # rates helper
  r <- severe_rates(lg)
  expect_equal(r$age_group, c("all", "0-5"))
  expect_true(all(r$severe > 0 & r$severe < r$clinical))

})

test_that("equilibrium and saved states without severe immunity still run", {
  skip_on_cran()

  hk <- housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)

  # an equilibrium state made before severe immunity was added
  ss_create <- equilibrium_ss_create
  local_mocked_bindings(equilibrium_ss_create = function(eqInit) {
    eqSS <- ss_create(eqInit)
    eqSS[setdiff(names(eqSS), c("IVAmat", "IVM_Init"))]
  })
  out <- pipeline(EIR = 10, N = 300, years = 0.5, num_loci = 2, seed = 4L,
                  housekeeping_list = hk)
  expect_gt(out$Loggers$Severe, 0)
  # severe immunity starts at zero and only builds from new infections
  expect_lt(mean(out$Loggers$IVA[out$Loggers$Ages > 20 * 365]), 2)

  # a saved state made before severe immunity was added
  fs <- pipeline(EIR = 10, N = 300, years = 0.5, num_loci = 2, full_save = TRUE,
                 seed = 5L, housekeeping_list = hk)
  new_pop <- c("IVA", "IVM", "IVM_init", "IVA_last_boost_time")
  fs$population_List <- fs$population_List[setdiff(names(fs$population_List), new_pop)]
  new_par <- paste0("g_", c("theta0", "theta1", "dVA", "IV0", "kV", "uVA", "fV0",
                            "aV", "gammaV", "PVM", "dVM",
                            "mean_maternal_severe_immunity",
                            "sum_maternal_severe_immunity"))
  expect_true(all(new_par %in% names(fs$parameters_List)))
  fs$parameters_List <- fs$parameters_List[setdiff(names(fs$parameters_List), new_par)]
  fs$saved_inputs$core_parameter_list[c("theta0", "theta1", "dVA", "IV0", "kV", "uVA",
                                        "fV0", "aV", "gammaV", "PVM", "dVM")] <- NULL
  tf <- tempfile(fileext = ".rds")
  saveRDS(fs, tf)
  re <- pipeline(EIR = 10, N = 300, years = 0.5, num_loci = 2,
                 saved_state_path = tf, seed = 6L, housekeeping_list = hk)
  expect_true(is.finite(re$Loggers$Severe))
  expect_true(all(re$Loggers$IVA >= 0))

})

test_that("maternal immunity at birth survives a save and reload", {
  skip_on_cran()

  hk <- housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)
  fs <- pipeline(EIR = 10, N = 300, years = 0.5, num_loci = 2, full_save = TRUE,
                 seed = 7L, housekeeping_list = hk)
  expect_true(all(is.finite(fs$population_List$ICM_init)))
  expect_gt(max(fs$population_List$ICM_init), 0)
  # ICM = ICM_init * exp(-age / dCM)
  dcm <- fs$parameters_List$g_dCM
  expect_equal(fs$population_List$ICM,
               fs$population_List$ICM_init * exp(-fs$population_List$Ages / dcm),
               tolerance = 1e-6)

  tf <- tempfile(fileext = ".rds")
  saveRDS(fs, tf)
  re <- pipeline(EIR = 10, N = 300, years = 0.5, num_loci = 2, full_save = TRUE,
                 saved_state_path = tf, seed = 8L, housekeeping_list = hk)
  # everyone ages by the run length, so only the people alive in both can be matched
  # by checking the reloaded values are finite and of the same scale as before
  expect_true(all(is.finite(re$population_List$ICM_init)))
  expect_lt(max(re$population_List$ICM_init), 10 * max(fs$population_List$ICM_init))
  expect_true(all(is.finite(re$Loggers$ICM)))
})
