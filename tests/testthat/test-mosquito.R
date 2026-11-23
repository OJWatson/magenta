test_that("mosquito biting intervals are geometric with mean 1/fv", {

  out <- pipeline(EIR = 4, years = 0.1, N = 500, full_save = TRUE, seed = 3L,
                  housekeeping_list = housekeeping_list_create(quiet = TRUE, quiet_test = TRUE,
                                                               cluster = TRUE))
  bites <- out$parameters_List$g_mosquito_next_biting_day_vector

  # never 0 (a mosquito would stall), spread out rather than a fixed 3 day cycle
  expect_gte(min(bites), 1)
  expect_gt(length(unique(bites)), 5)

  # no interventions so fv = 1/3 and the mean is 3 (sd of the mean over 1000 draws ~ 0.08)
  expect_lt(abs(mean(bites) - 3), 0.3)
  expect_lt(abs(mean(bites == 1) - 1/3), 0.06)

})
