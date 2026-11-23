test_that("mosquito biting intervals are geometric with mean 1/fv", {

  out <- pipeline(EIR = 4, years = 0.1, N = 500, full_save = TRUE, seed = 3L,
                  housekeeping_list = housekeeping_list_create(quiet = TRUE, quiet_test = TRUE,
                                                               cluster = TRUE))
  bites <- out$parameters_List$g_mosquito_next_biting_day_vector

  # never 0 (a mosquito would stall), and spread out rather than a fixed cycle
  expect_gte(min(bites), 1)
  expect_gt(length(unique(bites)), 5)

  # with no interventions the biting rate is 1 / (tau1 + tau2)
  mpl <- model_param_list_create()
  fv <- 1 / (mpl$tau1 + mpl$tau2)

  # geometric on 1, 2, ...: mean 1/fv, variance (1 - fv)/fv^2 and P(1) = fv,
  # checked to within 5 standard errors
  n <- length(bites)
  expect_lt(abs(mean(bites) - 1 / fv), 5 * sqrt((1 - fv) / fv^2 / n))
  expect_lt(abs(mean(bites == 1) - fv), 5 * sqrt(fv * (1 - fv) / n))

})
