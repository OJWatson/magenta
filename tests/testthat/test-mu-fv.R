# The closed-form mosquito rates against the deterministic model they come
# from. The model is run here (slowly) and its mu and fv compared with
# mu_fv_from_coverage() over every configuration the pipeline can produce.

eq_for <- function(EIR = 5, ft = 0.4, country = NULL, admin = NULL) {
  equilibrium_init_create(
    age_vector = magenta:::age_brackets(100, 40, TRUE), het_brackets = 5,
    ft = ft, EIR = EIR, model_param_list = model_param_list_create(),
    country = country, admin = admin, quiet = TRUE
  )
}

expect_same_rates <- function(eq, ft, itn_cov, irs_cov, years, int_times = NULL) {
  odin <- suppressMessages(magenta:::run_odin_model(eq, ft, itn_cov, irs_cov, years, int_times))
  direct <- magenta:::mu_fv_from_coverage(eq, ft, itn_cov, irs_cov, years, int_times)
  testthat::expect_equal(nrow(direct), round(years * 365))
  testthat::expect_identical(direct$mu, odin$mu)
  testthat::expect_identical(direct$fv, odin$fv)
  invisible(direct)
}

test_that("the closed-form rates equal the deterministic model's: single interventions", {
  skip_on_cran()
  eq <- eq_for()
  itn <- expect_same_rates(eq, ft = 0.4, itn_cov = 0.5, irs_cov = 0, years = 1)
  expect_gt(itn$mu[365], 0.132)
  expect_lt(itn$fv[365], 1 / 3)
  expect_same_rates(eq, ft = 0.4, itn_cov = 0, irs_cov = 0.3, years = 1)
  expect_same_rates(eq, ft = 0.4, itn_cov = 0.5, irs_cov = 0.3, years = 1)
})

test_that("the closed-form rates equal the deterministic model's: ramps and seasons", {
  skip_on_cran()
  eq <- eq_for()
  # a multi-year ramp in both coverages, with the last year padded
  expect_same_rates(eq,
    ft = seq(0.2, 0.6, by = 0.1), itn_cov = seq(0.1, 0.5, by = 0.1),
    irs_cov = c(0, 0.1, 0.2, 0.3), years = 5
  )
  # a part year
  expect_same_rates(eq, ft = 0.4, itn_cov = c(0.4, 0.6), irs_cov = 0.1, years = 1.5)
  # coverage changing at other times than the anniversaries
  expect_same_rates(eq,
    ft = c(0.4, 0.4), itn_cov = c(0.2, 0.7), irs_cov = c(0, 0.5), years = 2,
    int_times = c(-25, 0, 100)
  )
  # the seasonal Kinshasa run with its historic interventions
  cc <- "Democratic Republic of the Congo"
  ints <- suppressMessages(intervention_grab(cc, "Kinshasa", 2012:2013))
  eq_k <- eq_for(EIR = 1, ft = ints$ft[1], country = cc, admin = "Kinshasa")
  expect_same_rates(eq_k, ft = ints$ft, itn_cov = ints$itn_cov, irs_cov = ints$irs_cov, years = 2)
})

test_that("mu_fv_create uses the closed form, the shortcut, or the full model", {
  skip_on_cran()
  eq <- eq_for()
  direct <- magenta:::mu_fv_create(eq, ft = 0.4, itn_cov = 0.5, irs_cov = 0.2, years = 1)
  expect_identical(direct, magenta:::mu_fv_from_coverage(eq, 0.4, 0.5, 0.2, 1))
  none <- magenta:::mu_fv_create(NULL, ft = 0.4, itn_cov = 0, irs_cov = 0, years = 1)
  expect_equal(unique(none$mu), 0.132)
  full <- suppressMessages(magenta:::mu_fv_create(eq,
    ft = 0.4, itn_cov = 0.5, irs_cov = 0.2, years = 1, full = TRUE
  ))
  expect_true(all(c("t", "mu", "fv", "Sout", "prev") %in% names(full)))
  expect_identical(full$mu, direct$mu)

  det <- suppressMessages(run_deterministic_model(
    EIR = 5, ft = 0.4, years = 1, itn_cov = 0.5, irs_cov = 0.2, quiet = TRUE
  ))
  expect_identical(det$mu, full$mu)
  expect_error(
    magenta:::mu_fv_from_coverage(eq, 0.4, 0.5, 0, 1, int_times = c(-25, 0, 365)),
    "one entry per yearly value"
  )
  expect_error(
    magenta:::mu_fv_from_coverage(eq, 0.4, 0.5, 0, 1, int_times = c(10, 20)),
    "on or before day 1"
  )
})
