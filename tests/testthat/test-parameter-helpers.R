test_that("model_param_list_create honours the half-life arguments", {
  mpl <- model_param_list_create(irs_half_life = 100, itn_half_life = 200)
  expect_equal(mpl$irs_half_life, 100)
  expect_equal(mpl$itn_half_life, 200)
  expect_equal(mpl$irs_loss, log(2) / 100)
  expect_equal(mpl$itn_loss, log(2) / 200)

  default <- model_param_list_create()
  expect_equal(default$irs_half_life, 0.5 * 365)
  expect_equal(default$itn_half_life, 2.64 * 365)
})

test_that("model_param_list_create returns finite scalars", {
  mpl <- model_param_list_create()
  lens <- vapply(mpl, length, integer(1))
  expect_true(all(lens[names(lens) != "cov" & names(lens) != "pop_split"] == 1))
  expect_true(all(is.finite(unlist(mpl))))
})

test_that("microscopy_detected uses the supplied age bracket", {
  ages <- c(500, 1000, 5000)
  ids <- c(0, 0, 0)
  states <- c(1, 2, 0)

  # only the second individual is in the 2-10 year bracket; they are in state A
  # with zero detection immunity, so detection probability is 1
  expect_no_warning(res <- magenta:::microscopy_detected(ages, ids, states))
  expect_equal(res, 1)

  # a bracket that includes everyone: one clinical (state 1) and one A out of 3
  res_all <- magenta:::microscopy_detected(ages, ids, states, age_brackets = c(0, 1e5))
  expect_equal(res_all, 1 / 3 + 1 / 3)
})

test_that("pfpr_to_eir_heuristic accepts an age range", {
  eir_all <- suppressMessages(magenta:::pfpr_to_eir_heuristic(ft = 0, PfPR_micro = 0.2))
  eir_2_10 <- suppressMessages(magenta:::pfpr_to_eir_heuristic(
    ft = 0, PfPR_micro = 0.2, age_range = c(2, 10)
  ))
  expect_true(is.finite(eir_all) && eir_all > 0)
  expect_true(is.finite(eir_2_10) && eir_2_10 > 0)
  expect_false(isTRUE(all.equal(eir_all, eir_2_10)))
})

test_that("metapopulation spatial types give a clear error", {
  expect_error(
    magenta:::spl_create(spatial_type = 2),
    "metapopulation spatial lists are not implemented"
  )
  expect_error(
    suppressMessages(pipeline(
      EIR = 1, years = 1, N = 50, seed = 1L, spatial_type = "metapop",
      housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
    )),
    "metapopulation simulations are not implemented"
  )
})

test_that("intervention_grab returns per-year vectors without touching options()", {
  old <- options("quiet_admin")
  on.exit(options(old), add = TRUE)
  ints <- suppressMessages(intervention_grab(
    country = "Democratic Republic of the Congo", admin = "Kinshasa",
    year_range = 2010:2012
  ))
  expect_named(ints, c("itn_cov", "irs_cov", "ft"))
  expect_length(ints$itn_cov, 3)
  expect_length(ints$ft, 3)
  expect_null(getOption("quiet_admin"))
})
