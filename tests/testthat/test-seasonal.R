test_that("match_clean matches loosely and reports fuzzy matches when asked", {
  expect_equal(magenta:::match_clean("kinshasa", c("Bas-Congo", "Kinshasa")), 2)
  expect_equal(magenta:::match_clean("Kinshasa City", c("Bas-Congo", "Kinshasa")), 2)
  expect_output(
    magenta:::match_clean("Kinshasa City", c("Bas-Congo", "Kinshasa"), quiet = FALSE),
    "kinshasa"
  )
})

test_that("admin_match finds admin units by name and country", {
  drc <- "Democratic Republic of the Congo"
  expect_equal(magenta:::admin_match(NULL, NULL), 0)
  expect_message(i <- magenta:::admin_match("Kinshasa", drc), "Returned: Kinshasa")
  ads <- magenta::admin_units_seasonal
  expect_equal(ads$admin1[i], "Kinshasa")
  expect_silent(j <- magenta:::admin_match("Kinshasa", drc, quiet = TRUE))
  expect_equal(i, j)
  # without a country the admin name alone is used
  expect_equal(suppressMessages(magenta:::admin_match("Kinshasa", NULL)), i)
  # fuzzy spellings still land on some admin unit
  k <- suppressMessages(magenta:::admin_match("kinshasa", "Congo Democratic Republic"))
  expect_true(k >= 1 && k <= nrow(ads))
})

test_that("seasonal_profile is flat without an admin unit and periodic with one", {
  flat <- magenta:::seasonal_profile(NULL, NULL)
  expect_equal(flat, rep(1, 365))
  prof <- suppressMessages(
    magenta:::seasonal_profile("Kinshasa", "Democratic Republic of the Congo")
  )
  expect_length(prof, 365)
  expect_true(all(prof >= 0.001))
  expect_gt(max(prof), min(prof))
  expect_equal(mean(prof), 1, tolerance = 0.05)
})

test_that("pfpr_to_eir_heuristic finds an EIR for each target type", {
  eir_micro <- suppressMessages(magenta:::pfpr_to_eir_heuristic(ft = 0.2, PfPR_micro = 0.3))
  eir_pcr <- suppressMessages(magenta:::pfpr_to_eir_heuristic(ft = 0.2, PfPR = 0.3))
  eir_mv <- suppressMessages(magenta:::pfpr_to_eir_heuristic(ft = 0.2, mv = 5))
  expect_true(all(is.finite(c(eir_micro, eir_pcr, eir_mv))))
  expect_gt(eir_micro, eir_pcr)
  eir_all_a <- suppressMessages(
    magenta:::pfpr_to_eir_heuristic(ft = 0.2, PfPR_micro = 0.3, all_a = TRUE)
  )
  expect_lt(eir_all_a, eir_micro)
})
