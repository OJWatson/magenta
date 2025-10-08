test_that("pfpr_to_eir_heuristic returns finite EIR for age-restricted targets", {
  skip_on_cran()
  result <- pfpr_to_eir_heuristic(ft = 0.2, PfPR = 0.1, age_range = c(2, 10))
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))
})

test_that("pfpr_to_eir_heuristic validates malformed age ranges", {
  skip_on_cran()
  expect_message(
    pfpr_to_eir_heuristic(PfPR = 0.1, age_range = c(10, 2)),
    "Final prevalence difference"
  )
  expect_error(
    pfpr_to_eir_heuristic(PfPR = 0.1, age_range = c(1)),
    "length 2"
  )
})
