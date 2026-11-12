test_that("build_schedule splits the simulation into periods by first day", {
  s <- build_schedule(years = 1, update_length = 30, plaf = rep(0.5, 2))
  expect_s3_class(s, "magenta_schedule")
  expect_equal(nrow(s), 12)
  expect_equal(s$day_start, seq(1, 331, by = 30))
  expect_equal(s$day_end, seq(30, 360, by = 30))
  expect_true(all(s$year == 1))

  # a period belongs to the year its first day falls in
  s2 <- build_schedule(years = 2, update_length = 60, plaf = rep(0.5, 2))
  expect_equal(nrow(s2), 12)
  expect_equal(s2$year, c(rep(1, 7), rep(2, 5)))
  expect_equal(s2$day_end[12], 720)

  # the last period is shortened to the final day
  s3 <- build_schedule(years = 1, update_length = 100, plaf = rep(0.5, 2))
  expect_equal(nrow(s3), 4)
  expect_equal(s3$day_end, c(100, 200, 300, 365))

  # a single period covering everything
  s4 <- build_schedule(years = 3, update_length = 3 * 365, plaf = rep(0.5, 2))
  expect_equal(nrow(s4), 1)
  expect_equal(s4$day_end, 1095)

  # calendar year ranges mean their length
  s5 <- build_schedule(years = 2000:2003, update_length = 365, plaf = rep(0.5, 2))
  expect_equal(nrow(s5), 4)
  expect_equal(s5$year, 1:4)
})

test_that("yearly inputs are padded by repeating the last value", {
  s <- build_schedule(
    years = 4, update_length = 365, ft = c(0.2, 0.6), itn_cov = 0.3,
    irs_cov = c(0, 0.1, 0.2), resistance_flag = c(FALSE, TRUE),
    mutation_flag = TRUE, plaf = rep(0.5, 2)
  )
  expect_equal(s$ft, c(0.2, 0.6, 0.6, 0.6))
  expect_equal(s$itn_cov, rep(0.3, 4))
  expect_equal(s$irs_cov, c(0, 0.1, 0.2, 0.2))
  expect_equal(s$resistance_flag, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(s$mutation_flag, rep(TRUE, 4))
  annual <- attr(s, "annual")
  expect_equal(annual$ft, c(0.2, 0.6, 0.6, 0.6))

  # plaf rows and importation rows follow the same rule
  plaf <- matrix(c(0.5, 0.5, 0, 0), ncol = 2, byrow = TRUE)
  s2 <- build_schedule(
    years = 3, update_length = 365, plaf = plaf,
    spatial_incidence = c(1, 0), spatial_mosquitoFOI = matrix(c(2, 0), ncol = 1)
  )
  expect_equal(s2$plaf, list(c(0.5, 0.5), c(0, 0), c(0, 0)))
  expect_equal(s2$importation_human, list(1, 0, 0))
  expect_equal(s2$importation_mosquito, list(2, 0, 0))

  # a plaf vector is one row for every year
  s3 <- build_schedule(years = 2, update_length = 365, plaf = c(0.1, 0.9))
  expect_equal(s3$plaf, list(c(0.1, 0.9), c(0.1, 0.9)))

  # no spatial inputs give NULL rows
  expect_null(s3$importation_human[[1]])
})

test_that("inputs longer than the simulation are an error", {
  expect_error(
    build_schedule(years = 2, update_length = 365, ft = c(0.1, 0.2, 0.3)),
    "ft has 3 values but the simulation covers 2 year\\(s\\)"
  )
  expect_error(
    build_schedule(years = 1, update_length = 365, plaf = matrix(0.5, 2, 2)),
    "plaf has 2 rows but the simulation covers 1 year\\(s\\)"
  )
  expect_error(
    build_schedule(years = 1, update_length = 365, ft = 1.5),
    "ft must be a number vector between 0 and 1"
  )
  expect_error(
    build_schedule(years = 1, update_length = 365, mutation_flag = 1),
    "mutation_flag must be a logical vector"
  )
  expect_error(
    build_schedule(years = 1, update_length = 365, spatial_incidence = 1),
    "same dimensions"
  )
  expect_error(build_schedule(years = 1, update_length = 0), "update_length must be")
})

test_that("the schedule prints its period table", {
  s <- build_schedule(years = 1, update_length = 365, plaf = 0.5)
  expect_output(print(s), "1 period\\(s\\) over 365 days")
})

test_that("the padding helpers reject empty and non-numeric inputs", {
  expect_error(magenta:::pad_annual(numeric(0), 2, "ft"), "ft must have at least one value")
  expect_error(magenta:::pad_annual_rows("a", 2, "plaf"), "plaf must be numeric")
  expect_null(magenta:::pad_annual_rows(NULL, 2, "plaf"))
})
