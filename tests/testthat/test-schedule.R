test_that("build_schedule splits the simulation into periods by first day", {
  s <- build_schedule(years = 1, update_length = 30)
  expect_s3_class(s, "magenta_schedule")
  expect_equal(nrow(s), 12)
  expect_equal(s$day_start, seq(1, 331, by = 30))
  expect_equal(s$day_end, seq(30, 360, by = 30))
  expect_true(all(s$year == 1))
  # scalar columns only
  expect_false(any(vapply(s, is.list, logical(1))))
  expect_null(attr(s, "annual"))

  # a period belongs to the year its first day falls in
  s2 <- build_schedule(years = 2, update_length = 60)
  expect_equal(nrow(s2), 12)
  expect_equal(s2$year, c(rep(1, 7), rep(2, 5)))
  expect_equal(s2$day_end[12], 720)

  # the last period is shortened to the final day
  s3 <- build_schedule(years = 1, update_length = 100)
  expect_equal(nrow(s3), 4)
  expect_equal(s3$day_end, c(100, 200, 300, 365))

  # a single period covering a year
  s4 <- build_schedule(years = 1, update_length = 365)
  expect_equal(nrow(s4), 1)
  expect_equal(s4$day_end, 365)

  # calendar year ranges mean their length
  s5 <- build_schedule(years = 2000:2003, update_length = 365)
  expect_equal(nrow(s5), 4)
  expect_equal(s5$year, 1:4)

  # a part year beyond the last period is simply not simulated
  s6 <- build_schedule(years = 2.5, update_length = 365)
  expect_equal(nrow(s6), 2)
  expect_equal(s6$year, 1:2)
})

test_that("yearly inputs are padded by repeating the last value", {
  s <- build_schedule(
    years = 4, update_length = 365, ft = c(0.2, 0.6), itn_cov = 0.3,
    irs_cov = c(0, 0.1, 0.2), resistance_flag = c(FALSE, TRUE),
    mutation_flag = TRUE
  )
  expect_equal(s$ft, c(0.2, 0.6, 0.6, 0.6))
  expect_equal(s$itn_cov, rep(0.3, 4))
  expect_equal(s$irs_cov, c(0, 0.1, 0.2, 0.2))
  expect_equal(s$resistance_flag, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(s$mutation_flag, rep(TRUE, 4))

  # the first row of each year gives the yearly value back
  m <- build_schedule(years = 2, update_length = 30, ft = c(0.2, 0.6))
  expect_equal(m$ft[!duplicated(m$year)], c(0.2, 0.6))
})

test_that("plaf and importation rows follow the same rule, beside the schedule", {
  plaf <- matrix(c(0.5, 0.5, 0, 0), ncol = 2, byrow = TRUE)
  a <- magenta:::annual_inputs(
    years = 3, plaf = plaf,
    spatial_incidence = c(1, 0), spatial_mosquitoFOI = matrix(c(2, 0), ncol = 1)
  )
  expect_equal(a$plaf, matrix(c(0.5, 0.5, 0, 0, 0, 0), ncol = 2, byrow = TRUE))
  expect_equal(a$spatial_incidence, matrix(c(1, 0, 0), ncol = 1))
  expect_equal(a$spatial_mosquitoFOI, matrix(c(2, 0, 0), ncol = 1))
  expect_equal(magenta:::annual_row(a$spatial_incidence, 2), 0)

  # a plaf vector is one row for every year; no spatial inputs give NULL
  b <- magenta:::annual_inputs(years = 2000:2001, plaf = c(0.1, 0.9))
  expect_equal(b$plaf, matrix(c(0.1, 0.9, 0.1, 0.9), ncol = 2, byrow = TRUE))
  expect_null(b$spatial_incidence)
  expect_null(magenta:::annual_row(b$spatial_incidence, 1))
})

test_that("inputs longer than the simulation are an error", {
  expect_error(
    build_schedule(years = 2, update_length = 365, ft = c(0.1, 0.2, 0.3)),
    "ft has 3 values but the simulation covers 2 year\\(s\\)"
  )
  expect_error(
    magenta:::annual_inputs(years = 1, plaf = matrix(0.5, 2, 2)),
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
    magenta:::annual_inputs(years = 1, plaf = 0.5, spatial_incidence = 1),
    "same dimensions"
  )
  expect_error(magenta:::annual_inputs(years = 1, plaf = 2), "plaf must be a number vector")
  expect_error(build_schedule(years = 1, update_length = 0), "update_length must be")
})

test_that("a simulated year with no period starting in it is an error", {
  expect_error(
    build_schedule(years = 3, update_length = 3 * 365),
    "update_length = 1095 leaves no period starting in year 2"
  )
  expect_error(
    build_schedule(years = 2, update_length = 300),
    "leaves no period starting in year 2"
  )
  # 400-day periods over two years start one in each year
  s <- build_schedule(years = 2, update_length = 400)
  expect_equal(s$year, 1:2)
})

test_that("the schedule prints its period table", {
  s <- build_schedule(years = 1, update_length = 365)
  expect_output(print(s), "1 period\\(s\\) over 365 days")
  expect_output(print(s), "drug_choice")
  # a column subset keeps the class and prints without the day count
  sub <- s[, c("period", "year", "ft")]
  expect_s3_class(sub, "magenta_schedule")
  expect_output(print(sub), "1 period\\(s\\)\\n")
  expect_output(print(sub), "year")
})

test_that("the padding helpers reject empty and non-numeric inputs", {
  expect_error(magenta:::pad_annual(numeric(0), 2, "ft"), "ft must have at least one value")
  expect_error(magenta:::pad_annual_rows("a", 2, "plaf"), "plaf must be numeric")
  expect_null(magenta:::pad_annual_rows(NULL, 2, "plaf"))
})
