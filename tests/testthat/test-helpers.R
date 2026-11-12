test_that("check helpers describe what is wrong", {
  expect_error(magenta:::check_number("a", "x"), "x must be a single number$")
  expect_error(magenta:::check_number(c(1, 2), "x"), "x must have length 1, not 2")
  expect_error(magenta:::check_number(-1, "x", min = 0), "x must be a single number of at least 0")
  expect_error(magenta:::check_number(2, "x", max = 1), "x must be a single number of at most 1")
  expect_error(magenta:::check_number(1.5, "x", integer = TRUE), "x must be a single whole number")
  expect_error(magenta:::check_number(c(1, NA), "x", n = NULL), "x must be a number vector")
  expect_equal(magenta:::check_number(TRUE, "x", integer = TRUE), TRUE)
  expect_error(magenta:::check_flag(NA, "x"), "x must be TRUE or FALSE")
  expect_error(magenta:::check_flag(c(TRUE, FALSE), "x"), "x must be TRUE or FALSE")
  expect_error(magenta:::check_flag(1, "x", n = NULL), "x must be a logical vector without NAs")
  expect_error(magenta:::check_length(1:3, "x", 2), "x must have length 2, not 3")
  expect_error(
    magenta:::check_class(1, "magenta_drug_list", "x"),
    "must be created by drug_list_create\\(\\)"
  )
  expect_error(magenta:::check_fields(1, magenta:::cpp_contract$nmf, "x"), "x must be a list")
})

test_that("assert helpers and small utilities behave", {
  expect_error(magenta:::assert_scalar_character(1), "must be a scalar character")
  expect_error(magenta:::assert_scalar_logical(1), "must be a scalar logical")
  expect_error(magenta:::assert_scalar_numeric("a"), "must be a scalar numeric")
  expect_silent(magenta:::assert_null_and_func(NULL, magenta:::assert_scalar_numeric))
  expect_error(
    magenta:::assert_null_and_func("a", magenta:::assert_scalar_numeric), "scalar numeric"
  )
  expect_equal(magenta:::ranges(2, 6), list(1:2, 3:4, 5:6))
  df <- magenta:::rbind_list_base(list(list(a = 1, b = "x"), list(a = 2, b = "y")))
  expect_equal(df$a, c(1, 2))
  expect_equal(rownames(df), c("1", "2"))
})

test_that("progress_logging prints markers on the cluster and ticks a bar otherwise", {
  hk_cluster <- housekeeping_list_create(quiet = TRUE, cluster = TRUE)
  res <- vector("list", 50)
  expect_message(
    p <- magenta:::progress_logging(hk_cluster, res, NULL, initial = TRUE), "\\|\\|:\\|"
  )
  expect_equal(p, 2)
  expect_message(p <- magenta:::progress_logging(hk_cluster, res, NULL, i = 2, p_print = 2), "2\\|")
  expect_equal(p, 4)
  expect_silent(p2 <- magenta:::progress_logging(hk_cluster, res, NULL, i = 2, p_print = 10))
  expect_equal(p2, 10)
  expect_message(magenta:::progress_logging(hk_cluster, res, NULL, i = 50, p_print = 98), "98\\|")

  bar <- progress::progress_bar$new(total = 50, show_after = 1e6)
  hk_local <- housekeeping_list_create(quiet = TRUE, cluster = FALSE)
  expect_silent(magenta:::progress_logging(hk_local, res, bar, initial = TRUE))
  expect_silent(magenta:::progress_logging(hk_local, res, bar, i = 1, p_print = 2))
})

test_that("detection helpers give probabilities", {
  mpl <- model_param_list_create()
  fdv <- magenta:::fd(c(0, 365 * 5, 365 * 50), mpl$fD0, mpl$aD, mpl$gammaD)
  expect_equal(fdv[1], mpl$fD0)
  expect_true(all(diff(fdv) > 0) && all(fdv <= 1))
  q <- magenta:::q_fun(mpl$d1, c(0, 5, 50), mpl$ID0, mpl$kD, fdv)
  expect_equal(q[1], 1)
  expect_true(all(q >= mpl$d1 & q <= 1))
})

test_that("magenta_sim prints without timings and with a long call", {
  sim <- magenta:::new_magenta_sim(
    output = list(Loggers = list(Time = 1)),
    params = magenta:::new_magenta_params(list(N = 10, years = 1, EIR = 1)),
    schedule = NULL,
    meta = list(
      call = quote(pipeline(
        EIR = 1, years = 1, N = 10, a_very_long_argument_name = "x", another = "y"
      )),
      seed = 1L, version = "1.4.0", saver = "loggers", final = "loggers", update_save = FALSE
    )
  )
  expect_output(print(sim), "\\.\\.\\.")
  expect_output(print(summary(sim)), "Time")
  empty <- sim
  empty$output <- list(Loggers = list(Ages = 1:3))
  expect_output(print(summary(empty)), "No scalar per-period output")
})

test_that("the assert and check helpers return their input invisibly", {
  expect_identical(magenta:::assert_scalar_character("a"), "a")
  expect_identical(magenta:::assert_scalar_logical(TRUE), TRUE)
  expect_identical(magenta:::assert_scalar_numeric(1), 1)
  dl <- drug_list_create()
  expect_identical(magenta:::check_class(dl, "magenta_drug_list", "x"), dl)
  expect_equal(magenta:::scalar_fields(1), list())
  expect_equal(magenta:::scalar_fields(list(Loggers = list(a = 1, b = 1:2))), list(a = 1))
})
