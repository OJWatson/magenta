# Findings from R_CODE_REVIEW.md section 4.6 whose fix would change results.
#
# Each test pins what the code does today so that the oddity is visible and
# cannot change by accident. The comment above each test says what the
# obvious fix would make it expect instead. None of these fixes has been
# applied: they are waiting on a decision (see the review's section 7).

test_that("finding: every C++ call reseeds the generator, so updates replay the same stream", {
  skip_on_cran()
  # Two identical simulations updated for the same ten days give identical
  # loggers even though the generator was moved on between them, because
  # sim_update() calls set.seed(seed) itself. With seeding done once at set up
  # (review 6.6) the two would differ, and every update_save snapshot would
  # move.
  a <- tiny_state(seed = 3L)
  b <- tiny_state(seed = 3L)
  mu_fv <- constant_mu_fv(365)
  a <- magenta:::run_period(a, 1:10, 0.4, mu_fv, seed = 3L)
  runif(1000)
  b <- magenta:::run_period(b, 1:10, 0.4, mu_fv, seed = 3L)
  expect_identical(a$sim$Loggers, b$sim$Loggers)
  magenta:::sim_free(a$sim$Ptr, 3L)
  magenta:::sim_free(b$sim$Ptr, 3L)
})

test_that("finding: a single-period run uses the smallest non-zero ft for every year", {
  skip_on_cran()
  # With update_save = FALSE a two year run with ft = c(0.2, 0.8) is simulated
  # at ft = 0.2 throughout (the deterministic mosquito model does see both
  # years). An all-zero ft becomes Inf with a warning. The obvious fix, one
  # period per year when update_save = FALSE, would change multi-year
  # single-period runs; the one-year characterisation runs would not move.
  out <- suppressMessages(pipeline(
    EIR = 1, years = 2, N = 60, seed = 4L, ft = c(0.2, 0.8),
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(nrow(out$schedule), 1)
  expect_equal(out$schedule$ft, 0.2)
  expect_equal(attr(out$schedule, "annual")$ft, c(0.2, 0.8))

  expect_warning(
    out0 <- suppressMessages(pipeline(
      EIR = 1, years = 1, N = 60, seed = 4L, ft = 0,
      housekeeping_list = quiet_housekeeping_test()
    )),
    "no non-missing arguments to min"
  )
  expect_equal(out0$schedule$ft, Inf)
})

test_that("finding: a resumed run takes the barcode and spatial lists from the new call", {
  skip_on_cran()
  # saved_inputs carries all seven lists, but on resume only the drug,
  # housekeeping, nmf, vector adaptation and core lists are used; the barcode
  # and spatial lists are rebuilt from the resuming call's arguments. The fix
  # (use the saved ones) only changes resumed runs whose arguments differ from
  # the saved ones.
  first <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 5L, num_loci = 3, plaf = rep(0.1, 3),
    output = save_loggers(final = "full"),
    housekeeping_list = quiet_housekeeping_test()
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  sim_save(first, path)

  resumed <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 6L, num_loci = 3, plaf = rep(0.9, 3),
    saved_state_path = path,
    drug_list = drug_list_create(drug_choice = 0, number_of_drugs = 1,
      drugs = list(drug_create_al()), number_of_resistance_loci = 3,
      cost_of_resistance = rep(1, 3)
    ),
    housekeeping_list = quiet_housekeeping_test()
  ))
  expect_equal(resumed$params$barcode$plaf, rep(0.9, 3))
  expect_equal(resumed$params$drug, first$params$drug)
  expect_identical(resumed$params$drug$drugs[[1]], drug_create_default_no_resistance())
})

test_that("misspelt pipeline() arguments are reported as such", {
  expect_error(
    suppressMessages(pipeline(EIR = 1, years = 1, N = 50, update_lenght = 30)),
    "unknown argument\\(s\\) to pipeline\\(\\): update_lenght"
  )
})

test_that("finding: spl_grab drops the year 2000 from a range that ends there", {
  # The historic series runs 2000-2017 but the truncation branch only applies
  # to end years after 2000, so a range ending in 2000 keeps no data rows: a
  # range starting before 2000 comes back one row short (all zeros), and
  # 2000:2000 is an error. Changing `end_year > 2000` to `>= 2000` gives
  # six rows for 1995:2000 (the last one the 2000 data) and one row for
  # 2000:2000; nothing else moves.
  cc <- "Democratic Republic of the Congo"
  short <- suppressMessages(spl_grab(cc, "Kinshasa", 1995:2000))
  expect_equal(nrow(short$incidence), 5)
  expect_true(all(short$incidence == 0))
  expect_error(suppressMessages(spl_grab(cc, "Kinshasa", 2000:2000)), "argument of length 0")
})
