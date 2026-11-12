# pop_strains_df() and COI_df_create() on a tiny simulated population.

strains_df <- function(seed = 21L, ibd_length = 1, days = 60) {
  state <- get("tiny_state")(seed = seed, N = 80, ibd_length = ibd_length)
  mu_fv <- get("constant_mu_fv")(365)
  state <- magenta:::run_period(state, seq_len(days), 0.4, mu_fv, seed = seed)
  df <- magenta:::pop_strains_df(state$sim$Ptr,
    seed = seed, sample_size = 0, sample_states = 0:5, ibd = ibd_length > 1
  )
  magenta:::sim_free(state$sim$Ptr, seed)
  df
}

test_that("pop_strains_df returns one row per sampled individual", {
  skip_on_cran()
  df <- strains_df()
  expect_equal(nrow(df), 80)
  expect_true(all(
    c("age", "state", "coi", "coi_detected_pcr", "barcode_states", "nums") %in% names(df)
  ))
  expect_true(any(df$coi > 0))
  infected <- df$coi > 0
  expect_true(all(vapply(df$nums[infected], is.raw, logical(1))))
  expect_true(all(vapply(df$nums[infected], ncol, integer(1)) == 2))

  state <- tiny_state(seed = 22L, N = 80)
  state <- magenta:::run_period(state, 1:60, 0.4, constant_mu_fv(365), seed = 22L)
  some <- magenta:::pop_strains_df(state$sim$Ptr, seed = 22L, sample_size = 10, sample_states = 1:4)
  expect_lte(nrow(some), 10)
  expect_true(all(some$state %in% 1:4))
  magenta:::sim_free(state$sim$Ptr, 22L)
})

test_that("COI_df_create summarises everyone, subsamples, and tables barcodes", {
  skip_on_cran()
  df <- strains_df()
  all <- magenta:::COI_df_create(df, nl = 2)
  expect_named(all, c(
    "summary_coi", "summary_coi_detected", "summary_moi_detected", "summary_polygenom",
    "summary_ages", "summary_unique", "summary_cou", "clonality", "coi_table", "barcode_freq"
  ))
  expect_length(all$barcode_freq, 2)
  expect_true(all(all$barcode_freq >= 0 & all$barcode_freq <= 1))
  expect_equal(sum(all$summary_coi$N), nrow(df))

  with_barcodes <- magenta:::COI_df_create(df, nl = 2, barcodes = TRUE)
  expect_true("barcodes" %in% names(with_barcodes))

  full_stats <- magenta:::COI_df_create(df, nl = 2, mean_only = FALSE)
  expect_true(all(c("sd", "se", "ci") %in% names(full_stats$summary_coi)))

  sub <- magenta:::COI_df_create(df, nl = 2, n = 20)
  expect_equal(sum(sub$summary_coi$N), 20)
  reps <- magenta:::COI_df_create(df, nl = 2, n = 20, reps = 2)
  expect_length(reps, 2)
  expect_equal(sum(reps[[2]]$summary_coi$N), 20)
  more <- magenta:::COI_df_create(df, nl = 2, n = 200, reps = 2)
  expect_equal(sum(more[[1]]$summary_coi$N), 200)
  sizes <- magenta:::COI_df_create(df, nl = 2, n = c(10, 30))
  expect_length(sizes, 3)
  expect_equal(sum(sizes[[1]]$summary_coi$N), 10)
  expect_equal(sum(sizes[[3]]$summary_coi$N), nrow(df))
})

test_that("COI_df_create summarises IBD simulations", {
  skip_on_cran()
  df <- strains_df(seed = 23L, ibd_length = 4)
  expect_true(all(c("pibd", "pibd_within") %in% names(df)))
  res <- magenta:::COI_df_create(df, nl = 2, ibd = TRUE)
  expect_named(res, c(
    "summary_ibd", "summary_within_ibd", "summary_ages",
    "mean_ibd", "mean_ibd_d", "mean_within", "mean_within_d"
  ))
})
