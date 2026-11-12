quiet <- function() housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)

test_that("save_humans and save_full keep the population after every period", {
  skip_on_cran()
  h <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 31L, num_loci = 2, update_save = TRUE,
    update_length = 180, output = save_humans(), housekeeping_list = quiet()
  ))
  expect_named(h$output[[1]], c(
    "Infection_States", "Zetas", "Ages", "ID", "Strain_infection_state_vectors",
    "Strain_day_of_infection_state_change_vectors", "Strain_barcode_vectors"
  ))
  expect_named(h$output[[2]], "Loggers")

  f <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 31L, num_loci = 2, update_save = TRUE,
    update_length = 180, output = save_full(), housekeeping_list = quiet()
  ))
  expect_true(magenta:::is_full_state(f$output[[1]]))
  expect_true(magenta:::is_full_state(f$output[[2]]))
  expect_named(f$output[[1]]$saved_inputs, c(
    "barcode_list", "spatial_list", "housekeeping_list", "drug_list",
    "vector_adaptation_list", "nmf_list", "core_parameter_list", "spatial_settings"
  ))
  expect_equal(f$meta$saver, "full")

  hf <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 60, seed = 31L, num_loci = 2, update_save = TRUE,
    update_length = 180, output = save_loggers(final = "humans"), housekeeping_list = quiet()
  ))
  expect_true("Zetas" %in% names(hf$output[[2]]))
  expect_false("Loggers" %in% names(hf$output[[2]]))
})

test_that("save_summary keeps COI tables or the raw data frame", {
  skip_on_cran()
  coi <- suppressMessages(pipeline(
    EIR = 2, years = 1, N = 80, seed = 32L, num_loci = 2, update_save = TRUE,
    update_length = 120, output = save_summary(barcodes_every = 1, barcodes_last = 10),
    housekeeping_list = quiet()
  ))
  expect_true(all(
    c("summary_coi", "barcodes", "pcr_prev", "micro_2_10") %in% names(coi$output[[1]])
  ))
  expect_equal(coi$meta$saver, "summary (coi)")
  df <- as.data.frame(coi)
  expect_true(all(c("pcr_prev", "micro_2_10", "overall_treatment_failure") %in% names(df)))
  expect_equal(nrow(df), 3)
  expect_output(print(summary(coi)), "pcr_prev")

  raw <- suppressMessages(pipeline(
    EIR = 2, years = 1, N = 80, seed = 32L, num_loci = 2, update_save = TRUE,
    update_length = 120, output = save_summary(genetics = "df"), housekeeping_list = quiet()
  ))
  expect_s3_class(raw$output[[1]]$df, "data.frame")
  expect_equal(raw$meta$saver, "summary (df)")
})

test_that("set_up_only honours the saver's final choice", {
  skip_on_cran()
  humans <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 50, seed = 33L, num_loci = 2, set_up_only = TRUE,
    output = save_loggers(final = "humans"), housekeeping_list = quiet()
  ))
  expect_true("Zetas" %in% names(humans))
  expect_false("population_List" %in% names(humans))
})

test_that("clear_up = FALSE leaves the pointer for the caller to free", {
  skip_on_cran()
  out <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 50, seed = 34L, num_loci = 2,
    housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE, clear_up = FALSE)
  ))
  expect_type(out$meta$Ptr, "externalptr")
  expect_type(out$output$Ptr, "externalptr")
  state <- magenta:::sim_get_state(out$meta$Ptr)
  expect_equal(state$parameters_List$g_current_time, 366)
  magenta:::sim_free(out$meta$Ptr)
})

test_that("a calendar year range with historic interventions runs", {
  skip_on_cran()
  out <- suppressMessages(pipeline(
    EIR = 1, years = 2012:2013, N = 60, seed = 35L, num_loci = 2,
    country = "Democratic Republic of the Congo", admin = "Kinshasa",
    use_historic_interventions = TRUE, spatial_type = "island",
    update_save = TRUE, update_length = 365, housekeeping_list = quiet()
  ))
  expect_equal(nrow(out$schedule), 2)
  expect_equal(out$params$years, 2)
  expect_gt(out$schedule$itn_cov[1], 0)
  expect_equal(out$output[[2]]$Loggers$Time, 731)
})

test_that("sim_save fails cleanly on the wrong input and sim_restore rejects bad files", {
  expect_error(sim_save(list(a = 1), tempfile()), "must be a magenta_sim")
  expect_error(sim_save(list(), 1), "path must be a single string")
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(list(population_List = 1), tf)
  expect_error(sim_restore(tf), "not a saved magenta state")
  full <- as.list(setNames(rep(1, 5), magenta:::cpp_contract$saved_state))
  saveRDS(full, tf)
  expect_error(sim_restore(tf), "missing saved_inputs")
  newer <- magenta:::new_saved_state(c(full, list(saved_inputs = list())), 99L, "9.9.9")
  saveRDS(newer, tf)
  expect_error(sim_restore(tf), "newer than this version")
})
