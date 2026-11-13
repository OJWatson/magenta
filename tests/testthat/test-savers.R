quiet <- function() housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)

test_that("savers are classed objects with a name and a final choice", {
  s <- save_loggers()
  expect_s3_class(s, "magenta_saver")
  expect_equal(s$final, "loggers")
  expect_output(print(s), "loggers, final = loggers")
  expect_equal(save_full()$final, "full")
  expect_equal(save_summary(final = "humans")$final, "humans")
  expect_error(save_summary(sample_reps = 0), "sample_reps must be")
  expect_error(magenta:::new_saver(function(x) x, "bad"), "function\\(state, period, params\\)")
})

test_that("the deprecated flags map onto the same saver with a warning", {
  flags <- function(...) {
    f <- list(
      full_save = FALSE, full_update_save = FALSE, human_only_full_save = FALSE,
      human_update_save = FALSE, genetics_df_without_summarising = FALSE,
      summary_saves_only = FALSE, only_allele_freqs = TRUE, mean_only = TRUE,
      save_lineages = FALSE, sample_size = Inf, sample_states = 0:5,
      age_breaks = c(-0.001, 5, 15, 100.1), sample_reps = 1, update_save_func = NULL
    )
    utils::modifyList(f, list(...))
  }
  m <- magenta:::saver_from_flags(flags(full_save = TRUE))
  expect_equal(m$saver$name, "loggers")
  expect_equal(m$saver$final, "full")
  expect_equal(m$replacement, 'save_loggers(final = "full")')

  m <- magenta:::saver_from_flags(flags(human_update_save = TRUE, summary_saves_only = TRUE))
  expect_equal(m$saver$name, "summary (coi)")
  m <- magenta:::saver_from_flags(flags(
    human_update_save = TRUE, summary_saves_only = TRUE,
    genetics_df_without_summarising = TRUE, save_lineages = TRUE
  ))
  expect_equal(m$saver$name, "summary (allele_freqs)")
  expect_match(m$replacement, "lineages = TRUE")
  m <- magenta:::saver_from_flags(flags(
    human_update_save = TRUE, summary_saves_only = TRUE,
    genetics_df_without_summarising = TRUE, only_allele_freqs = FALSE
  ))
  expect_equal(m$saver$name, "summary (df)")
  m <- magenta:::saver_from_flags(flags(human_update_save = TRUE, human_only_full_save = TRUE))
  expect_equal(m$saver$name, "humans")
  expect_equal(m$saver$final, "humans")
  m <- magenta:::saver_from_flags(flags(human_update_save = TRUE, full_update_save = TRUE))
  expect_equal(m$saver$name, "full")

  expect_error(
    magenta:::resolve_saver(save_loggers(), flags(full_save = TRUE), c(full_save = TRUE)),
    "deprecated argument\\(s\\) full_save"
  )
  expect_error(magenta:::resolve_saver("loggers", flags(), c(full_save = FALSE)), "output must be")
})

test_that("pipeline warns about the old flags and gives the same output", {
  skip_on_cran()
  new <- suppressMessages(pipeline(
    EIR = 2, years = 1, N = 100, seed = 7L, update_save = TRUE, update_length = 90,
    output = save_summary(genetics = "allele_freqs"), housekeeping_list = quiet()
  ))
  expect_warning(
    old <- suppressMessages(pipeline(
      EIR = 2, years = 1, N = 100, seed = 7L, update_save = TRUE, update_length = 90,
      human_update_save = TRUE, summary_saves_only = TRUE,
      genetics_df_without_summarising = TRUE, housekeeping_list = quiet()
    )),
    "deprecated; use output = save_summary\\("
  )
  expect_equal(sim_output(old), sim_output(new))
  expect_equal(sim_meta(new)$saver, "summary (allele_freqs)")
})

test_that("a user function is accepted as the saver", {
  skip_on_cran()
  mine <- function(state, period, params) {
    list(period = period$i, day_end = period$day_end, time = state$sim$Loggers$Time)
  }
  out <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 8L, update_save = TRUE, update_length = 90,
    output = mine, housekeeping_list = quiet()
  ))
  expect_equal(out[[1]], list(period = 1, day_end = 90, time = 91))
  expect_equal(out[[3]]$period, 3)
  expect_named(out[[4]], "Loggers")
  expect_equal(sim_meta(out)$saver, "user function")

  # the old update_save_func signature still works, with a warning
  legacy <- function(res, i, sim.out, mpl, num_loci) {
    res[[i]] <- list(i = i, time = sim.out$Loggers$Time)
    res
  }
  expect_warning(
    out2 <- suppressMessages(pipeline(
      EIR = 1, years = 1, N = 100, seed = 8L, update_save = TRUE, update_length = 90,
      update_save_func = legacy, housekeeping_list = quiet()
    )),
    "update_save_func"
  )
  expect_equal(out2[[2]], list(i = 2, time = 181))
})

test_that("magenta_sim prints, summarises and converts to a data frame", {
  skip_on_cran()
  out <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 9L, update_save = TRUE, update_length = 90,
    housekeeping_list = quiet()
  ))
  expect_s3_class(out, "magenta_sim")
  # the 1.3.x shape: the periods, with the run's inputs as attributes
  expect_length(out, 4)
  expect_null(names(out))
  expect_identical(attr(out, "meta"), sim_meta(out))
  expect_identical(attr(out, "params"), sim_params(out))
  expect_null(attributes(sim_output(out)))
  expect_s3_class(sim_params(out), "magenta_params")
  expect_s3_class(sim_schedule(out), "magenta_schedule")
  expect_equal(sim_params(out)$N, 100)
  expect_equal(sim_params(out)$seed, 9L)
  expect_null(sim_meta(out)$Ptr)
  expect_output(print(out), "N: 100   years: 1   EIR: 1")
  expect_output(print(out), "loggers, 4 period\\(s\\), final = loggers")
  expect_output(print(sim_params(out)), "<magenta_params>")

  df <- as.data.frame(out)
  expect_equal(nrow(df), 4)
  expect_true(all(
    c("period", "day_end", "year", "Time", "S", "D", "A", "U", "T", "P", "Mos_I") %in% names(df)
  ))
  expect_equal(df$Time, c(91, 181, 271, 361))
  expect_output(print(summary(out)), "Per-period output \\(4 period\\(s\\)\\)")

  single <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 9L, housekeeping_list = quiet()
  ))
  expect_named(single, "Loggers")
  expect_named(sim_output(single), "Loggers")
  expect_equal(nrow(as.data.frame(single)), 1)
  expect_output(print(single), "1 period\\(s\\)")
})

test_that("sim_save writes a versioned state that sim_restore reads back", {
  skip_on_cran()
  out <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 10L, num_loci = 3,
    output = save_loggers(final = "full"), housekeeping_list = quiet()
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  expect_invisible(sim_save(out, path))
  saved <- sim_restore(path)
  expect_s3_class(saved, "magenta_saved_state")
  expect_equal(saved$format_version, 1L)
  expect_equal(saved$version, as.character(utils::packageVersion("magenta")))
  expect_true(all(magenta:::cpp_contract$saved_state %in% names(saved)))
  expect_named(saved$saved_inputs, c(
    "barcode_list", "spatial_list", "housekeeping_list", "drug_list",
    "vector_adaptation_list", "nmf_list", "core_parameter_list", "spatial_settings"
  ))
  expect_s3_class(saved$params, "magenta_params")
  expect_s3_class(saved$schedule, "magenta_schedule")
  expect_equal(saved$meta$seed, 10L)
  expect_output(print(saved), "format 1")

  # continuing from it uses the saved drug list and runs on
  cont <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 11L, num_loci = 3,
    saved_state_path = path, update_save = TRUE, update_length = 180,
    housekeeping_list = quiet()
  ))
  expect_gt(cont[[1]]$Time, 365)

  # a run without a full final state cannot be saved
  loggers_only <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 10L, housekeeping_list = quiet()
  ))
  expect_error(sim_save(loggers_only, path), "did not keep its final state")
  expect_error(sim_restore(tempfile()), "no file at")

  # a set_up_only state can be saved too
  state <- suppressMessages(pipeline(
    EIR = 1, years = 1, N = 100, seed = 10L, num_loci = 3, set_up_only = TRUE,
    housekeeping_list = quiet()
  ))
  sim_save(state, path)
  expect_equal(sim_restore(path)$format_version, 1L)

  # older plain-list saves are accepted and upgraded on read
  saveRDS(unclass(state), path)
  old <- sim_restore(path)
  expect_equal(old$format_version, 0L)
  expect_error(sim_restore(system.file("DESCRIPTION", package = "magenta")))
})
