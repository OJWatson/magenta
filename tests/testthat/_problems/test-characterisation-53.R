# Extracted from test-characterisation.R:53

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
out <- suppressMessages(pipeline(
  EIR = 5, years = 1, N = 300, seed = 3L, num_loci = 6,
  update_save = TRUE, update_length = 60,
  human_update_save = TRUE, summary_saves_only = TRUE,
  sample_size = Inf, sample_reps = 1, mean_only = TRUE,
  housekeeping_list = quiet_housekeeping()
))
first <- out[[1]]
expect_true(all(c(
  "summary_coi", "summary_coi_detected", "summary_polygenom",
  "summary_unique", "summary_cou", "clonality", "barcode_freq",
  "overall_treatment_failure", "pcr_prev", "micro_2_10"
) %in% names(first)))
expect_snapshot_value(snapshot_summary(out), style = "json2", tolerance = 1e-6)
