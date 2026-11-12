# Extracted from test-characterisation.R:33

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
out <- suppressMessages(pipeline(
  EIR = 2, years = 1, N = 200, seed = 2L,
  update_save = TRUE, update_length = 30,
  ft = 0.5, itn_cov = c(0.2), irs_cov = 0,
  housekeeping_list = quiet_housekeeping()
))
expect_length(out, round(365 / 30))
expect_snapshot_value(snapshot_summary(out), style = "json2", tolerance = 1e-6)
