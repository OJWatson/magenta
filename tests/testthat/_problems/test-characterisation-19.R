# Extracted from test-characterisation.R:19

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
out <- suppressMessages(pipeline(
  EIR = 1, years = 1, N = 200, seed = 1L,
  housekeeping_list = quiet_housekeeping()
))
expect_named(out, c("Loggers"))
expect_type(attr(out, "meta"), "list")
expect_equal(attr(out, "meta")$seed, 1L)
expect_snapshot_value(snapshot_summary(out), style = "json2", tolerance = 1e-6)
