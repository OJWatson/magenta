# Extracted from test-characterisation.R:136

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
state <- suppressMessages(pipeline(
  EIR = 1, years = 1, N = 100, seed = 6L, num_loci = 3,
  set_up_only = TRUE, full_save = TRUE,
  housekeeping_list = quiet_housekeeping()
))
expect_true(all(c(
  "population_List", "populations_event_and_strains_List", "scourge_List",
  "parameters_List", "parameter_housekeeping_List", "saved_inputs"
) %in% names(state)))
expect_snapshot_value(snapshot_summary(state), style = "json2", tolerance = 1e-6)
