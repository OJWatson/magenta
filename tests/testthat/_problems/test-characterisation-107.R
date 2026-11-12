# Extracted from test-characterisation.R:107

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
ints <- suppressMessages(intervention_grab(
  country = "Democratic Republic of the Congo", admin = "Kinshasa",
  year_range = 2012:2013
))
spl <- suppressMessages(spl_grab(
  "Democratic Republic of the Congo", "Kinshasa", 2012:2013
))
expect_snapshot_value(
  list(ints = ints, spl = snapshot_summary(spl)),
  style = "json2", tolerance = 1e-6
)
