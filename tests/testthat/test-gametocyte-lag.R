test_that("saved states keep the infection state history and cleared strains for the gametocyte lag", {

  out <- pipeline(EIR = 20, years = 0.2, N = 500, ft = 0.4, full_save = TRUE, seed = 3L,
                  housekeeping_list = housekeeping_list_create(quiet = TRUE, quiet_test = TRUE,
                                                               cluster = TRUE))
  pl <- out$populations_event_and_strains_List
  states <- out$population_List$Infection_States

  # the last day simulated, and the oldest day infectivity looks back to from it
  last_day <- out$parameters_List$g_current_time - 1
  oldest <- last_day - ceiling(out$parameters_List$g_delay_gam)

  days <- pl$State_history_days
  hist <- pl$State_history_states

  # everyone has a history that ends in their current state, with a new entry only when the state changes
  expect_true(all(lengths(days) >= 1))
  expect_identical(lengths(days), lengths(hist))
  expect_identical(vapply(hist, function(x) x[length(x)], integer(1)), as.integer(states))
  expect_false(any(vapply(hist, function(x) any(diff(x) == 0), logical(1))))
  expect_false(any(vapply(days, is.unsorted, logical(1))))

  # only the state on the oldest day and any changes since are kept
  expect_true(all(vapply(days, function(d) length(d) == 1 || d[2] > oldest, logical(1))))

  # cleared strains are kept until they are older than the lag, with one barcode each
  cleared <- unlist(pl$Cleared_strain_days)
  expect_gt(length(cleared), 0)
  expect_true(all(cleared > oldest & cleared <= last_day))
  expect_identical(lengths(pl$Cleared_strain_days), lengths(pl$Cleared_strain_barcode_vectors))
  expect_identical(lengths(pl$Cleared_strain_days), lengths(pl$Cleared_strain_day_of_acquisition_vectors))

  # people in S or P can still have recently cleared strains, so are not all uninfectious
  expect_true(any(lengths(pl$Cleared_strain_days)[states %in% c(0, 5)] > 0))

})
