# The five bindings over the C++ entry points, driven directly.

test_that("sim_init returns a state pointer and loggers for a tiny population", {
  skip_on_cran()
  state <- tiny_state(seed = 11L, N = 40)
  expect_type(state$sim$Ptr, "externalptr")
  expect_type(state$sim$Loggers, "list")
  expect_length(state$sim$Loggers$Ages, 40)
  expect_null(state$sim$Loggers$Time)
  magenta:::sim_free(state$sim$Ptr, 11L)
})

test_that("sim_update advances the clock by the days it is given", {
  skip_on_cran()
  state <- tiny_state(seed = 12L, N = 40)
  mu_fv <- constant_mu_fv(365)
  # the loggers report the day after the last simulated one
  state <- magenta:::run_period(state, 1:30, 0.4, mu_fv, seed = 12L)
  expect_equal(state$sim$Loggers$Time, 31)
  expect_equal(state$sim$Loggers$Log_Counter, 30)
  state <- magenta:::run_period(state, 31:45, 0.4, mu_fv, seed = 12L)
  expect_equal(state$sim$Loggers$Time, 46)
  magenta:::sim_free(state$sim$Ptr, 12L)
})

test_that("sim_get_state returns the five lists and sim_init_saved restarts from them", {
  skip_on_cran()
  state <- tiny_state(seed = 13L, N = 40)
  state <- magenta:::run_period(state, 1:20, 0.4, constant_mu_fv(365), seed = 13L)
  saved <- magenta:::sim_get_state(state$sim$Ptr, 13L)
  expect_named(saved, magenta:::cpp_contract$saved_state)
  expect_equal(saved$parameters_List$g_current_time, 21)
  expect_length(saved$population_List$Infection_States, 40)
  magenta:::sim_free(state$sim$Ptr, 13L)

  restored <- magenta:::sim_init_saved(saved, 13L)
  expect_type(restored$Ptr, "externalptr")
  again <- magenta:::sim_get_state(restored$Ptr, 13L)
  expect_equal(again$population_List$Infection_States, saved$population_List$Infection_States)
  expect_equal(again$parameters_List$g_current_time, 21)
  expect_invisible(magenta:::sim_free(restored$Ptr, 13L))
})

test_that("the bindings check their inputs", {
  expect_error(magenta:::sim_get_state("ptr", 1L), "statePtr must be an external pointer")
  expect_error(magenta:::sim_free(NULL, 1L), "statePtr must be an external pointer")
  expect_error(magenta:::sim_init_saved("x", 1L), "saved_state must be a list")
  expect_error(
    magenta:::sim_init_saved(list(population_List = 1), 1L),
    "missing the elements: populations_event_and_strains_List"
  )
})

test_that("simulation_R still dispatches, with a deprecation warning", {
  skip_on_cran()
  state <- tiny_state(seed = 14L, N = 40)
  ptr <- state$sim$Ptr
  expect_warning(
    got <- magenta:::simulation_R(list(statePtr = ptr), 14L),
    "simulation_R\\(\\) is deprecated"
  )
  expect_named(got, magenta:::cpp_contract$saved_state)

  pl <- magenta:::param_list_simulation_update_create(
    years = 10 / 365, ft = 0.4, statePtr = ptr,
    spatial_list = state$spatial_list, drug_list = state$drug_list,
    barcode_list = state$barcode_list
  )
  expect_warning(upd <- magenta:::simulation_R(pl, 14L), "deprecated")
  expect_equal(upd$Loggers$Time, 11)

  expect_warning(magenta:::simulation_R(list(statePtr = ptr, finalizer = TRUE), 14L), "deprecated")
  expect_warning(
    expect_error(magenta:::simulation_R(list(foo = 1), 14L), "No matching parameter list"),
    "deprecated"
  )
})

test_that("param_list_simulation_update_create fills in constant rates and checks lengths", {
  skip_on_cran()
  state <- tiny_state(seed = 15L, N = 40)
  pl <- magenta:::param_list_simulation_update_create(
    years = 10 / 365, statePtr = state$sim$Ptr,
    spatial_list = state$spatial_list, drug_list = state$drug_list,
    barcode_list = state$barcode_list
  )
  expect_identical(names(pl), magenta:::cpp_contract$update)
  expect_equal(pl$mu_vec, rep(0.132, 10))
  expect_equal(pl$fv_vec, rep(1 / 3, 10))
  expect_error(
    magenta:::param_list_simulation_update_create(
      years = 10 / 365, statePtr = state$sim$Ptr, mu_vec = rep(0.1, 3),
      spatial_list = state$spatial_list, drug_list = state$drug_list,
      barcode_list = state$barcode_list
    ),
    "mu_vec must have length 10, not 3"
  )
  magenta:::sim_free(state$sim$Ptr, 15L)
})
