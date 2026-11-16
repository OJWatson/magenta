# The five bindings over the C++ entry points, driven directly.

test_that("sim_init returns a state pointer and loggers for a tiny population", {
  skip_on_cran()
  state <- tiny_state(seed = 11L, N = 40)
  expect_type(state$sim$Ptr, "externalptr")
  expect_type(state$sim$Loggers, "list")
  expect_length(state$sim$Loggers$Ages, 40)
  expect_null(state$sim$Loggers$Time)
  magenta:::sim_free(state$sim$Ptr)
})

test_that("sim_update advances the clock by the days it is given", {
  skip_on_cran()
  state <- tiny_state(seed = 12L, N = 40)
  mu_fv <- constant_mu_fv(365)
  # the loggers report the day after the last simulated one
  state <- magenta:::run_period(state, 1:30, 0.4, mu_fv)
  expect_equal(state$sim$Loggers$Time, 31)
  expect_equal(state$sim$Loggers$Log_Counter, 30)
  state <- magenta:::run_period(state, 31:45, 0.4, mu_fv)
  expect_equal(state$sim$Loggers$Time, 46)
  magenta:::sim_free(state$sim$Ptr)
})

test_that("sim_get_state returns the five lists and sim_init_saved restarts from them", {
  skip_on_cran()
  state <- tiny_state(seed = 13L, N = 40)
  state <- magenta:::run_period(state, 1:20, 0.4, constant_mu_fv(365))
  saved <- magenta:::sim_get_state(state$sim$Ptr)
  expect_named(saved, magenta:::cpp_contract$saved_state)
  expect_equal(saved$parameters_List$g_current_time, 21)
  expect_length(saved$population_List$Infection_States, 40)
  magenta:::sim_free(state$sim$Ptr)

  restored <- magenta:::sim_init_saved(saved)
  expect_type(restored$Ptr, "externalptr")
  again <- magenta:::sim_get_state(restored$Ptr)
  expect_equal(again$population_List$Infection_States, saved$population_List$Infection_States)
  expect_equal(again$parameters_List$g_current_time, 21)
  expect_invisible(magenta:::sim_free(restored$Ptr))
})

test_that("the bindings check their inputs", {
  expect_error(magenta:::sim_get_state("ptr"), "statePtr must be an external pointer")
  expect_error(magenta:::sim_free(NULL), "statePtr must be an external pointer")
  expect_error(magenta:::sim_init_saved("x"), "saved_state must be a list")
  expect_error(
    magenta:::sim_init_saved(list(population_List = 1)),
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
  magenta:::sim_free(state$sim$Ptr)
})

test_that("simulation_R dispatches init and saved-state lists, nested or not", {
  skip_on_cran()
  mpl <- model_param_list_create()
  eq_init <- equilibrium_init_create(
    age_vector = age_brackets(100, 40, TRUE), het_brackets = 5, ft = 0.4,
    EIR = 1, model_param_list = mpl
  )
  pl <- magenta:::param_list_simulation_init_create(
    N = 40, eqSS = equilibrium_ss_create(eqInit = eq_init),
    barcode_list = barcode_list_create(num_loci = 2, plaf = c(0.5, 0.5)),
    spatial_list = magenta:::spl_create(plaf = c(0.5, 0.5)),
    housekeeping_list = quiet_housekeeping_test(), drug_list = drug_list_create(),
    nmf_list = nmf_list_create(), vector_adaptation_list = vector_adaptation_list_create(2),
    mpl = mpl
  )
  # the old pipeline wrapped the init list in another list
  expect_warning(sim <- magenta:::simulation_R(list(pl), 16L), "deprecated")
  expect_type(sim$Ptr, "externalptr")
  expect_warning(saved <- magenta:::simulation_R(list(statePtr = sim$Ptr), 16L), "deprecated")
  expect_warning(
    magenta:::simulation_R(list(statePtr = sim$Ptr, finalizer = TRUE), 16L), "deprecated"
  )
  expect_warning(restored <- magenta:::simulation_R(list(savedState = saved), 16L), "deprecated")
  expect_type(restored$Ptr, "externalptr")
  magenta:::sim_free(restored$Ptr)
})
