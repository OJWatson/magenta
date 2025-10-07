run_saved_state_simulation <- function(saved_state_path = NULL, set_up_only = FALSE) {
  housekeeping <- housekeeping_list_create(
    quiet = TRUE,
    quiet_test = TRUE,
    cluster = TRUE,
    clear_up = TRUE
  )
  
  base_args <- list(
    EIR = 0.1,
    years = 1,
    N = 200,
    itn_cov = rep(0, 1),
    ft = rep(0.4, 1),
    num_loci = 2,
    update_save = TRUE,
    update_length = 30,
    human_update_save = TRUE,
    summary_saves_only = TRUE,
    genetics_df_without_summarising = TRUE,
    save_lineages = TRUE,
    full_save = TRUE,
    seed = 42L,
    housekeeping_list = housekeeping,
    saved_state_path = saved_state_path,
    set_up_only = set_up_only
  )
  
  do.call(pipeline, base_args)
}



test_that("pipeline full saves contain restart metadata", {
  skip_on_cran()
  
  full_run <- run_saved_state_simulation()
  expect_gt(length(full_run), 0)
  
  final_state <- full_run[[length(full_run)]]
  expect_true(is.list(final_state$saved_inputs))
  expect_true(all(c(
    "barcode_list",
    "spatial_list",
    "housekeeping_list",
    "drug_list",
    "vector_adaptation_list",
    "nmf_list",
    "core_parameter_list"
  ) %in% names(final_state$saved_inputs)))
})


test_that("saved state can be reloaded for continuation or inspection", {
  skip_on_cran()
  
  full_run <- run_saved_state_simulation()
  final_state <- full_run[[length(full_run)]]
  
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(final_state, tf)
  
  resumed <- run_saved_state_simulation(saved_state_path = tf)
  meta <- attr(resumed, "meta")
  expect_true(is.list(meta))
  expect_equal(meta$seed, 42L)
  expect_equal(meta$call$saved_state_path, tf)
  expect_gt(length(resumed), 0)
  
  setup_only <- run_saved_state_simulation(saved_state_path = tf, set_up_only = TRUE)
  expect_equal(setup_only$population_List, final_state$population_List)
  expect_equal(
    setup_only$populations_event_and_strains_List,
    final_state$populations_event_and_strains_List
  )
  expect_equal(setup_only$saved_inputs$barcode_list, final_state$saved_inputs$barcode_list)
})