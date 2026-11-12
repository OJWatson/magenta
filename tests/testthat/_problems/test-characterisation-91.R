# Extracted from test-characterisation.R:91

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
nl <- 4
plaf <- matrix(c(0.5, 0.5, 0, 0, rep(0, nl)), ncol = nl, byrow = TRUE)
drug_list <- drug_list_create(
  resistance_flag = c(FALSE, TRUE),
  number_of_resistance_loci = nl,
  artemisinin_loci = 3,
  cost_of_resistance = rep(0.99, nl),
  drugs = list(drug_create(
    prob_of_lpf = seq(1, 0.5, length.out = 2^nl),
    barcode_res_pos = seq_len(nl) - 1,
    prophylactic_pos = nl - 1
  )),
  partner_drug_ratios = 1
)
out <- suppressMessages(pipeline(
  EIR = 2, years = 2, N = 200, seed = 4L, num_loci = nl, ft = 0.5,
  plaf = plaf, mutation_flag = c(FALSE, TRUE), mutation_rate = rep(1e-3, nl),
  mutation_treated_modifier = 10,
  spatial_type = "island", island_imports_plaf_linked_flag = TRUE,
  spatial_incidence_matrix = c(1, 0), spatial_mosquitoFOI_matrix = c(1, 0),
  update_save = TRUE, update_length = 60,
  human_update_save = TRUE, summary_saves_only = TRUE,
  genetics_df_without_summarising = TRUE, save_lineages = TRUE,
  sample_size = c(50, 100), sample_states = c(1, 2, 4),
  drug_list = drug_list,
  nmf_list = nmf_list_create(nmf_flag = TRUE, prob_of_testing_nmf = 0.5),
  housekeeping_list = quiet_housekeeping()
))
expect_true(all(c("af", "lineage", "mutations") %in% names(out[[1]])))
expect_length(out[[1]]$lineage, 2^nl)
expect_snapshot_value(snapshot_summary(out), style = "json2", tolerance = 1e-6)
