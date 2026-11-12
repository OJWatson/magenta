# magenta 1.4.0

Restructuring of the R layer (see `R_CODE_REVIEW.md` in the repository for
the review and plan). Simulation results are unchanged for every
configuration covered by the characterisation tests; the few inputs whose
handling has changed are listed under "Schedule" below.

## Parameter lists

* Every parameter list constructor (`model_param_list_create()`,
  `drug_list_create()`, `drug_create()` and the preset drugs,
  `housekeeping_list_create()`, `nmf_list_create()`,
  `vector_adaptation_list_create()` and `barcode_list_create()`) is now
  exported, returns a classed list and is checked by a matching
  `validate_*()` function with messages that name the list and field, e.g.
  `drug_list$partner_drug_ratios must have length 2, not 3`.
  The lists are still plain lists underneath.
* The field names the C++ simulation reads are recorded once, in
  `cpp_contract`, and a test checks that each constructor produces exactly
  that set and that the C++ sources read no others.
* `model_param_list_create()` is written as its arguments plus the derived
  quantities; the output is the same list as before, in the same order.
* `barcode_list_create()` defaults `plaf`, `prob_crossover` and
  `mutation_rate` to `num_loci` entries rather than 24.

## Schedule

* `build_schedule()` turns the yearly inputs to `pipeline()` (`ft`,
  `itn_cov`, `irs_cov`, `resistance_flag`, `mutation_flag`, `plaf` and the
  two importation matrices) into one data frame with a row per update
  period. `pipeline()` runs from it and records the drug used in each
  period. The schedule is returned in the `meta` attribute.
* One padding rule replaces the five that were in use: a yearly input shorter
  than the simulation repeats its *last* value to the end, and one longer
  than `ceiling(years)` is an error. This changes results only for inputs
  that relied on the old behaviour: importation matrices with fewer rows than
  years were padded at the *front* with their first row; `plaf` matrices were
  padded with their *second* row; `ft`, `itn_cov` and `irs_cov` shorter than
  `years` were recycled from the start in the deterministic mosquito model
  and gave `NA` in the stochastic one. Inputs with one value or one value per
  year, which is every configuration in the tests and vignettes, are
  unaffected.
* Each period's inputs, including the last period's, are now applied. Before,
  a year change that fell on the final period (for example the second year of
  a two year run with `update_length = 365`) was silently skipped and the
  previous year's `ft`, flags and importation rates were used.
* `update_save = TRUE` with a single period (for example `years = 1,
  update_length = 365`) now runs instead of failing with
  `object 'i' not found`.
* A vector of calendar years (`years = 2000:2015`) is accepted again, with or
  without `use_historic_interventions`.
* The constant mosquito death and biting rates used when there are no vector
  interventions apply whenever every `itn_cov` and `irs_cov` value is zero,
  rather than only when each is the single value `0`.

## Output, results and saved states

* `pipeline()` takes one `output` argument in place of the seven output
  booleans. It is a *saver*: `save_loggers()` (the default),
  `save_summary()`, `save_humans()` or `save_full()`, or a
  `function(state, period, params)` of your own. Each saver has a `final`
  argument saying what the last element holds (`"loggers"`, `"full"` or
  `"humans"`). The old flags (`full_save`, `full_update_save`,
  `human_only_full_save`, `human_update_save`, `summary_saves_only`,
  `genetics_df_without_summarising`, `only_allele_freqs`, `save_lineages`,
  `mean_only`, `sample_size`, `sample_states`, `sample_reps`, `age_breaks`
  and `update_save_func`) still work for this release: they are mapped onto
  the equivalent saver with a warning that spells out the replacement.
  `update_saves()` is gone.
* `pipeline()` returns a `magenta_sim` with elements `output` (what the saver
  kept; with `update_save = FALSE` the final slot alone), `params`,
  `schedule` and `meta`, with `print()`, `summary()` and `as.data.frame()`
  methods. `output` holds exactly what `pipeline()` used to return, so
  `out[[i]]` becomes `out$output[[i]]`, and the `meta` attribute becomes
  `out$meta`. `set_up_only = TRUE` still returns the initialised state.
* `sim_save(sim, path)` writes the full state at the end of a run as one
  versioned object together with the run's parameters, schedule and
  metadata; `sim_restore(path)` reads and checks it. `pipeline()` accepts
  either the path or the restored object as `saved_state_path`, and still
  accepts the plain full saves written by earlier versions.
* A human-only final save from a single-period run (`update_save = FALSE`)
  now includes `ID`, as it always did from an update run. Per-period full
  saves from `save_full()` carry `saved_inputs`, so any of them can be
  restarted from.
* The vignettes are written against the new interface.

## Simulation bindings

* `simulation_R()` is deprecated and no longer exported. `pipeline()` calls
  the C++ entry points through the internal bindings `sim_init()`,
  `sim_update()`, `sim_get_state()`, `sim_free()` and `sim_init_saved()`.
  The raw `Simulation_*_cpp()` wrappers, the `param_list_simulation_*()`
  helpers and the `%>%` re-export are no longer exported.
* `pipeline()` is now a short wrapper over `sim_setup()`, `run_schedule()`,
  `prepare_state_save()` and `collect_output()`.
* `Person::m_cA` is initialised, so the asymptomatic infectivity reported for
  a never-updated individual in a saved state is `0` rather than whatever was
  in memory.

## Bug fixes

* `model_param_list_create()` now honours its `irs_half_life` and
  `itn_half_life` arguments; previously both were silently overwritten with
  the defaults.
* `pfpr_to_eir_heuristic(age_range = )` works; it previously referenced the
  equilibrium before it was computed.
* `microscopy_detected()` uses the lower bound of `age_brackets` rather than
  calling the `age_brackets()` function.
* Asking for `spatial_type = "metapop"` now stops with a clear message
  instead of failing inside the C++ initialisation. Metapopulation
  simulations were never implemented.

## Removed

* The research-era analysis code (`bootstrap()`, `kl()`, `pool_res()`,
  `create_res_list()`, the `plot_*()` functions, `heatmap_ordered_binary_plot()`,
  `sample_coi()`, `coi_age_plot_sample_x()`, `summary_data_frames_from_sims()`,
  `country_seasonal_profiles()` and unused IBD helpers) has been removed. It
  was never called by the package, tests or vignettes. `convert_barcode_vectors()`
  and `microscopy_detected()` are kept.
* Dropped dependencies: ggplot2, reshape2, redux, cowplot and RColorBrewer.
  odin and dde are taken from CRAN. The package now requires R >= 4.1.

## Development

* Fixed-seed characterisation tests pin `pipeline()` output for several
  small configurations so later refactors can be checked as
  behaviour-preserving. The package uses the third testthat edition.
* Sources are formatted with styler and held to the linters in `.lintr`,
  which run in CI.

## magenta 1.3.8

* New vignettes: `Loggers`, `DrugStrategies` and `NewDrugs`. `Resistance` and
  `Introduction` rewritten, and vignette sources moved to `vignettes_src/`.
* Bug fix to initialisation that caused resistance simulations to crash during
  the LPF calculation.
* Bug fix to the asymptomatic discount in the LPF calculation, which was
  computed with integer arithmetic and so never applied.
* Gametocyte sterilisation now applies to artemisinin sensitive strains, as
  documented, rather than resistant ones.

## magenta 1.3.7

* Fix to saved state initialisation of barcode bitset sizing on Windows

## magenta 1.3.6

* Reloading models saved to file working again
* New function `update_saved_state_barcode_plaf` to update plaf in saved simulation

## magenta 1.3.5

* Seasonality patch so that admin and country name specification to `pipeline`
correctly creates seasonal simulations.

## magenta 1.3.4

* Non-malarial fevers were not going through correct drug allocation so could not
track who was getting what drug when it was an NMF.

## magenta 1.3.3

* Prophylaxis fixes. Individuals who failed treatment and recrudesce were not getting
full duration of prophylaxis previously (short by a couple days) 

## magenta 1.3.2

* NMF fixes. Individuals not infected who are treated due to NMF go straight to 
an elongated prophylaxis. 

## magenta 1.3.1

* Importations can have same barcode as PLAF via `island_imports_plaf_linked_flag`

## magenta 1.3.0

* Major changes inplace for simulating resistance

## magenta 1.2.0

* Documentation for base model
* `dde::dopri` fix that caused deterministic model to fail below very small EIR

## magenta 1.1.0

* New `param_list_simulation_finalizer_create` used within `pipeline` to return Rcpp memory

## magenta 1.0.0

* New major version for stable base model for publication. Versioning now incremental. 
* New version fixes minor changes to historic strain clearing in resistance

## magenta < 1.0.0

* Iterative testing with versioning only used for distinguishing testing on cluster
