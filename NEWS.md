# magenta 1.4.0

Restructuring of the R layer (see `R_CODE_REVIEW.md` in the repository for
the review and plan). The model is the same; the R side around it is
rewritten, and the changes that move results are listed first.

## Results change

* **One seed per run.** `pipeline()` seeds R's random number generator once,
  at the start, and every C++ call continues from wherever the generator
  is. Previously every call into the simulation (initialisation, each
  update, each state copy, each genetics sample) re-seeded with the same
  seed, so every update period replayed the same random stream. Every
  stochastic run with a given seed therefore gives different numbers from
  1.3.x: runs with periodic saves because the periods no longer repeat the
  stream, and single-period runs because the initialisation no longer
  re-seeds after the importation vectors have been drawn. The equilibrium
  initialisation is unchanged. Two different seeds now give different
  results, and the same seed the same, which the tests check.
* The generator state at the end of a run is kept in `meta$rng_state` and in
  every full state save. A run resumed with `saved_state_path =` restores it
  and so continues exactly where the saved run stopped; passing `seed`
  explicitly reseeds instead. Saves made before 1.4.0 have no generator
  state and are reseeded.
* **Yearly inputs apply year by year even without periodic saving.** With
  `update_save = FALSE` a multi-year run used to be one C++ call at the
  smallest non-zero `ft` for the whole run (and `Inf`, with a warning, when
  every `ft` was 0); `plaf`, `resistance_flag`, `mutation_flag` and the
  importation matrices beyond year one were never applied, and temporal
  drug cycling never happened. The run is now split into one period per
  year, each with that year's inputs; only the final slot is returned.
  `itn_cov` and `irs_cov` were already time-varying (they go through the
  deterministic mosquito model, which is solved over the whole run) and
  are unaffected. One-year runs are unchanged by this.
* **A resumed run continues with the saved run's parameter lists.** The
  core, drug, housekeeping, nmf, vector adaptation and barcode lists and
  the spatial settings (`spatial_type`, `island_imports_plaf_linked_flag`,
  `survival_percentage`, `oocyst_mean`, `oocyst_shape`) are taken from the
  save; the barcode and spatial settings used to be rebuilt from the
  resuming call. Supplying one of them that differs from the saved value is
  an error rather than a silent override. The yearly inputs are for the new
  years and may be given; `plaf`, `mutation_flag` and the importation
  vectors default to the saved run's last values.
* `spl_grab()` keeps the year 2000 when a year range ends there: 1995:2000
  gives six rows rather than five all-zero rows, and 2000:2000 works. A
  range ending in any other year is unchanged.

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

## Smaller changes

* The mosquito death and feeding rates under ITN/IRS are evaluated in closed
  form by `mu_fv_from_coverage()` instead of by running the deterministic
  model: they depend only on the coverage and the parameters, and the
  closed form is bit-identical to the model's output (tested against it).
  A two-year seasonal run's rates take under a millisecond rather than two
  seconds. The deterministic model is still in the package;
  `run_deterministic_model()` runs it and returns all of its outputs, for
  checking the simulation against it, and `mu_fv_create(full = TRUE)` does
  the same.

* An argument in `...` that `model_param_list_create()` does not take is
  reported as an unknown `pipeline()` argument, which is what a misspelt
  argument almost always is.
* `admin_match()` honours `quiet`, so `equilibrium_init_create(quiet = TRUE)`
  is quiet.
* Every exported function has a runnable example.

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
* Sequential drug cycling reads the treatment failure from the loggers every
  period returns, so it works with every saver and without periodic saving.
  It used to see it only through `save_summary()` with `update_save = TRUE`,
  never switching otherwise (and erroring on the missing value in some
  configurations). Runs that used `save_summary()` with periodic saving
  switch at the same times as before.
* `drug_list_create(number_of_resistance_loci = 1)` works; it previously
  failed building the resistance cost table because the one-column matrix
  of resistance barcodes was dropped to a vector. The table for two or more
  loci is unchanged.

## Removed

* `generate_default_model()`, which nothing called (the deterministic model
  is driven from `mu_fv_create()`).

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
