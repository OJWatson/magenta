# magenta 1.4.0

Tidy-up of the R layer ahead of a larger refactor (see `R_CODE_REVIEW.md` in
the repository for the review and plan). No change to simulation results.

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
