## magenta 1.4.0

* Results change for every seeded stochastic run: the random number generator
  is seeded once per run instead of before every C++ call, and 1.4.0 carries
  the 1.3.8 model fixes.
* `pipeline()` returns the same list as 1.3.x, now classed `magenta_sim` with the run's
  `params`, `schedule` and `meta` as attributes (`sim_params()`, `sim_schedule()`,
  `sim_meta()`) and `print()`, `summary()` and `as.data.frame()` methods.
* `output =` takes a saver: `save_loggers()`, `save_summary()`, `save_humans()`,
  `save_full()` or your own function. The old output flags still work, with a warning.
* With `update_save = FALSE` a multi-year run steps year by year, so yearly inputs
  (`ft`, `plaf`, flags, importation) apply and temporal drug cycling happens.
* Yearly inputs shorter than the run repeat their last value; longer ones are an
  error. The schedule is built by the exported `build_schedule()`.
* Sequential drug cycling reads treatment failure from the loggers, so it works
  with every saver and without periodic saving.
* `sim_save()` and `sim_restore()` write and read a versioned saved state; old
  full saves still load.
* A resumed run uses the saved run's parameter lists and generator state; a
  differing explicit list is an error, and passing `seed` reseeds.
* Every parameter list constructor is exported, returns a classed list and is
  checked by a `validate_*()` naming the list and field; `cpp_contract` records
  the fields the C++ reads.
* `barcode_list_create()` defaults `plaf`, `prob_crossover` and `mutation_rate`
  to `num_loci` entries.
* Mosquito death and feeding rates under ITN/IRS come from the closed form
  `mu_fv_from_coverage()`, identical to the odin model output and far faster;
  `run_deterministic_model()` runs the model itself.
* `spl_grab()` keeps the year 2000 in a range that ends there.
* An all-zero `itn_cov` or `irs_cov` vector takes the no-intervention shortcut
  (`fv = 0.333`, `mu = 0.132`) as a scalar `0` always did, instead of running the model.
* `drug_list_create(number_of_resistance_loci = 1)` works.
* `model_param_list_create()` honours `irs_half_life` and `itn_half_life`;
  `pfpr_to_eir_heuristic(age_range = )` works; `microscopy_detected()` uses the
  age bracket bound; `spatial_type = "metapop"` stops with a clear message.
* A misspelt argument to `pipeline()` is reported as such;
  `equilibrium_init_create(quiet = TRUE)` is quiet.
* `simulation_R()` is deprecated; the `Simulation_*_cpp()` wrappers, the
  `param_list_simulation_*()` helpers and `%>%` are no longer exported.
* Removed `bootstrap()`, `kl()` and `generate_default_model()`. The other
  research-era helpers stay as internal functions pending a workflows phase;
  ggplot2, reshape2, cowplot and RColorBrewer are Suggests and redux is dropped.
* odin and dde come from CRAN; R >= 4.1 is required.
* Vignettes rebuilt against the new interface; every export has an example;
  fixed-seed characterisation tests, styler and lintr in CI; testthat edition 3.
* Two new vignettes: `Schedule` (schedules, parameter lists, the result
  object and reproducibility) and `Deterministic` (checks against the
  deterministic model).
* `Person::m_cA` is initialised, so set-up-only saves no longer carry an
  uninitialised value.

## magenta 1.3.8

* New vignettes `Loggers`, `DrugStrategies` and `NewDrugs`; `Resistance` and
  `Introduction` rewritten; sources in `vignettes_src/`.
* Fixed an initialisation crash in resistance simulations during the LPF calculation.
* The asymptomatic discount in the LPF calculation now applies; integer
  arithmetic had truncated it to zero.
* Gametocyte sterilisation applies to artemisinin sensitive strains, as
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
