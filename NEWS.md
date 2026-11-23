## magenta 1.3.9

* Severe malaria. Each person now carries acquired (`IVA`) and maternal
  (`IVM`) severe disease immunity, following Griffin et al. 2015 Proc R Soc B
  282:20142657, with the parameter values used in malariasimulation. They are
  boosted and decay in the same way as `ICA` and `ICM`, and are initialised
  from equilibrium.
* New loggers `Severe` and `Severe_05`: the expected number of severe cases
  per day (the sum of the probability that each new infection is severe),
  averaged over the logging period like `Incidence`. No new random draws are
  made, so all existing outputs are unchanged. `IVA` and `IVM` are logged per
  person and are kept in saved states. Older saved states and equilibrium
  states still load, starting with no severe immunity.
* New functions `severe_rates()`, to give clinical and severe disease per
  person per year, and `severe_outcomes()`, to give hospitalisations and
  deaths using the same assumptions as postie.
* Maternal clinical immunity at birth (`ICM_init`) is now saved and restored with the model state; it was previously left uninitialised after a reload.

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
