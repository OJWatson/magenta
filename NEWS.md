## magenta 1.5.0

* Severe malaria: acquired and maternal severe immunity (`IVA`, `IVM`) as in
  Griffin et al. 2015, `Severe` and `Severe_05` loggers, and `severe_rates()`
  and `severe_outcomes()`.
* Infection, immunity boosting, maternal immunity and mosquito biting brought
  in line with malariasimulation (see the Equivalence vignette). This changes
  model outputs.
* Run from Imperial's malariasimulation site files with `site_to_pipeline()`,
  `site_ft()` and `site_seasonality()` (see the "Running from a site file"
  vignette). Only treatment, bed nets and IRS are used.
* `pipeline()` and `equilibrium_init_create()` take `seasonality` directly.
* New `itn_decay_clock` in `model_param_list_create()`: "breakpoint" (the
  default, as before), "campaign" or "constant".
* `drug_create()` and the `drug_create_*()` functions are exported.

Behaviour changes and bug fixes:

* `itn_half_life` and `irs_half_life` now take effect (defaults unchanged).
* Non-malarial fevers in uninfected people only give drugs to those who seek
  treatment (`ft`).
* `ft = 0` no longer runs as `ft = Inf`.
* `Q0`, `mu0` and the infection and immunity parameters are now read by the
  individual model on a fresh start, not only from saved states.
* A `seasonality` argument now reaches the individual model.
* Initially infected people can realise new infections, and `ICM_init` is
  kept in saved states.

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
