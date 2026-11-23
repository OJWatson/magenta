## magenta 1.4.0

* Severe malaria: acquired and maternal severe immunity (`IVA`, `IVM`) as in
  Griffin et al. 2015, new `Severe` and `Severe_05` loggers, and new functions
  `severe_rates()` and `severe_outcomes()`.
* Infection, immunity boosting and mosquito biting brought in line with
  malariasimulation, which changes model outputs:
  * Bug fix to initialisation, which stopped initially infected people from
    realising new infections and lowered incidence for the first 2-3 years.
  * Bites no longer fail to infect because `IB` was boosted recently.
  * Infection, clinical and severe disease use `IB`, `ICA` and `IVA` before
    the day's boost, plus 0.5 if non-zero.
  * `ICA`, `ID` and `IVA` are only boosted, and severe cases only counted, by
    infections that are not cleared by prophylaxis.
  * Mosquito feeding intervals are geometric with mean 1/fv rather than a fixed
    3 days, and the EIP is 10 or 11 days.
  * Infectivity to mosquitoes lags the infection state by 12.5 days
    (`delay_gam`), as in the Griffin et al. model, so people stay infectious
    for 12.5 days after recovery or treatment. Saved states now keep each
    person's recent infection states and cleared strains.
* Bug fix: maternal immunity at birth (`ICM_init`) is now kept in saved states.

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
