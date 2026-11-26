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
  * Reinfection in A no longer lengthens the time spent in A. The new
    infection sets when the person moves to U, rather than the later of the
    new and current times, as A -> U is memoryless in the Griffin et al. model.
    Strains in A are kept from outlasting the person's move to U.
  * Whether an infectious bite infects is decided when the infection starts,
    after the latent period, as in the Griffin et al. model, rather than on
    the day of the bite. Pending infections are kept through recovery and
    treatment, and `ICA`, `ID` and `IVA` boosts, severe cases and clinical
    incidence are counted on the day the infection starts. This raises
    transmission by about 3-6%.
  * `pipeline()` starts from an equilibrium with fine age groups (0.1 years
    to age 10, then 0.5 years), which raises the mosquito density by about
    3% to match the deterministic model.
  * Newborns take maternal immunity from a random 20-year-old in the same
    biting heterogeneity group, as in malariasimulation and the Griffin et al.
    model, rather than the mean immunity of 15-25 year olds scaled by their
    own biting rate. This raises infant maternal immunity by 4-6% and lowers
    severe incidence by about 2% at EIR 32 and 5% at EIR 64.
* Bug fix: maternal immunity at birth (`ICM_init`) is now kept in saved states.
* Bug fix: a new treated strain no longer sets the day of the next strain
  state change, as treated strains have no scheduled change.
* Bug fix: non-malarial fevers in uninfected people only give drugs (and so
  prophylaxis) to people who seek treatment (`ft`) and are not tested, rather
  than to everyone.
* Bug fix: `pipeline()` with `ft = 0` ran the simulation with `ft = Inf`, so
  every clinical case and non-malarial fever was treated.

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
