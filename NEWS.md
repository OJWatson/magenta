# magenta 1.3.8

## Bug fixes

* Fixed a crash (`SIGFPE`, integer division by zero) in
  `Person::late_paristological_failure_boolean()` that terminated the R session
  whenever a resistance simulation was run with treatment, i.e. any
  `drug_list_create(resistance_flag = TRUE)` with `ft > 0`. A strain whose
  status change day equalled its acquisition day divided by zero. The window is
  now guarded; behaviour is unchanged wherever it previously did not crash.
* `drug_list_update()` now says which drug sequential cycling has switched to,
  rather than printing a bare number.

## Documentation

* New vignette **Model output, loggers and saved state**, documenting every
  logger and the save and restart workflow. This replaces the "Loggers"
  vignette that the README and `Introduction` had linked to but which had never
  existed.
* New vignette **Drug deployment strategies**, covering multiple first line
  therapies, sequential and temporal cycling, and non-malarial fevers.
* New vignette **Adding a new drug: single low-dose primaquine**, covering how
  to define a therapy the package does not ship with, and what would be needed
  to represent a transmission-blocking drug.
* **Simulating antimalarial resistance** substantially expanded: resistance
  loci and their barcode positions, the anatomy of a drug, prophylaxis curves,
  standing variation, importation, mutation and fitness costs.
* `Introduction` now gives the correct installation command. It previously told
  users to install `rdhs`, as did the README and the package website.
* Broken vignette cross-references fixed; figures are referenced relatively
  rather than through raw GitHub URLs.
* Vignette build moved to `vignettes_src/build_vignettes.R`, which no longer silently
  deletes the whole vignette when the simulation progress bar is absent.

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
