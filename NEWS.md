# magenta 1.3.8

## Bug fixes

* Fixed a crash (`SIGFPE`, integer division by zero) in
  `Person::late_paristological_failure_boolean()` that terminated the R session
  whenever a resistance simulation was run with treatment, i.e. any
  `drug_list_create(resistance_flag = TRUE)` with `ft > 0`.

  The cause was in initialisation, not in the simulation loop. Individuals
  starting the simulation in the `TREATED` state had their strains given a
  status change day hardcoded to `0` and an acquisition day drawn as
  `0 - runiform_int_1(0, dur_T)`. That draw is inclusive of both bounds, so
  roughly one strain in six got an acquisition day of `0` as well, leaving a
  zero length course that divided by zero as soon as that individual's
  treatment outcome was evaluated. The draw now starts at `1`, so the
  acquisition day is always strictly before the status change day, which is
  what every in-run state change already guaranteed by drawing
  `current_time + 1` or later. The division is also guarded, for saved states
  written by earlier versions.
* The discounting of an asymptomatic strain's probability of surviving
  treatment, by how far through its course that strain is, now takes effect. It
  was computed with integer arithmetic and truncated to zero for any strain
  partway through its course, so it had silently done nothing. Which strains
  survive a late parasitological failure is unchanged.
* Gametocyte sterilisation is now applied to artemisinin *sensitive* strains
  rather than resistant ones, matching the documented intent that artemisinin
  sterilises the gametocytes of parasites sensitive to it. Switching
  `gametocyte_sterilisation_flag` on therefore now gives the artemisinin
  resistant parasite a transmission advantage, where previously it penalised
  it. Only affects simulations with that flag set, which is not the default.
* `drug_list_update()` now says which drug sequential cycling has switched to,
  rather than printing a bare number.

## Known issues

* `oocyst_reduction_by_artemisinin` is written back into the shared pre-drawn
  oocyst frequency vector in `src/mosquito.cpp` rather than applied to a copy,
  so a reduced value persists for later mosquitoes reusing that element until
  the vector is regenerated the following year. Only affects simulations with
  `gametocyte_sterilisation_flag` set.
* Gametocyte sterilisation cannot reduce a host's infectiousness. The
  multiplier is applied to `m_active_strain_contribution`, which `sample1()`
  treats as relative weights among the strains one host carries, and a host
  carrying a single strain returns before those weights are read at all. A host
  whose strains are all sensitive therefore has every weight scaled equally and
  nothing changes. Sterilisation only reweights which strain a mixed host
  transmits. There is consequently no population level signal from switching it
  on, in either direction, at the scale the test suite can run.
* In `Person::sample_two_barcodes()`, `m_female_contribution_sum` is computed
  by summing `m_active_strain_contribution` rather than
  `m_active_female_strain_contribution`, and is then used to sample from the
  female vector. The female vector also never receives the resistance fitness
  modifier or the sterilisation multiplier, so the two sums diverge whenever
  either applies and the female gamete draw is biased.
* The gate deciding which additional strains survive a late parasitological
  failure reads the loop variable after its loop has ended, so one strain's
  value is applied to all of them. Instrumented runs show it open on 98% of
  late parasitological failures, so it is close to a no-op, but it is not the
  per-strain test it appears to be.

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
