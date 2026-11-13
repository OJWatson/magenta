# magenta: R-side code review and tidy-up plan

Date: 2026-11-10. Reviewed at commit `dba6c5c` (v1.3.7), `R/` only (the C++ is out
of scope except where it defines the contract the R code must honour).

How this was produced: every file in `R/` and `tests/` was read in full; function
lengths were measured from parsed source; `lintr` was run over `R/`; `covr` was run
twice (with and without `NOT_CRAN`, because two tests are guarded by
`skip_on_cran()`); `R CMD check --no-tests` was run; and the git history was used to
date each file. Numbers below come from those runs, not estimates.

---

## 1. Executive summary

The package has a sound shape. The idea that grew over the PhD, of building named
parameter lists with constructor functions and validating them before they cross into
C++, is the right one, and the newest code (`update_saved_state_barcode_plaf`, 2025) shows
what the whole package should look like: small functions, explicit validation, tests
for the error paths. The problem is that the shape is only half-built, and the
pipeline grew around it rather than on top of it.

The ten things that matter most, in order:

1. **`pipeline()` is a 685-line function with 48 arguments plus `...`**
   ([pipeline.R:131](R/pipeline.R:131)). It contains a ~75-line block duplicated
   verbatim between the "fresh start" and "resume from saved state" branches
   ([pipeline.R:221-297](R/pipeline.R:221) vs [pipeline.R:325-401](R/pipeline.R:325)),
   a "save the full state" block written three times with small inconsistencies, and
   ~150 lines of commented-out metapopulation/redis code. It should become a thin
   wrapper over five or six named steps (section 6).

2. **The random seed is reset before every simulation chunk.** `simulation_R` calls
   `set.seed(seed)` before each Rcpp call ([mainFunctions.R:288](R/mainFunctions.R:288),
   [317](R/mainFunctions.R:317), [345](R/mainFunctions.R:345)), and the C++ draws from
   R's RNG (`R::unif_rand`, `src/probability.cpp`). With `update_save = TRUE` every
   30-day (or whatever `update_length`) chunk therefore replays the identical random
   stream, and `pop_strains_df` re-seeds again before each population sample. This is
   not a style issue; it is worth checking what it does to within-run variance before
   anything else is touched, and it needs its own NEWS entry when fixed because it
   changes results.

3. **Output shape is controlled by seven booleans** (`full_save`,
   `full_update_save`, `human_only_full_save`, `update_save`, `human_update_save`,
   `summary_saves_only`, `genetics_df_without_summarising`, plus `set_up_only` and
   `only_allele_freqs`). The returned object is an unclassed list whose structure
   depends on the combination, with metadata hidden in an attribute. One `output =`
   argument choosing a named "saver", and a classed result with `print()`, would
   remove most of the branching in both `pipeline()` and `update_saves()`.

4. **The parameter-list layer is inconsistent about what is public.**
   `drug_list_create` is exported but `drug_create_*`, `housekeeping_list_create`,
   `nmf_list_create` and `vector_adaptation_list_create` are not, even though they are
   `pipeline()` defaults; the vignette itself has to write `magenta:::drug_create_dhappq()`.
   Meanwhile the raw Rcpp wrappers (`Simulation_Init_cpp` etc.) and `%>%` are exported.

5. **Per-year schedules are handled ad hoc.** `ft`, `itn_cov`, `irs_cov`,
   `resistance_flag`, `mutation_flag`, the `plaf` matrix and the two spatial matrices
   are each stretched to `ceiling(years)` by a different piece of code with different
   padding rules (`spl_matrix_check` prepends row 1; `plaf_matrix_check` appends row 2;
   `mu_fv_create` has its own `length_check`). One `schedule` data frame built and
   validated once would replace all of it.

6. **R-side test coverage is 47.6%**, and what is covered is covered by four
   end-to-end simulations. `COI_df_create` (the genetics summariser the vignettes rely
   on, 248 lines) has 0% coverage. The default call path in the Introduction vignette
   (`update_save = FALSE`) is untested. There are no unit tests of any pure function.
   Three tests use a random seed, so they can fail by chance.

7. **About 1,600 lines are research-era code that no longer belongs in the package**:
   `bootstrap.R`, `kl.R`, `cluster_utils.R` (hard-coded `M:/OJ/...` paths, didehpc task
   objects), and most of `plotting.R`. None is referenced from `R/`, `tests/` or the
   vignettes, all is at 0% coverage, and some of it has latent bugs
   (`microscopy_detected` calls `age_brackets(1)` where it meant `age_brackets[1]`).

8. **`model_param_list_create` silently ignores two of its arguments**
   (`irs_half_life`, `itn_half_life` are overwritten at
   [model_param_list_create.R:281-282](R/model_param_list_create.R:281)), and its 72
   arguments are copied into a list one line at a time for 115 lines.

9. **Style is inconsistent in every dimension lintr measures**: 3,949 lints
   (with the naming, indentation and commented-code linters already switched off),
   four naming conventions in use (`sim.out`, `mp.list`, `statePtr`, `drug_list`),
   `T`/`F` for `TRUE`/`FALSE`, `&` for `&&`, `class(x) == "list"`, `1:length(x)`, and
   mixed `=`/`<-` assignment. This is the cheapest thing to fix (one `styler` commit
   plus a `.lintr` file) and the one that most changes how the package reads.

10. **Package metadata is a decade old**: `Depends: R (>= 2.10)`,
    `SystemRequirements: C++11` (now a check NOTE), `redux` imported but only used in
    commented-out code, `dde` imported but never referenced, `Remotes:` for `odin`
    which is on CRAN, README says to install `rdhs`, a Travis badge, no `@examples`
    anywhere (`R CMD check` reports "checking examples ... NONE"), `codemeta.json` at
    1.3.0, and vignette figure files with spaces in their names (a portability
    WARNING).

None of this is unusual for a model that has carried a research programme for ten
years. The good news is that the C++ boundary is narrow and well defined (five entry
points, each taking one named list), which makes the R side cheap to restructure
without touching the model.

---

## 2. The package by numbers

| Measure | Value |
|---|---|
| R source files | 20 (`RcppExports.R` and `odin.R` generated) |
| Lines in `R/` | 6,109 |
| Functions defined | 107 (90 hand-written) |
| Functions over 100 lines | 11 |
| Longest function | `pipeline()`, 685 lines, 48 formal args + `...` |
| Most arguments | `model_param_list_create()`, 72 |
| Exported functions | 23 (of which 7 are raw Rcpp wrappers, 1 is `%>%`) |
| Test files / `test_that()` blocks | 5 / 14 |
| Test suite wall-clock (`NOT_CRAN=true`) | 19.6 s, all passing |
| Line coverage (tests, `NOT_CRAN=true`) | total 58.3%, R 47.6%, C++ 67.6% |
| Line coverage without `NOT_CRAN` | total 57.7%, R 43.7% (saved-state tests skipped) |
| Functions at 0% coverage | 41 of 90 |
| lintr findings (`R/`, default linters minus naming/indent/commented-code) | 3,949 |
| `R CMD check --no-tests` | 1 WARNING (non-portable vignette file names), 3 NOTEs |
| Commits | 181, 2016 to 2025; 60 of them in 2017-2019 |
| `TODO`/`FIXME` markers in `R/` | 2 |

Function lengths, top of the list:

| Function | File | Lines | Args |
|---|---|---|---|
| `pipeline` | pipeline.R | 685 | 48 + `...` |
| `equilibrium_init_create` | equilibrium_init_create.R | 312 | 8 |
| `COI_df_create` | universe_utils.R | 248 | 9 |
| `summary_data_frames_from_sims` | plotting.R | 211 | 6 |
| `pool_res` | cluster_utils.R | 207 | 9 |
| `model_param_list_create` | model_param_list_create.R | 205 | 72 |
| `simulation_R` | mainFunctions.R | 164 | 2 |
| `plot_resistance_lineages` | plotting.R | 130 | 3 |
| `update_saves` | pipeline_utils.R | 128 | 19 |
| `generate_default_model` | equilibrium_state_create.R | 112 | 5 |
| `pfpr_to_eir_heuristic` | pfpr_to_eir.R | 108 | 7 |

Lint findings by type (top 10 of 20):

| Linter | Count | What it means here |
|---|---|---|
| trailing_whitespace | 1,174 | mechanical; `styler` fixes it |
| infix_spaces | 736 | `a=b`, `x+y`, `nums%>%unlist` |
| commas | 735 | `f(a,b,c)` |
| object_usage | 328 | mostly "no visible binding"; includes two real bugs (section 4.6) |
| spaces_left_parentheses | 241 | `if(` and `for(` |
| brace | 213 | `else` on its own line, `{` on its own line |
| paren_body | 144 | `function(x){` |
| return | 94 | explicit `return()` at end of function (fine either way, but pick one) |
| line_length (>100) | 84 | some lines are 180+ characters |
| pipe_consistency | 74 | `%>%` with and without `()` after the function |

---

## 3. How the R layer is wired today

Reading the code, the data flow is:

```
pipeline(...)                                   [pipeline.R]
  |
  |-- model_param_list_create(...)               72 scalar epi params -> mpl
  |-- equilibrium_init_create(age, het, ft, EIR, country, admin, mpl)   -> eqInit (deterministic equilibrium)
  |-- equilibrium_ss_create(eqInit)              -> eqSS  (the 19 matrices C++ wants)
  |-- barcode_list_create / spl_create / drug_list_create / nmf_list_create /
  |   vector_adaptation_list_create / housekeeping_list_create          -> sub-lists
  |-- param_list_simulation_init_create(N, eqSS, <sub-lists>)           -> pl (9 named elements)
  |-- simulation_R(pl, seed)  ->  Simulation_Init_cpp(pl)                -> list(Ptr = <externalptr>, ...)
  |
  |-- mu_fv_create(eqInit, ft, itn, irs, years)  runs the odin ODE -> daily mu, fv vectors
  |
  |-- loop over update periods:
  |     param_list_simulation_update_create(years, ft, mu, fv, statePtr, spatial, drug, barcode)
  |     simulation_R(pl2) -> Simulation_Update_cpp -> list(Ptr, Loggers)
  |     update_saves(res, i, sim.out, <11 flags>)   -> res[[i]]  (loggers | summary | humans | full)
  |        |-- pop_strains_df(Ptr) -> population_get_genetics_df_n (C++)
  |        |-- COI_df_create(df) | pop_alf | lineages
  |     drug_list_update(drug_list, year, treatment_failure)
  |
  |-- optional: param_list_simulation_get_create -> Simulation_Get_cpp -> full state dump
  |-- optional: param_list_simulation_finalizer_create -> Simulation_Finalizer_cpp (free memory)
  |
  '-- attr(res, "meta") <- list(call, seed, Ptr, version, times)
```

The C++ boundary is clean: five entry points, each taking one named list, with the
field names read by `Rcpp::as<...>(list["name"])` in `src/main_init.cpp`,
`src/main_update.cpp`, `src/main_saved_init.cpp`. That list of names *is* the
contract, and it is currently enforced only by `identical(names(param_list), c(...))`
checks inside `simulation_R` and by whatever `Rcpp::as` does when a name is missing
(an uninformative "index out of bounds" error).

The saved-state path is the same diagram with `Simulation_Saved_Init_cpp` in place of
`Simulation_Init_cpp`, except that `pipeline()` still recomputes `equilibrium_init_create`
and `equilibrium_ss_create` on resume solely to get an `eqInit` for `mu_fv_create`.

---

## 4. Findings in detail

### 4.1 `pipeline()`

**Duplication.** Lines [221-297](R/pipeline.R:221) (fresh start) and
[325-401](R/pipeline.R:325) (resume) are the same ~75 lines: age brackets, `ft_odin`,
`equilibrium_init_create`, `equilibrium_ss_create`, `plaf_matrix_check`,
`barcode_list_create`, mutation flag expansion, spatial type coercion,
`spl_matrix_check` x2, `spl_create`, resistance flag expansion. The only difference is
that the resume branch *overwrites* `barcode_list`, `spatial_list` and friends that it
just read from `saved_inputs`, so the saved inputs are partially discarded (the
barcode and spatial lists are rebuilt from the *current* call's arguments; the drug,
nmf, vector adaptation and housekeeping lists come from the file). That asymmetry is
not documented and is almost certainly not intended.

The "grab the full state and attach `saved_inputs`" block appears three times
([422](R/pipeline.R:422), [717](R/pipeline.R:717), [765](R/pipeline.R:765)).
In two of them `human_vars` includes `"ID"`; in the third
([777](R/pipeline.R:777)) it is redefined without it, so a `human_only_full_save` from a
non-`update_save` run has a different shape from one from an `update_save` run.

**Dead code.** [456-533](R/pipeline.R:456) and [640-676](R/pipeline.R:640) are
commented-out redis/metapopulation code; [264-270](R/pipeline.R:264) and
[369-375](R/pipeline.R:369) likewise. The
`spatial_type == "metapop"` branch never converts the string to `2`
([264](R/pipeline.R:264)), so a user passing `"metapop"` gets `spl_create("metapop")`,
which falls through its `if/else` and returns `NULL` ([equilibrium_state_create.R:153-155](R/equilibrium_state_create.R:153)),
and the eventual failure is an Rcpp subscript error rather than the intended
`stop("metapopulation simulation not implemented yet")` at [458](R/pipeline.R:458).

**Latent errors.**

- If `round(years * 365 / update_length) == 1` (e.g. `years = 1, update_length = 365`
  with `update_save = TRUE`), the `for` loop at [576](R/pipeline.R:576) is skipped and
  `i` is never assigned, so [695](R/pipeline.R:695) fails with "object 'i' not found".
- `use_historic_interventions = TRUE` is documented as "provide the years as a year
  range, e.g. 2000:2015", but everything after [202](R/pipeline.R:202) treats `years` as
  a scalar count (`ceiling(years)`, `round((years * 365) / update_length)`, `length(res) <- ...`).
  A vector `years` therefore errors. The feature as documented cannot work; the test
  that sets the flag passes a scalar, so the branch is never entered (0% coverage on
  [204-210](R/pipeline.R:204)).
- In the non-`update_save` branch, `ft = min(ft[ft > 0])` ([747](R/pipeline.R:747))
  collapses a per-year `ft` vector to its smallest positive value for the whole run,
  and if all `ft` are 0 returns `Inf` with a warning.
- `...` is forwarded only to `model_param_list_create()` ([218](R/pipeline.R:218)).
  Any misspelled `pipeline()` argument therefore becomes an unknown argument to that
  function and errors with a message about the wrong function; and because `args`
  is built from `environment()` before the `...` are consumed, the saved `call` is fine
  but there is no check that the dots are all valid epi parameters.

**Seed handling.** See item 2 of the summary. In addition, `set.seed(seed)` at
[239](R/pipeline.R:239) and [344](R/pipeline.R:344) with the comment "reset seed here as there is some randomness
in equilibrium (Need?)": there is no randomness in `equilibrium_ss_create`, so this
is a no-op that should go.

**Return value.** `res` is a bare list; when `update_save = TRUE` it is a list of
per-period lists whose last element has a different structure from the others; when
`FALSE` it is the raw C++ return. Metadata is stashed in `attr(res, "meta")`, which
is lost by any `[` subsetting, `lapply`, or `saveRDS`-then-modify. A class
(`"magenta_sim"`) with a stored `params`, `schedule` and `meta` field, and `print`,
`summary` and `as.data.frame` methods, would make the output self-describing.

### 4.2 The parameter-list layer

This is the part you said felt right, and it is. Keep it, finish it:

- **Make the C++ contract explicit and checked in one place.** Each `*_list_create`
  should return a classed list (`structure(l, class = c("magenta_drug_list", "list"))`),
  and a matching `validate_*()` should check names, types and lengths *against what
  C++ reads*. Today `drug_list_create` does no checking of `drugs` at all (there is a
  `TODO` to that effect at [pipeline_param_lists.R:72](R/pipeline_param_lists.R:72)),
  `barcode_list_create` checks one length, and the equilibrium list is checked by
  `identical(names(...))` inside `param_list_simulation_init_create`. A missing field
  currently surfaces as an Rcpp index error.
- **Decide what is public.** Everything a user can pass to `pipeline()` should be
  exported and documented together: `drug_list_create`, `drug_create`,
  `drug_create_al/asaq/dhappq/default_no_resistance`, `housekeeping_list_create`,
  `nmf_list_create`, `vector_adaptation_list_create`, `barcode_list_create`. The raw
  `Simulation_*_cpp` wrappers, `simulation_R`, the five `param_list_simulation_*`
  constructors, `population_get_genetics_*` and `%>%` should not be exported.
- **Hidden vector semantics.** `drug_list_create(resistance_flag=)` and
  `barcode_list_create(mutation_flag=)` are documented as Booleans but the tests pass
  per-year vectors; `pipeline()` then quietly takes `[1]` for C++ and cycles through
  the rest in the loop ([257-260](R/pipeline.R:257), [292-296](R/pipeline.R:292)).
  That is a schedule, not a parameter, and belongs with `ft`, `itn_cov`, `plaf` etc.
  (section 6).
- **Type inconsistency.** `drug_create()` returns a list;
  `drug_create_default_no_resistance()` returns a one-row `data.frame`
  ([pipeline_param_lists.R:339](R/pipeline_param_lists.R:339)). C++ reads both by name
  so it works, but `drug_list$drugs[[1]]$lpf <- ...` behaves differently on each.
- **`barcode_list_create` defaults** hard-code 24 (`plaf = rep(0.5, 24)`,
  `prob_crossover = rep(0.5, 24)`) rather than `num_loci`, so a direct call with
  `num_loci = 6` produces length-24 vectors that C++ will read past.
- **`model_param_list_create`** ([model_param_list_create.R:84](R/model_param_list_create.R:84)):
  72 arguments, then 115 lines of `mp.list$x <- x`. The two half-life arguments are
  overwritten at [281-282](R/model_param_list_create.R:281) so cannot be changed by a
  user. This can be `mpl <- as.list(environment())` followed by the ~12 derived
  quantities, with a validator that every element is a finite scalar.
- **`simulation_R`** ([mainFunctions.R:249](R/mainFunctions.R:249)) dispatches on the
  *shape* of the list through five sequential `if` blocks, each re-checking
  `is.list`, length and exact names. Replace with five plainly named internal
  functions (`sim_init(params)`, `sim_update(state, period)`, `sim_get_state(state)`,
  `sim_free(state)`, `sim_restore(saved)`) that call the Rcpp wrapper directly. The
  one-element "unpacked list" special case at [271](R/mainFunctions.R:271) is for a
  cluster use that no longer exists.

### 4.3 State saving and reload

- The saved object is the raw C++ dump (`population_List`, `scourge_List`, ...) with
  `saved_inputs` bolted on afterwards in three places. On reload, `pipeline()` reads
  `saved_inputs`, then overwrites half of them (4.1), then recomputes the equilibrium
  it does not need, then deletes `saved_inputs` and calls `gc()`.
- There is no version field in the saved object (the 1.3.7 fix for Windows bitset
  sizing shows why one is needed), and no single function a user can call to save or
  load; the test does `saveRDS(final_state, tf)` by hand.
- Memory ownership is implicit: `housekeeping_list$clear_up` decides whether the
  external pointer is freed. If `FALSE`, the pointer lives in `attr(res, "meta")$Ptr`
  and in `res[[i]]$Ptr` for each period, and nothing documents that the user must call
  the finaliser. `return_sim_memory` ([pipeline.R:819](R/pipeline.R:819)) walks every
  element removing `Ptr`.
- Recommended API: `sim_save(run, path)` writes `list(magenta_version, params,
  schedule, period_index, cpp_state)`; `sim_restore(path)` validates the version and
  returns a run object the loop can continue. `pipeline(saved_state_path=)` becomes a
  wrapper over that.
- `update_saved_state_barcode_plaf` and its helpers
  ([stat_genetics_utils.R:215-447](R/stat_genetics_utils.R:215)) are the model for the
  rest of the package: validate inputs with specific messages, small single-purpose
  helpers, tests for both the happy path and the error paths. The one thing it does
  not do is check or record the saved-state version.

### 4.4 The update loop and summaries

- `update_saves()` ([pipeline_utils.R:430](R/pipeline_utils.R:430)) takes 19
  arguments, 11 of them flags threaded through from `pipeline()`, and chooses what to
  store via a four-deep `if`. The magic condition `i %% 12 == 0 && i >= (length(res) - 180)`
  at [474](R/pipeline_utils.R:474) ("keep barcodes every 12th period in the last 180")
  is undocumented and only correct for `update_length = 30`.
- `COI_df_create()` ([universe_utils.R:92](R/universe_utils.R:92)) builds seven
  summary tables by copy-paste (`summary_unique1/2/3`, `summary_cou1/2/3`) and
  handles `n` being a scalar, `Inf`, or a vector with three different sampling
  schemes. It is the main genetics output used in the vignettes and has no test.
- `pop_alf()` and `lineages()` test `class(nums %>% unlist) == "raw"` to decide whether
  there is any data ([pipeline_utils.R:400](R/pipeline_utils.R:400)); `vapply`/`is.raw`
  would be clearer and safe against `class()` returning a vector.
- `drug_list_update()` ([pipeline_utils.R:361](R/pipeline_utils.R:361)) is the only
  place state mutates between periods on the R side; it `message()`s the drug choice
  unconditionally and is at 29% coverage (sequential and temporal cycling untested).
- `progress_logging()` ([utils.R:13](R/utils.R:13)): with `cluster = TRUE`, `quiet = TRUE`
  *prints* progress markers and `quiet = FALSE` prints nothing. The name is inverted
  relative to its effect. The C++ side uses the same flag with the opposite sense.

### 4.5 Research code that should leave the package

Zero references from `R/`, `tests/` or `vignettes_src/`; 0% coverage; last edited
2019-2020:

| File | Lines | Why it should go |
|---|---|---|
| `bootstrap.R` | 121 | Power analysis for a specific figure; default `sim = "M:/OJ/magenta_Results/scripts/fig3_r1.rds"`. |
| `kl.R` | 92 | KL divergence for one paper; reads an RDS by path. |
| `cluster_utils.R` | 236 | Pools didehpc task results (`grp$tasks[[i]]$status()`); didehpc is not a dependency. |
| `plotting.R` (most of it) | ~700 | `summary_data_frames_from_sims` hard-codes `Low/Medium/High` intervention labels and `rep(c(rep(1,10),rep(2,10),rep(3,10)))`; `plot_x_loggers`/`plot_mean_summary` take an `extra` string and `eval(parse(text = extra))` it; `heatmap_ordered_binary_plot` hard-codes 24 loci and "Times New Roman". `microscopy_detected` has the `age_brackets(1)` bug. |
| `utils.R`: `person_make`, `binary`, `func_call`, `assert_*` | ~50 | Unused. The `assert_*` helpers are good and should be *used*, not deleted. |
| `stat_genetics_utils.R`: `population_ibd_barcodes`, `population_ibd_barcodes_c`, `pibd_from_barcode_list`, `summarySE_mean_only_max_mean` | ~50 | Unused. |
| `seasonal_utils.R`: `country_seasonal_profiles` | 33 | Plotting helper, `par(mfrow)` side effects. |

Dropping `plotting.R`, `bootstrap.R`, `kl.R` and `cluster_utils.R` removes `ggplot2`,
`reshape2`, `cowplot`, `grDevices::tiff` and most of the `stats` imports from the core.
If the plots are still wanted, they belong in a companion analysis repository (or a
`magenta.plots` package) that depends on magenta, where they can be rewritten against
the classed result object.

`convert_barcode_vectors` (exported, used once in the Introduction vignette) and
`sample_coi` / `coi_age_plot_sample_x` are the borderline cases: keep
`convert_barcode_vectors` if the vignette keeps using it, but it needs a test.

### 4.6 Correctness issues found while reading

These are things to fix, each small, each wanting a regression test. None were
hunted for; they are what fell out of reading the code once.

| Where | Issue |
|---|---|
| [model_param_list_create.R:281-282](R/model_param_list_create.R:281) | `irs_half_life` and `itn_half_life` arguments are overwritten with the defaults; the arguments have no effect. |
| [mainFunctions.R:288,317,345,371,398](R/mainFunctions.R:288) | `set.seed(seed)` before every C++ call; C++ consumes R's RNG, so every update chunk replays the same stream. |
| [pipeline.R:695](R/pipeline.R:695) | `i` undefined when there is exactly one update period. |
| [pipeline.R:202-210](R/pipeline.R:202) | `use_historic_interventions` with a year-range `years` (as documented) errors downstream. |
| [pipeline.R:747](R/pipeline.R:747) | `min(ft[ft > 0])` collapses a per-year `ft`; `Inf` if all zero. |
| [pipeline.R:777](R/pipeline.R:777) | `human_vars` redefined without `"ID"` in one of three copies. |
| [pipeline.R:264](R/pipeline.R:264) | `spatial_type = "metapop"` is never mapped to `2`; fails with an Rcpp error instead of the intended `stop()`. |
| [pipeline_utils.R:142-145](R/pipeline_utils.R:142) | `plaf_matrix_check` pads with `matrix[2, ]` (row 2, not the last row); errors on a one-row matrix. |
| [pipeline_utils.R:151](R/pipeline_utils.R:151) | Vector `plaf` with non-integer `years`: `rep(matrix, years - 1)` truncates, then `matrix(nrow = ceiling(years))` recycles. |
| [pipeline_utils.R:114](R/pipeline_utils.R:114) vs [142](R/pipeline_utils.R:142) | `spl_matrix_check` *prepends* copies of row 1; `plaf_matrix_check` *appends* copies of row 2. Two different padding semantics for the same concept, neither documented. |
| [pipeline_utils.R:191-195](R/pipeline_utils.R:191) | `is.null(options("quiet_admin"))` is never `TRUE` (`options()` returns a list); `quiet = options("quiet_admin")` passes a list; `admin_match()` ignores `quiet` anyway and always messages. |
| [pfpr_to_eir.R:51-52](R/pfpr_to_eir.R:51) | `eqInit$age` is used before `eqInit` exists; `age_range` can never work (lintr flags it as a missing binding). |
| [plotting.R:824](R/plotting.R:824) | `age_brackets(1)` calls the function instead of indexing the argument. |
| [pipeline_param_lists.R:454-455](R/pipeline_param_lists.R:454) | `barcode_list_create` defaults are length 24 regardless of `num_loci`. |
| [magenta.R:27](R/magenta.R:27) | `globalVariables(c(..., ""))` contains an empty string. |
| [pipeline.R:239](R/pipeline.R:239) | `set.seed()` "reset" before a deterministic step; no-op. |

### 4.7 Style

Concrete conventions the code currently mixes, with the recommended choice (this is
the tidyverse style guide, which is what `styler` and `lintr` default to; the point is
to pick one and automate it, not which one):

| Dimension | Seen in the code | Pick |
|---|---|---|
| Object names | `sim.out`, `mp.list`, `age.vector`, `chosen.country` / `statePtr`, `eqSS`, `eqInit`, `EIRY_eq` / `drug_list`, `num_loci` | `snake_case` everywhere; rename the handful of dot and camel names |
| Assignment | `=` and `<-` both used for assignment (`output = sum(...)`, `age_rate[na] = 0`, `IB_eq = array(...)`) | `<-` |
| Logical constants | `T`, `F` (`byrow = T`, `header = T`) | `TRUE`, `FALSE` |
| Scalar logic | `&`, `\|` in `if()` conditions (`is.null(x) & !is.null(y)`) | `&&`, `\|\|` |
| Class checks | `class(x) == "list"`, `class(x) != "externalptr"` | `is.list()`, `inherits()` |
| Sequences | `1:length(x)`, `1:nrow(x)` | `seq_along()`, `seq_len()` |
| Braces | `}\n else {`, `if (x)\n {`, `function(x){` | `} else {`, `function(x) {` |
| Return | explicit `return()` in most functions, implicit in others | implicit (last value), `return()` only for early exit |
| Pipes | `x %>% unlist`, `x %>% unlist()`, nested `lapply`s | drop `magrittr`; use `\|>` (needs `R >= 4.1`) or plain calls |
| Comments | section banners `# ---- ####`, `#---`, `#----------`, `###` | one style; roxygen `@section` for docs |
| Namespacing | `dplyr::group_by` in some places, `importFrom` in others | `pkg::fun()` for anything outside base |
| Validation messages | `stop("eqSS is not of class list")`, `stop("mu_vec not long enough")`, `stop(sprintf("...%s...", name))` | one `assert_*` family with the argument name in the message (the one in `utils.R` is a good start) |
| roxygen | `@param` descriptions like "Default = 0.4" repeated in prose; `@keywords internal` on things that are exported; `@noRd` on things that are not | `@keywords internal` + `@noRd` only for truly internal; defaults live in `\usage` not prose |

Mechanically: run `styler::style_pkg()` once as its own commit, add a `.lintr` file
with the project's chosen linters, and add `lintr::lint_package()` to the
`check-standard` workflow so it cannot drift again.

### 4.8 Package hygiene

| Item | Now | Recommend |
|---|---|---|
| `Depends: R (>= 2.10)` | 2010-era | `R (>= 4.1)` (native pipe, `\(x)`) |
| `SystemRequirements: C++11` | check NOTE ("support has been removed") | delete the field; `src/Makevars` needs no `CXX_STD` |
| `Imports: redux` | only in commented-out code | remove |
| `Imports: dde` | never referenced from R (odin uses it internally when `use_dde = TRUE`) | keep only if `mu_fv_create` really needs `use_dde = TRUE`; then `@importFrom dde dopri` to make it explicit, otherwise drop |
| `Imports: ggplot2, reshape2, magrittr, dplyr, rlang` | plotting and `COI_df_create` only | `ggplot2`/`reshape2` go with the plotting code; `dplyr`/`rlang` can be replaced by base `aggregate`/`tapply` in `COI_df_create` or kept; `magrittr` goes once `\|>` is used |
| `Remotes: mrc-ide/odin, mrc-ide/dde` | odin is on CRAN | remove `odin`; `dde` is on CRAN too |
| `Maintainer:` | "Oliver Watson <o.watson@imperial.ac.uk>" conflicts with `Authors@R` ("OJ Watson", same email) | delete the `Maintainer` line; `Authors@R` generates it |
| `README.Rmd` | "install `rdhs`", Travis badge, WIP badge, mentions a `Loggers` vignette that does not exist | rewrite; add GitHub Actions and codecov badges (the workflows already exist) |
| Examples | none in any Rd (`checking examples ... NONE`) | `@examples` on every exported function, `\donttest{}` the slow ones |
| `vignettes/` | prebuilt `.html` and `*_files/` committed, figure names with spaces (check WARNING) | keep the `vignettes_src/` precompute pattern but fix chunk labels (no spaces) and ignore the `.html` |
| `codemeta.json` | version 1.3.0 | regenerate or delete |
| `src/` | 51 MB on disk: `.o`/`.so` from local builds (gitignored, but `R CMD build` from this directory picks them up) | `R CMD build` from a clean clone or add `^src/.*\.(o\|so)$` to `.Rbuildignore` |
| `.github/.gitignore` | stray untracked file containing `*.html` | delete or move to root `.gitignore` |
| `NAMESPACE` / `man/magenta.Rd` | modified but uncommitted (roxygen 7.3 format change) | commit |
| `tests/testthat.R` | fine; but no `Config/testthat/edition: 3` and tests use `context()` | move to edition 3; drop `context()` |
| NEWS.md | good, keep | |
| CI | `check-standard`, `test-coverage`, `pkgdown` present | add lintr; add `NOT_CRAN=true` explicitly in the coverage job if r-lib's action ever stops setting it (today it does, which is why the saved-state tests run on CI but not in a bare `covr::package_coverage()`) |

---

## 5. Test coverage

### 5.1 What exists

Five files, 14 `test_that()` blocks, 661 lines, all passing. Four of the fourteen are full
multi-year `pipeline()` runs with N = 1000 (the "simulations" and "seasonality" tests),
two are saved-state round-trips with N = 200, two test
`update_saved_state_barcode_plaf`, and five test C++ helpers through the `test_*`
Rcpp exports (barcode generation, IBD recombination, bitset serialisation).


### 5.2 Coverage, R side, by function (from `covr`, `NOT_CRAN=true`)

Fully covered (100%): `model_param_list_create`, `drug_create`, `drug_create_asaq`,
`drug_create_default_no_resistance`, `seasonal_profile`, `equilibrium_ss_create`,
`spl_create`, `housekeeping_list_create`, `nmf_list_create`,
`vector_adaptation_list_create`, `lineages`, `ztrnbinom`, `rebuild_barcode_vectors`,
`set_nested_element`, `get_nested_element`, `ranges`, `fd`, `q_fun`.

Mostly covered (85-97%): `equilibrium_init_create` 97%, `drug_list_create` 96%,
`update_saved_state_barcode_plaf` 95%, `progress_logging` 94%, `mu_fv_create` 93%,
`return_sim_memory` 92%, `barcode_list_create` 89%, `match_clean` 89%,
`flatten_barcode_vectors` 87%, `simulation_R` 86%, `pop_alf` 86%, `pipeline` 85%,
`pop_strains_df` 84%.

Partially covered (29-83%): `age_brackets` 83%, `param_list_simulation_init_create`
81%, `spl_grab` 77%, `admin_match` 68%, `pfpr_to_eir_heuristic` 68%,
`intervention_grab` 65%, `param_list_simulation_update_create` 65%,
`adjust_matrix_to_plaf` 66%, `plaf_matrix_check` 60%, `update_saves` 56%,
`spl_matrix_check` 55%, `drug_list_update` 29%.

Zero coverage (41 functions): everything in `bootstrap.R`, `kl.R`, `cluster_utils.R`;
all of `plotting.R` except the two one-liners; `COI_df_create` (160 instrumented
lines); `generate_default_model`; `drug_create_al`, `drug_create_dhappq`;
`summarySE`, `summarySE_mean_only`, `summarySE_mean_only_max_mean`;
`clonality_from_barcode_list`, `cou_from_barcode_list`, `convert_ibd_barcode`,
`population_ibd_barcodes*`, `pibd_from_barcode_list`; `rbind_list_base`,
`person_make`, `binary`, `func_call`, the `assert_*` family;
`country_seasonal_profiles`.

### 5.3 What the numbers hide

The 85% on `pipeline()` is flattering: it is reached only through four end-to-end
runs, so the uncovered 15% is exactly the branches nobody set up:

- the non-`update_save` path ([743-790](R/pipeline.R:743)), which is what
  `pipeline(EIR = 1, years = 10, N = 10000)` in the Introduction vignette uses;
- `use_historic_interventions` ([204-210](R/pipeline.R:204));
- `human_only_full_save` on resume ([442-449](R/pipeline.R:442)).

Nothing tests the *values* the pipeline produces. The seasonality test asserts that
variance with seasonality exceeds variance without; the fitness test asserts one
allele-frequency sum exceeds another; both with `seed = as.integer(runif(1, 1, 1e9))`,
so both can fail by chance and neither can detect a regression that preserves the
inequality. There is no fixed-seed snapshot of any output, so a refactor of
`pipeline()` today cannot be verified as behaviour-preserving.

No pure function has a unit test: `age_brackets`, `spl_matrix_check`,
`plaf_matrix_check`, `intervention_grab` (year-range edge cases), `mu_fv_create`
(`length_check`), `drug_list_update` (cycling logic), `pop_alf`, `lineages`,
`bitsToInt`, `ztrgeomintp`, `match_clean`/`admin_match` (fuzzy matching), the
`param_list_*` validators' error paths, `simulation_R`'s "no handler" error. These are
all fast and all deterministic; they are where the bugs in 4.6 live.

`COI_df_create` is the biggest single gap: 248 lines, 0% coverage, the function that
turns a population into the COI/clonality/COU tables that the vignettes and papers
report.

Test code quality: `test-simulations.R` contains three ~95-line `pipeline()` calls
that differ in two or three arguments; `test-reload.R` defines its helper inside the
test file rather than in `helper-*.R`; `context()` is deprecated in testthat 3e.

### 5.4 Test timings

`NOT_CRAN=true`, local run via `pkgload::load_all()` + `testthat::test_dir()`, 19.6 s
in total. The suite is already fast; the point of the unit tests in 5.5 is coverage
and determinism, not speed.

| File | Test | Expectations | Seconds |
|---|---|---|---|
| test-simulations.R | Seasonality Test | 1 | 14.6 |
| test-simulations.R | Fitness Test | 1 | 2.1 |
| test-reload.R | saved state can be reloaded for continuation or inspection | 7 | 0.9 |
| test-simulations.R | Asymptomatic Prophylaxis Test | 1 | 0.6 |
| test-simulations.R | dependent plaf draws with 100% importation fixes plaf | 2 | 0.6 |
| test-strain.R | random from plaf | 2 | 0.3 |
| test-reload.R | pipeline full saves contain restart metadata | 3 | 0.2 |
| test-plaf-update.R | update_saved_state_barcode_plaf adjusts barcode PLAF | 5 | 0.06 |
| test-util.R | bitset roundtrip | 40 | 0.03 |
| test-plaf-update.R | validates inputs | 2 | 0.02 |
| test-strain.R | depedent plaf, ibd recombinants, ibd next generation, ibd conversion | 4 | <0.01 |

The seasonality test is 75% of the runtime because it runs two 6-year, N = 1000
simulations with the odin ODE for interventions; it would be the natural candidate
for `skip_on_cran()` once there are fast unit tests behind it.

### 5.5 Priority order for new tests

1. **Characterisation snapshots first** (before any refactor): four or five
   `pipeline()` configurations at N = 100-200, 1-2 years, fixed seed, with
   `testthat::expect_snapshot_value()` on a digest of the output (or key summaries).
   These exist to prove a refactor changed nothing. Note that fixing the seed
   replay (item 2) will legitimately change them; do that fix as its own PR and
   re-snapshot.
2. Unit tests for the schedule/matrix helpers (`spl_matrix_check`,
   `plaf_matrix_check`, `intervention_grab`, `spl_grab`, `mu_fv_create$length_check`),
   including the edge cases in 4.6.
3. `COI_df_create` and `pop_strains_df` on a small synthetic data frame (no C++ needed).
4. `drug_list_update` cycling: sequential and temporal, boundary at `number_of_drugs`.
5. Validator error paths for every `*_list_create` and `param_list_simulation_*`.
6. `equilibrium_init_create` numerical regression against stored values for two or
   three (EIR, ft) pairs; `pfpr_to_eir_heuristic` round-trip.
7. Make the stochastic tests deterministic (fixed seed) and replace inequality
   assertions with tolerance checks on known-good values.

Target: 80%+ on the R side with the research code removed. The suite already runs
in 20 s, so there is no speed excuse for not running it on every push.

---

## 6. Proposed target architecture for the R side

The aim is that `pipeline()` still exists with the same signature for existing
scripts, but is a 40-line wrapper over pieces that can each be read, tested and
replaced on their own.

### 6.1 Four layers

```
R/params-*.R        parameter objects:  new_*() constructors, validate_*(), print methods
R/schedule.R        build_schedule(): one row per update period, validated once
R/sim-bindings.R    sim_init(), sim_update(), sim_get_state(), sim_free(), sim_restore()
R/run.R             sim_setup(), sim_run(), sim_save(), sim_restore(); the loop
R/savers.R          save_loggers(), save_summary(), save_humans(), save_full()
R/results.R         magenta_sim class: print, summary, as.data.frame, meta accessors
R/equilibrium-*.R   unchanged maths, split so the 312-line function is 3-4 pieces
R/pipeline.R        pipeline(): back-compatible wrapper, ~40 lines, deprecating old flags
```

### 6.2 Parameters

`magenta_params(...)` returns a classed list holding the sub-lists exactly as C++
wants them (`core`, `barcode`, `drug`, `nmf`, `vector_adaptation`, `housekeeping`,
`spatial`, `equilibrium`). Each sub-list has `new_<x>_list()` (today's
`*_list_create`, renamed or aliased) and `validate_<x>_list()`, and the top-level
validator checks the set of names C++ reads. Put the list of field names C++ reads in
one R object (`.cpp_contract`) and test that each constructor produces exactly that
set, so the contract is visible from R and a change on either side fails a test.

`model_param_list_create()` becomes `as.list(environment())` plus derived values plus
`validate_core_params()`, and the half-life bug disappears with it.

### 6.3 Schedule

Everything that varies by period is one data frame built by
`build_schedule(years, update_length, ft, itn_cov, irs_cov, resistance_flag,
mutation_flag, plaf, spatial_incidence, spatial_mosquitoFOI, drug_choice, ...)`:

| period | day_start | day_end | year | ft | itn_cov | irs_cov | resistance_flag | mutation_flag | plaf (list col) | importation_human | importation_mosquito |
|---|---|---|---|---|---|---|---|---|---|---|---|

Built once, with one padding rule (recycle the *last* value forward, error if a vector
is longer than `years`), validated once (lengths, ranges, types). The main loop then
reads row `i` and never calls `rep`, `ceiling`, `%%` or a `*_matrix_check` again.
`mu_fv_create` consumes the schedule's `ft`/`itn`/`irs` columns. `spl_matrix_check`,
`plaf_matrix_check` and the ad-hoc flag expansion in `pipeline()` are deleted.

### 6.4 Run loop and savers

```r
sim_run <- function(params, schedule, saver = save_loggers(), seed, progress = TRUE) {
  state <- sim_setup(params, seed)                     # or sim_restore(path)
  out   <- vector("list", nrow(schedule))
  for (i in seq_len(nrow(schedule))) {
    period <- schedule[i, ]
    state  <- sim_update(state, period, mu_fv[period$day_start:period$day_end, ])
    out[[i]] <- saver(state, period, params)
    params$drug <- update_drug_choice(params$drug, period, out[[i]])
  }
  new_magenta_sim(out, params, schedule, meta = list(seed, version, call, timings))
}
```

A saver is a function `function(state, period, params) -> list`. The four current
output modes become four savers (`save_loggers()`, `save_summary(sample_size,
sample_states, age_breaks, barcodes_every = 12)`, `save_humans()`, `save_full()`),
and `update_save_func` becomes "pass your own saver". The seven booleans collapse to
one argument. `update_saves()` is deleted.

`set_up_only = TRUE` becomes `sim_setup()` plus `sim_get_state()`; the non-update
path becomes a schedule with one row.

### 6.5 State

`sim_save(sim, path)` and `sim_restore(path)` as in 4.3, with a `version` field and
the params and schedule stored alongside the C++ dump, so resuming does not need to
recompute anything and cannot silently mix saved and current inputs.

### 6.6 Seeds

Seed once in `sim_setup()`. If per-chunk reproducibility is wanted (so that a resumed
run matches an uninterrupted one), derive a per-period seed from the master seed and
the period index and store it in the schedule; never reuse the master seed.

### 6.7 Results

`magenta_sim` object with `$output` (list of per-period saver results),
`$params`, `$schedule`, `$meta`; `print()` shows the call, N, years, EIR, saver, seed,
version, elapsed; `as.data.frame()` for the loggers saver; `summary()` for the summary
saver. The attribute-based `meta` goes away.

---

## 7. Phased plan

Each phase is one or two PRs; each PR leaves `R CMD check` and the tests green.

Phase 6 starts only after `docs-revamp`, Phase 0 and `rewrite` are merged; it
is scoped to the five PRs listed under it and nothing else.

**Status (2026-11-13).** `rewrite` is rebased onto `master` at 88ffded (the
`docs-revamp` merge) and carries Phase 0 and everything after it as one
branch for one PR; it supersedes `worktree-r-code-review`. The characterisation
snapshots moved once more with the three 1.3.8 model fixes and were accepted.

**Status (2026-11-12).** Phase 0 is done on branch `worktree-r-code-review`, cut
from the `docs-revamp` tip so it merges after the documentation work: snapshot
tests (`tests/testthat/test-characterisation.R`, Linux only), styler + `.lintr` +
a lint workflow (zero lints), research code removed, DESCRIPTION brought up to
date, the plain bug fixes in 4.6 with tests, version 1.4.0 with NEWS.

Phases 1-3 are done on branch `rewrite`, cut from the
Phase 0 tip, as one commit per step with the tests and lints green after each.
The characterisation snapshot file is unchanged throughout; the tests read
`out$output` instead of `out` because of the new return object. What landed,
and where it departs from the plan above:

- *Phase 1.* `R/sim_bindings.R` holds `sim_init()`, `sim_update()`,
  `sim_get_state()`, `sim_free()` and `sim_init_saved()`; the last is named
  that rather than `sim_restore()` because 6.5 reserves `sim_restore(path)` for
  the user-facing function. `simulation_R()` is a deprecated, unexported
  alias. `pipeline()` is a wrapper over `sim_setup()`, `run_schedule()`,
  `prepare_state_save()` and `collect_output()` (~120 lines of body).
- *Phase 2.* `R/params_contract.R` holds `cpp_contract` and the `check_*()`
  helpers; every constructor returns a classed list and has a `validate_*()`;
  `tests/testthat/test-params.R` checks the constructors against the contract
  and the contract against the C++ sources. `model_param_list_create()` is
  `as.list(environment())` plus derived values, reordered to the contract so
  the list is unchanged. `R/schedule.R` holds `build_schedule()` with the
  padding rule of 6.3 (last value forward, longer than `ceiling(years)` is an
  error); the behaviour changes it implies for odd inputs are in NEWS and none
  touched a snapshot. The year change that fell on the final period, the
  single-period `update_save` error and calendar-year `years` were fixed in
  passing because the schedule made the old code paths unreachable; they are
  in NEWS too. `Person::m_cA` is initialised in C++ because an uninitialised
  read showed up in the set_up_only snapshot once allocations changed.
  Exports: constructors and `build_schedule()` in, Rcpp wrappers,
  `simulation_R()`, `param_list_simulation_*()` and `%>%` out.
- *Phase 3.* `R/savers.R` (`save_loggers()`, `save_summary()`,
  `save_humans()`, `save_full()`, each with `final = "loggers" | "full" |
  "humans"` for the last slot; a user `function(state, period, params)` is
  accepted; the old flags map onto a saver with a warning naming the
  replacement; `update_saves()` deleted). `R/results.R` (`magenta_sim` with
  `print()`, `summary()`, `as.data.frame()`; `sim_save()` / `sim_restore()`
  writing a flat, classed, `format_version`ed object that is also accepted
  by `pipeline(saved_state_path = )`, as are the old plain full saves). The
  saver's `final` replaces `full_save` / `human_only_full_save` for
  `set_up_only` too. The vignette sources are rewritten against the new
  interface.

Phase 4 is applied, after OJ's review of each item: one seed per run with
the generator state kept in `meta$rng_state` and in full saves and restored
on resume; one period per year when `update_save = FALSE` so that every
yearly input applies (`itn_cov`/`irs_cov` already did, through the odin
solve); a resumed run continues with all of the saved run's lists and errors
on a differing explicit argument; `spl_grab()` keeps the year 2000. Every
stochastic characterisation snapshot moved with the seed change and was
accepted deliberately; the equilibrium snapshot did not.
`tests/testthat/test-findings.R` now tests the fixes. Sequential drug
cycling reads the per-period loggers (`period_treatment_failure()`), so it
works with every saver and without periodic saving; the `save_summary()`
path switches at the same times as before (pinned).

Deferred, to take with the next snapshot-moving change: the no-intervention
shortcut in `mu_fv_create()` uses `fv = 0.333`, `mu = 0.132` where the
closed form gives `1/(tau1 + tau2)` and `-fv * log(p10 * p2)` exactly;
removing it moves every no-intervention run by a rounding error. The other 4.6 items are done (Phase 0, or absorbed in
Phases 2-3 as noted above). Phase 5: unit tests for every `sim_*`, parameter,
schedule, saver and `magenta_sim` function and for the utilities (R-side
line coverage 86.9% -> 98.7% over this step; 47.6% at Phase 0), an example
on every export, `R CMD check --as-cran` clean (one pre-existing NOTE).
`drug_list_create(number_of_resistance_loci = 1)` no longer crashes (a
`drop = FALSE`; the table for two or more loci is unchanged).
`generate_default_model()` is deleted.

*The deterministic model.* The pipeline only ever read `mu` and `fv` from
the odin run, and both are closed-form in `t` (the intervention block of
`inst/extdata/odin_itn_irs.R`: ITN/IRS decay, the per-category feeding and
repellency probabilities, `zbar`, `wbar`, `p1`); no state feeds into them.
`R/mu_fv.R` evaluates that block directly in the generated C's order of
operations, and `tests/testthat/test-mu-fv.R` shows it bit-identical to the
model (`expect_identical`) for ITN only, IRS only, both, a five-year ramp, a
part year, off-anniversary change times and the Kinshasa 2012-2013 seasonal
run; the characterisation snapshots did not move. The `t < ITN_IRS_on`
branch is implemented but cannot be reached from the pipeline, which fixes
`ITN_IRS_on = 0` and starts output on day 1, so it is not covered by the
comparison. The rates for that seasonal run take 0.5 ms instead of 2.0 s.
The model stays in the package: `run_deterministic_model()` (exported) runs
it from equilibrium for an EIR, treatment rate, age and heterogeneity setup
and returns every output, which is how to check the individual-based
simulation against it (prevalence and incidence for a large population with
no resistance; its example shows the comparison). `R/odin.R` was generated
with odin 1.5.10 (r-universe) while CRAN odin is 1.2.7; regenerate with a
pinned version later, not now, and do not migrate to odin2. The model is
written against odin 1, is an early copy of the Imperial deterministic
malaria model, and is a candidate for replacement by the current MRC-IDE
packages (under review separately).

Still open: the deferred `mu_fv_create()` shortcut above, and a coverage
gate in CI (the pkgdown reference regrouping is done).

*Parked for Phase 6.* The research-era helpers that Phase 0 removed are back
on the branch, parked and unchanged (styled; ggplot2 calls qualified;
`requireNamespace()` guards): `R/cluster_utils.R` (`pool_res()`,
`create_res_list()`), `R/plotting.R` (the `plot_*()` functions,
`heatmap_ordered_binary_plot()`, `sample_coi()`, `coi_age_plot_sample_x()`,
`summary_data_frames_from_sims()`, `gg_color_hue()`) and
`country_seasonal_profiles()` in `R/seasonal_utils.R`. `bootstrap()` and
`kl()` stay removed. They are written against the 1.3.x return value of
`pipeline()`, are excluded from lintr and are left out of the coverage figure:
R-side line coverage is 98.6% without `R/cluster_utils.R` and `R/plotting.R`
and 78.3% with them (whole package 79.9% and 75.9%).

**Phase 6: workflows (not started).** Give the six recurring steps (Appendix B)
package functions, written against `magenta_sim`, and three vignettes that show
them end to end. ggplot2 and hipercow go in `Suggests`; nothing in `Imports`
changes. Each PR leaves check and the characterisation snapshots green. The
parked helpers in `R/cluster_utils.R`, `R/plotting.R` and
`country_seasonal_profiles()` are the raw material and are deleted by the PR
that replaces each. The collectors in 6.1 are checked against three code bases:
arms, peru_hrp23 and the Lancet Microbe 2022 ("MMC") analysis, whose research
compendium is `compendia/art_resistance_consensus_modelling` (local; its
`analysis/` reads pre-summarised CSVs and `R/utils.R`, `R/plot_utils.R` hold only
`rbind_list_base()` and `save_figs()`; the pooling from magenta runs happened in
the MMC_WP2 work-package repository, which the summaries' columns describe).

**6.1 Tidy collectors (one PR, pure R, no new dependencies).**
`R/tidy.R`. One function per quantity, each taking a `magenta_sim` and returning
a data frame with `period`, `day`, `year` (from `out$schedule`) first:

- `tidy_periods(sim)` — what `as.data.frame.magenta_sim()` does now, kept as the
  alias so the family has one name.
- `tidy_prevalence(sim, age_range = NULL)` — `S, D, A, U, T, P`, `pcr_prev`
  (`D+A+U+T`), `dat`, `micro_2_10` (via `microscopy_detected()`), from loggers
  or from a `save_summary()` slot; replaces `plot_prevalence()`'s two branches
  and `pool_res`'s `make_prevs_df`.
- `tidy_treatment(sim)` — successful, LPF, not-treated counts, `treatment_failure`,
  `overall_treatment_failure`, `cases = successful + lpf + not_treated` (arms
  `get_df`, peru `get_cases`, by name not position).
- `tidy_allele_freqs(sim)` — long (period, locus, af) from `af`; `tidy_mutations()`
  likewise from `mutations`.
- `tidy_lineages(sim, genotypes = NULL)` — (period, lineage as a `"011"`
  barcode string, count, freq) from `lineage`; the `intToBits` labelling from
  `plot_resistance_lineages`. With `genotypes = magenta::drug_table` the
  six-locus lineages are labelled with the `Genotype` names (`"TNY0C0"`), which
  is what the MMC `lineage_df` (`Strains, Time, rep, Genotype`) carried.
- `tidy_mutation_events(sim)` — mutation events per locus over the run (the
  MMC `NUMBER_OF_MUTATION_EVENTS_1..6` columns) from `mutations`, and the run
  totals of treatments and treatment failures (`NUMBER_OF_TREATMENTS`,
  `NUMBER_OF_TREATMENT_FAILURES`) as sums of `tidy_treatment()`.
- `time_to(df, var, threshold)` — the first `year` at which a tidy column
  crosses a threshold, `NA` if never: the MMC `TIME_TO_1p_580Y`,
  `TIME_TO_10p_TF` ... milestones, computed once from `tidy_allele_freqs()` and
  `tidy_treatment()` rather than per project.
- `tidy_coi(sim, detected = FALSE)` — (period, metric in {COI, COU, % polygenomic,
  % unique, pIBD, iIBD}, group_by, group, mean, N) from the `summary_*` tables.
- `selection_coefficient(af, t, window = c(0.05, 0.975))` — arms `get_s()`, as a
  plain function on two vectors.
- `bind_sims(sims, grid = NULL)` — `rbind` of any `tidy_*` over a list of sims,
  adding `run`, `seed` (from `meta$seed`) and the grid's columns. This is
  `pool_res`'s outer loop, `create_case_df` and 02_sim_processing's `cbind(res, parms[i,])`.

Tests: run the four savers once in `helper-sim.R` (already there), check column
names and row counts per saver, check `tidy_prevalence` against the loggers by
hand, check `selection_coefficient` recovers a known `s` from a logistic curve,
and that `bind_sims` carries grid columns. Delete `cluster_utils.R` in this PR;
`pool_res` has no caller once `tidy_*` exist.

**6.2 Grid and local batch (one PR).** `R/grid.R`, `R/batch.R`.

- `sim_grid(..., reps = 1, seed = NULL, args = history_policy_args)` —
  `expand.grid` over whatever columns the user likes (scalars, or list-columns
  wrapped with `I()` for a `plaf` matrix, a `drug_list`, a yearly `ft` vector),
  then `rep` and a deterministic `seed` column (`seed + row`), class
  `magenta_grid`. The grid never guesses how its columns map onto `pipeline()`:
  the user supplies an args function that turns one row into the argument
  list, kept with the grid, and `grid_args(grid, i)` calls it on row `i`, so
  the worker function is always `do.call(pipeline, grid_args(grid, i))`,
  locally and on the cluster. The package ships one default,
  `history_policy_args()`, for the common history-then-policy shape (columns
  `years_history`, `years_policy`, `ft`, `EIR`, `seed`, plus the policy's
  `drug_list` and starting `plaf`): it pads `resistance_flag`, `mutation_flag`,
  `plaf` and the importation vectors to history years then policy years, so
  the standard case (arms, the MMC runs) is a one-liner; anything beyond that
  shape is the user's own function. This replaces the `for` loops in
  `01_dideweb.R` and `didewin.R` and the `split(pardf, pardf$seed)` in peru.
  Because 1.4.0 validates every parameter list, a bad row fails at
  `grid_args()` time, locally, not on the cluster.
- `run_grid(grid, dir = NULL, parallel = c("none", "mclapply"), progress = TRUE)`
  — runs every row locally, returns a list of `magenta_sim` or, with `dir`,
  writes `dir/run_<row>.rds` and returns the paths; skips rows whose file exists
  (restartable). `read_grid(dir)` reads them back and `bind_sims()` consumes
  either. `parallel::mclapply` is base R; no `future` dependency.
- `shrink_sim(sim, drop_final = TRUE)` — drops the population-length vectors
  (`InfectionStates, Ages, IB, ICA, ICM, ID`, barcode lists) from the final slot,
  which both arms and peru did by hand before saving. A plain helper the worker
  function can call, not a saver option: the default stores everything, and
  what to drop is decided afterwards (decision 3).

Tests: `sim_grid` is pure (row counts, list-columns, seeds unique and
reproducible); `run_grid` on a two-row grid with `N = 100, years = 1`; `read_grid`
round trip; `shrink_sim` sizes.

**6.3 hipercow runner and the cluster vignette (one PR).**

- `run_grid(grid, backend = "hipercow", bundle_name = , resources = , environment = )`
  — `hipercow::task_create_bulk_call(function(i) do.call(magenta::pipeline, grid_args(grid, i)), seq_len(nrow(grid)), ...)`,
  saves `grid` as `hipercow/magenta/<bundle_name>_grid.rds` next to the bundle
  name, returns the bundle. Results stay in the hipercow store (decision 4):
  `collect_grid(bundle_name, tidy = list(tidy_prevalence, tidy_treatment, ...), dir = NULL)`
  loads the bundle by name, reads each task's result with `task_result()`,
  applies the `tidy_*` collectors and `bind_sims()` with the saved grid, and
  returns the tidy frames; per-row `run_<row>.rds` files are written only when
  `dir` is given. Two-stage runs are two grids: the second has a
  `saved_state_path` column pointing at the first's files (or `sim_save()`
  objects), which `grid_args()` passes through unchanged (`pipeline()` accepts
  the path).
- Vignette `vignettes/cluster.Rmd`, "Running magenta in batches": local
  `run_grid()` first, then the hipercow section (`hipercow_init`, a
  `pkgdepends.txt` line `github::OJWatson/magenta@v1.4.0`, `hipercow_provision()`,
  `run_grid(backend = "hipercow")`, status, retry, `collect_grid()`), then the
  setup-then-continue pattern from peru_hrp23 with `sim_save()` and
  `update_saved_state_barcode_plaf()`. Chunks that touch the DIDE driver are
  `eval = FALSE`; the rest run under hipercow's `"example"` driver, which runs
  tasks locally in the same layout (decision 7).
- `DESCRIPTION`: `Suggests: hipercow`, `Additional_repositories: https://mrc-ide.r-universe.dev`;
  the `check-standard` and `test-coverage` workflows add the r-universe repo to
  `setup-r-dependencies` so hipercow installs in CI. Before writing the test,
  read hipercow's own CI and vignette setup (how its vignettes run under the
  example driver, what it mocks) and record what it does here.

Tests: `skip_if_not_installed("hipercow")`, `hipercow_init(driver = "example")`
in a temp dir, two-row grid, `collect_grid()` returns the tidy frames with the
grid columns attached.

**6.4 Plot helpers and the plotting vignette (one PR).** `R/plot.R`, ggplot2
in `Suggests`, `requireNamespace()` guard once in `check_ggplot2()`.

- `plot_prevalence(sim | tidy df)`, `plot_loggers(sim, vars)`, `plot_allele_freqs()`,
  `plot_lineages()`, `plot_treatment()`, `plot_coi(by = c("all", "age_bin", "clinical"))`,
  each returning a ggplot built from its `tidy_*` frame with `period`/`day` on x.
- `plot_replicates(df, y, group = "run", colour = NULL, facet = NULL)` — the
  one plot every project draws: lines per replicate at low alpha plus the mean
  line, faceted by scenario columns from the grid. `gen_plot()`, `cases_plot()`,
  paper-1 Fig 3 and `faking_it.R` are all this function.
- Vignette `vignettes/plotting.Rmd`: one `pipeline()` run per saver, each
  `tidy_*` printed, each `plot_*` shown, then a nine-run `sim_grid()` and
  `plot_replicates()` faceted by EIR, then the selection-coefficient fit (arms).
  Delete `plotting.R` in this PR: `summary_data_frames_from_sims`, the `plot_*`
  and `gg_color_hue` have replacements; `heatmap_ordered_binary_plot()` and
  `sample_coi()` (the age-matched sampler behind paper 1's Fig 1a) are paper-one
  only and move, unchanged, to helper scripts under
  `submission_data/genomic_tools_data/` next to `bootstrap()` and `kl()`, out of
  `R/` (decision 5). Remove `reshape2`, `cowplot`, `RColorBrewer` from
  `Suggests` and the lintr exclusions from `.lintr`.

Tests: `skip_if_not_installed("ggplot2")`; each helper returns a `"ggplot"` with
the expected layers and mapping; no snapshot images.

**6.5 Seasonality and site selection (one PR).** `R/seasonal_utils.R`.

- Export `seasonal_profile(admin, country)` (already there) and add
  `seasonal_profiles(country = NULL)` returning a data frame (country, admin1,
  day, value, plus the `ft`, `a0..b3`, `theta_c` columns from
  `admin_units_seasonal`), and `plot_seasonal_profiles(country, admins = NULL)`
  (`facet_wrap(~ admin1)`, shared y). `admin_match()` stays internal; a
  `find_admin(admin, country)` wrapper that returns the matched row as a data
  frame makes the fuzzy match visible before a 20-year run depends on it.
- Vignette `vignettes/seasonality.Rmd`, "From site to simulation": list the
  countries and admin units in `admin_units_seasonal`, plot a country's
  profiles, pick an admin, show how `equilibrium_init_create()` uses the curve
  and how `use_historic_interventions = TRUE` pulls `itn_2000_2015`/`irs_2000_2015`
  and `ft` for it, `PfPR_to_EIR_heuristic()` from a MAP prevalence (the Uganda
  pattern in `didewin.R`), then a short run with `save_loggers()` and
  `plot_prevalence()` showing the seasonal peak. Delete `country_seasonal_profiles()`.
- `_pkgdown.yml`: new reference groups "Batches and grids", "Tidy results",
  "Plotting", "Seasonality"; articles group "Workflows" with the three vignettes.

Tests: `seasonal_profiles()` row counts and that `seasonal_profile()` matches
the per-admin slice; `plot_seasonal_profiles()` guarded as in 6.4; one snapshot
of `find_admin("Tororo", "Uganda")`.

**Order.** 6.1 first (everything else consumes it, and it kills `cluster_utils.R`),
then 6.2, then 6.3 and 6.4 in either order, 6.5 last. 6.1 and 6.2 are a few days
each; 6.3 depends on access to the cluster to check the DIDE driver path once;
6.4 and 6.5 are mostly vignette writing. The vignettes should be built with
`vignettes_src/build_vignettes.R` as the others are, since 6.3 and the nine-run
grid in 6.4 are too slow for `R CMD check`.

**Dependencies after Phase 6.** `Imports` unchanged (dplyr is already there for
`tidy_*`; `tidyr` is not needed, `stats::reshape` or a `lapply`/`rbind` does the
melts). `Suggests` gains `hipercow` and keeps `ggplot2`; loses `reshape2`,
`cowplot`, `RColorBrewer`. `parallel` is base; no `future` (decision 2).

**Decisions (OJ, 2026-11-13).**

1. The Lancet Microbe 2022 analysis is the MMC code; its compendium is
   `compendia/art_resistance_consensus_modelling`. 6.1 is checked against it as
   the third code base (the lineage, mutation-event, treatment-total and
   time-to-threshold collectors above come from its summaries).
2. Local parallelism: base `parallel::mclapply` only, no `future`.
3. Final-slot stripping: a plain `shrink_sim()` helper the worker can call;
   the default stores everything, so other outputs can be kept and dropped
   later.
4. Cluster results stay in the hipercow store, found by bundle name;
   `collect_grid()` returns tidy frames and the grid is saved alongside the
   bundle name. Per-row rds files are optional, not the default.
5. The barcode heatmap and the age-matched sampler are paper-one only and move
   to helper scripts under `submission_data/genomic_tools_data/` with
   `bootstrap()` and `kl()` in the 6.4 PR; there will be no paper-one
   compendium.
6. Grid encoding of the "history years then treatment years" split: a
   row-to-args function. `sim_grid()` holds whatever columns the user likes;
   the user supplies the function that turns one row into the `pipeline()`
   argument list and `run_grid()` calls it on every row, locally and on the
   cluster. The package ships one default for the history-then-policy shape
   (`history_policy_args()`: history years, policy years, `ft`, `EIR`, `seed`
   to padded flags and vectors) so the standard case is a one-liner, and never
   guesses beyond that.
7. hipercow in CI: yes, through the example driver for tests and vignette
   chunks, DIDE-driver chunks `eval = FALSE`, the r-universe repository in
   `setup-r-dependencies`; confirm against hipercow's own CI and vignette
   setup first.
8. ggplot2 in `Suggests`: yes.
The `docs-revamp` model fixes (the asymptomatic discount, the sterilisation
direction and the zero-length strain course at initialisation) are on the
branch since the rebase; their snapshot movement was reviewed and accepted.

**Phase 0: safety net and hygiene (small PRs, no behaviour change).**
- Characterisation snapshot tests (5.5 item 1).
- `styler::style_pkg()`; `.lintr`; lintr in CI. One commit, no other changes.
- Delete `bootstrap.R`, `kl.R`, `cluster_utils.R`, dead `plotting.R` functions, unused
  utils; move anything wanted to `scripts/` or a companion repo. Remove `redux`.
- DESCRIPTION, README, NAMESPACE exports, `codemeta.json`, vignette chunk names,
  `.Rbuildignore`, testthat edition 3, commit the pending roxygen changes.
- Fix the plain bugs that have no behaviour impact on existing tests
  (`model_param_list_create` half-lives, `microscopy_detected`, `pfpr_to_eir`
  `age_range`, `globalVariables`), each with a test.

**Phase 1: extract without changing behaviour.**
- Pull the duplicated blocks out of `pipeline()` into `setup_fresh()`,
  `setup_from_saved()`, `prepare_state_save()`, `collect_output()`. Snapshots must not
  move. `pipeline()` shrinks to ~150 lines.
- Replace `simulation_R`'s shape dispatch with the five `sim_*` bindings; keep
  `simulation_R` as a deprecated alias for one release.

**Phase 2: parameters and schedule.**
- Classed parameter objects with validators and the `.cpp_contract` test.
- `build_schedule()`; delete the three matrix/flag expanders. This is where the
  padding-rule inconsistencies (4.6) get a single documented answer, so expect a NEWS
  entry and possibly snapshot updates for the affected inputs.
- Export the constructors; un-export the Rcpp wrappers.

**Phase 3: savers, results, state.**
- Savers and the `output =` argument; `magenta_sim` class; `sim_save`/`sim_restore`.
- `pipeline()` maps old flags to a saver with deprecation warnings.

**Phase 4: behaviour-changing fixes, one PR each, with NEWS.**
- Seed replay.
- `min(ft[ft > 0])`, `human_vars` consistency, `metapop` error, single-period `i`.

**Phase 5: tests and docs to production level.**
- Unit tests per 5.5; `COI_df_create` tests; coverage gate at 80% R-side.
- `@examples` on every export; pkgdown reference regrouped around params / run /
  savers / results / data; vignettes rebuilt against the new API.

Rough effort: Phase 0 is a couple of days; Phases 1-3 are the real work, perhaps two
to three weeks of focused time; Phases 4-5 a week. The order matters: nothing in
Phases 1-3 should start before the snapshot tests exist.

---

## Appendix A. Function inventory

Lines = source lines of the definition; Cov = line coverage from `covr`
(`NOT_CRAN=true`); Exp = exported. Verdict: **keep** (as is or styled), **refactor**,
**move** (out of the package), **delete**, **test** (keep but needs tests).

| File | Function | Lines | Args | Cov | Exp |
|---|---|---|---|---|---|
| bootstrap.R | `bootstrap` | 73 | 12 | 0% |  |
| bootstrap.R | `cou_function` | 12 | 9 | 0% |  |
| bootstrap.R | `pibd_function` | 12 | 9 | 0% |  |
| bootstrap.R | `unique_function` | 9 | 9 | 0% |  |
| bootstrap.R | `mean_function` | 3 | 9 | 0% |  |
| cluster_utils.R | `pool_res` | 207 | 9 | 0% |  |
| cluster_utils.R | `create_res_list` | 23 | 8 | 0% |  |
| equilibrium_init_create.R | `equilibrium_init_create` | 312 | 8 | 97% | yes |
| equilibrium_state_create.R | `generate_default_model` | 112 | 5 | 0% |  |
| equilibrium_state_create.R | `spl_create` | 86 | 9 | 100% |  |
| equilibrium_state_create.R | `equilibrium_ss_create` | 29 | 1 | 100% |  |
| kl.R | `kl` | 92 | 7 | 0% |  |
| mainFunctions.R | `simulation_R` | 164 | 2 | 86% | yes |
| mainFunctions.R | `param_list_simulation_init_create` | 47 | 9 | 81% | yes |
| mainFunctions.R | `param_list_simulation_update_create` | 45 | 8 | 65% | yes |
| mainFunctions.R | `param_list_simulation_saved_init_create` | 23 | 1 | 82% | yes |
| mainFunctions.R | `param_list_simulation_finalizer_create` | 17 | 1 | 80% | yes |
| mainFunctions.R | `param_list_simulation_get_create` | 16 | 1 | 75% | yes |
| model_param_list_create.R | `model_param_list_create` | 205 | 72 | 100% | yes |
| pfpr_to_eir.R | `pfpr_to_eir_heuristic` | 108 | 7 | 68% |  |
| pipeline_param_lists.R | `drug_list_create` | 79 | 14 | 96% | yes |
| pipeline_param_lists.R | `drug_create` | 65 | 8 | 100% |  |
| pipeline_param_lists.R | `barcode_list_create` | 39 | 8 | 89% |  |
| pipeline_param_lists.R | `drug_create_dhappq` | 26 | 0 | 0% |  |
| pipeline_param_lists.R | `drug_create_al` | 24 | 0 | 0% |  |
| pipeline_param_lists.R | `drug_create_asaq` | 22 | 0 | 100% |  |
| pipeline_param_lists.R | `nmf_list_create` | 17 | 4 | 100% |  |
| pipeline_param_lists.R | `vector_adaptation_list_create` | 17 | 6 | 100% |  |
| pipeline_param_lists.R | `drug_create_default_no_resistance` | 14 | 0 | 100% |  |
| pipeline_param_lists.R | `housekeeping_list_create` | 13 | 4 | 100% |  |
| pipeline_utils.R | `update_saves` | 128 | 19 | 56% |  |
| pipeline_utils.R | `intervention_grab` | 88 | 6 | 65% | yes |
| pipeline_utils.R | `mu_fv_create` | 80 | 8 | 93% |  |
| pipeline_utils.R | `spl_grab` | 66 | 3 | 77% | yes |
| pipeline_utils.R | `drug_list_update` | 36 | 3 | 29% |  |
| pipeline_utils.R | `spl_matrix_check` | 20 | 2 | 55% |  |
| pipeline_utils.R | `plaf_matrix_check` | 19 | 2 | 60% |  |
| pipeline_utils.R | `age_brackets` | 14 | 3 | 83% |  |
| pipeline_utils.R | `lineages` | 13 | 2 | 100% |  |
| pipeline_utils.R | `pop_alf` | 13 | 3 | 86% |  |
| pipeline.R | `pipeline` | 685 | 48 | 85% | yes |
| pipeline.R | `return_sim_memory` | 19 | 1 | 92% |  |
| plotting.R | `summary_data_frames_from_sims` | 211 | 6 | 0% |  |
| plotting.R | `plot_resistance_lineages` | 130 | 3 | 0% |  |
| plotting.R | `convert_barcode_vectors` | 96 | 6 | 0% | yes |
| plotting.R | `heatmap_ordered_binary_plot` | 87 | 5 | 0% |  |
| plotting.R | `plot_mean_summary` | 52 | 7 | 0% |  |
| plotting.R | `plot_x_loggers` | 41 | 5 | 0% |  |
| plotting.R | `coi_age_plot_sample_x` | 33 | 6 | 0% |  |
| plotting.R | `plot_prevalence` | 33 | 3 | 0% |  |
| plotting.R | `sample_coi` | 33 | 8 | 0% |  |
| plotting.R | `microscopy_detected` | 14 | 4 | 0% |  |
| plotting.R | `gg_color_hue` | 4 | 1 | 0% |  |
| plotting.R | `fd` | 3 | 4 | 100% |  |
| plotting.R | `q_fun` | 3 | 5 | 100% |  |
| seasonal_utils.R | `admin_match` | 55 | 3 | 68% |  |
| seasonal_utils.R | `country_seasonal_profiles` | 33 | 1 | 0% |  |
| seasonal_utils.R | `seasonal_profile` | 31 | 2 | 100% |  |
| seasonal_utils.R | `match_clean` | 13 | 3 | 89% |  |
| stat_genetics_utils.R | `update_saved_state_barcode_plaf` | 58 | 3 | 95% | yes |
| stat_genetics_utils.R | `adjust_matrix_to_plaf` | 47 | 2 | 66% |  |
| stat_genetics_utils.R | `flatten_barcode_vectors` | 27 | 1 | 87% |  |
| stat_genetics_utils.R | `summarySE_mean_only_max_mean` | 23 | 5 | 0% |  |
| stat_genetics_utils.R | `rebuild_barcode_vectors` | 22 | 2 | 100% |  |
| stat_genetics_utils.R | `validate_plaf` | 16 | 3 | 30% |  |
| stat_genetics_utils.R | `ztrgeomintp` | 16 | 3 | 70% |  |
| stat_genetics_utils.R | `flat_to_matrix` | 15 | 2 | 80% |  |
| stat_genetics_utils.R | `summarySE` | 15 | 4 | 0% |  |
| stat_genetics_utils.R | `summarySE_mean_only` | 14 | 4 | 0% |  |
| stat_genetics_utils.R | `convert_ibd_barcode` | 13 | 2 | 0% |  |
| stat_genetics_utils.R | `ztrnbinom` | 13 | 3 | 100% |  |
| stat_genetics_utils.R | `validate_plaf_vector` | 12 | 1 | 57% |  |
| stat_genetics_utils.R | `clonality_from_barcode_list` | 11 | 1 | 0% |  |
| stat_genetics_utils.R | `population_ibd_barcodes` | 11 | 2 | 0% |  |
| stat_genetics_utils.R | `population_ibd_barcodes_c` | 11 | 4 | 0% |  |
| stat_genetics_utils.R | `has_nested_element` | 10 | 2 | 83% |  |
| stat_genetics_utils.R | `pibd_from_barcode_list` | 10 | 2 | 0% |  |
| stat_genetics_utils.R | `set_nested_element` | 10 | 3 | 100% |  |
| stat_genetics_utils.R | `cou_from_barcode_list` | 8 | 1 | 0% |  |
| stat_genetics_utils.R | `get_nested_element` | 7 | 2 | 100% |  |
| universe_utils.R | `COI_df_create` | 248 | 9 | 0% |  |
| universe_utils.R | `pop_strains_df` | 53 | 5 | 84% |  |
| utils.R | `progress_logging` | 36 | 6 | 94% |  |
| utils.R | `person_make` | 17 | 2 | 0% |  |
| utils.R | `binary` | 10 | 2 | 0% |  |
| utils.R | `ranges` | 9 | 2 | 100% |  |
| utils.R | `rbind_list_base` | 8 | 1 | 0% |  |
| utils.R | `assert_scalar_character` | 6 | 2 | 0% |  |
| utils.R | `assert_scalar_logical` | 6 | 2 | 0% |  |
| utils.R | `assert_scalar_numeric` | 6 | 2 | 0% |  |
| utils.R | `bitsToInt` | 6 | 2 | 67% |  |
| utils.R | `func_call` | 6 | 2 | 0% |  |
| utils.R | `assert_null_and_func` | 5 | 3 | 0% |  |

Verdicts by file:

- `pipeline.R`: `pipeline` refactor (section 6); `return_sim_memory` fold into `sim_free`.
- `mainFunctions.R`: `simulation_R` replace with `sim_*` bindings; five `param_list_simulation_*` become internal validators.
- `pipeline_param_lists.R`: all keep; add validators; export; `barcode_list_create` defaults fix.
- `model_param_list_create.R`: refactor to `as.list(environment())` + validator; fix half-life bug.
- `equilibrium_init_create.R`, `equilibrium_state_create.R`: keep the maths; split `equilibrium_init_create` into demography / immunity / human states / mosquito+larval / seasonality pieces; `generate_default_model` is unused by the pipeline (odin is driven from `mu_fv_create`) and can go or be tested.
- `pipeline_utils.R`: `age_brackets` keep+test; `spl_grab`, `intervention_grab` keep+test (and fix `options()` bug); `spl_matrix_check`, `plaf_matrix_check` delete after `build_schedule`; `mu_fv_create` keep+test; `drug_list_update` keep+test; `pop_alf`, `lineages` keep+test; `update_saves` replace with savers.
- `universe_utils.R`: `pop_strains_df` keep+test; `COI_df_create` refactor (one summariser applied over a list of groupings) + test.
- `stat_genetics_utils.R`: `summarySE*` keep one, delete two; `ztr*` keep+test; `clonality_from_barcode_list`, `cou_from_barcode_list` keep+test; `convert_ibd_barcode`, `population_ibd_barcodes*`, `pibd_from_barcode_list` move or delete; `update_saved_state_barcode_plaf` and helpers keep (the style template).
- `seasonal_utils.R`: `match_clean`, `admin_match`, `seasonal_profile` keep+test (fix `quiet`); `country_seasonal_profiles` move.
- `pfpr_to_eir.R`: keep, fix `age_range`, replace the hand-rolled step search with `uniroot()`, test.
- `utils.R`: `progress_logging` rename/keep; `rbind_list_base`, `ranges`, `bitsToInt` keep+test; `binary`, `person_make`, `func_call` delete; `assert_*` keep and use everywhere.
- `plotting.R`: `convert_barcode_vectors` keep+test (exported, in vignette); `q_fun`, `fd` move next to the equilibrium code; everything else move out.
- `bootstrap.R`, `kl.R`, `cluster_utils.R`: move out.
- `data.R`, `magenta.R`: keep; fix `globalVariables`; `@docType package` is deprecated in favour of `"_PACKAGE"` alone.


## Appendix B. Workflow background for Phase 6

What the analysis code around magenta does, from the old `scripts/` folder
(2017-2018), the paper-one plotting script in `submission_data/genomic_tools_data/`,
the arms (2024) and peru_hrp23 (2025) repositories and the MMC compendium.

**The code bases.** The old scripts drove the DIDE cluster through `didehpc` +
`context` + `queuer::enqueue_bulk()` over hand-built `data.frame` grids, with
continuation runs as a second grid whose `saved_state_path` column pointed at
the first's files; the rest of the folder was manual update loops, COI-by-age
summaries and ggplot figures. Paper one ran with `summary_saves_only = TRUE`,
pooled with `summary_data_frames_from_sims()`, then `bootstrap()` / `kl()`, then
a plotting-only script. arms builds a Latin hypercube grid, turns each row into
`pipeline()` arguments in a `for` loop (drug lists, a `plaf` matrix, history then
treatment years), submits ten replicate bundles through didehpc with the seed
drawn on the worker, fetches each result from the context db, strips the
population-length vectors from the final slot and writes one rds per row; its
processing script builds a per-period frame (`cases`, `pcr_prev`, `micro_2_10`,
`tf`, `a_1..a_6`) and fits a selection coefficient per locus. peru_hrp23 is the
only hipercow user: a grid with a seed column, `task_create_bulk_call()`, a
two-stage setup-then-continue through `full_save` and
`update_saved_state_barcode_plaf()`, per-period case and lineage frames, and
replicate-lines-plus-mean plots. The MMC compendium reads per-run milestone
summaries (time to 1/10/25% 580Y, time to 10/25% treatment failure, mutation
events per locus, treatment and failure totals) and monthly lineage frames
(`Strains, Time, rep, Genotype`), produced in the MMC_WP2 work package. Six
steps recur in all of them, each written from scratch: build a grid, submit in
batches, pull down, per-period data frame, bind across replicates and grid,
plot replicate lines with a mean.

**The parked helpers.** `pool_res()` reads one didehpc bundle and builds five
tidy frames (allele frequency, mutations, treatment failure, prevalences,
lineages) for a fixed three-strategy design; its builders are generic, its
fetch and labels are not. `create_res_list()` recovers parameters from bundle
names. `summary_data_frames_from_sims()` and `plot_mean_summary()` are weighted
means of the `summary_*` tables over time; `plot_prevalence()`,
`plot_x_loggers()` and `plot_resistance_lineages()` are the prevalence, logger
and five-panel resistance plots; `heatmap_ordered_binary_plot()`, `sample_coi()`
and `coi_age_plot_sample_x()` served paper one's figures; `gg_color_hue()` is the
ggplot palette; `country_seasonal_profiles()` plots every admin unit's seasonal
curve for a country. Phase 6 maps each onto a `tidy_*` or `plot_*` function or
moves it out, as 6.1-6.5 say.

**hipercow today.** `hipercow_init(driver = "dide-windows")`,
`windows_authenticate()`, `hipercow_provision()` from a `pkgdepends.txt` line
such as `github::OJWatson/magenta@v1.4.0`, `hipercow_environment_create()` for
the worker functions, then `task_create_bulk_expr()` or
`task_create_bulk_call()` over a grid (one task per row), `hipercow_bundle_status()`,
`hipercow_bundle_wait()`, `hipercow_bundle_result()`, `hipercow_bundle_retry()`
for failures and `task_log_show()`; resources through `hipercow_resources()`.
The `"example"` driver runs tasks locally in the same layout. What hipercow
does not do, and where magenta helps, is keep the grid next to the bundle and
turn a bundle's results into tidy frames without holding every `magenta_sim`
in memory at once.
