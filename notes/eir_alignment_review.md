# Why magenta and malariasimulation give different EIR -> PfPR / clinical incidence

Review of the gap seen in `scripts/severe_comparison.R` (magenta clinical incidence 22-36% below
malariasimulation (ms) at the same EIR, PfPR2-10 lower at low EIR). This is a review only: no model
code was changed. The fixes that followed are listed at the end. Numbers below come from a scratch copy of magenta with
environment-variable toggles and event counters, and from ms 3.0.0 with functions swapped via
`assignInNamespace`. Line numbers are for magenta at `cc35c34` and ms at `d315a6c`.

## Short answer

Three things explain most of the gap, in this order:

1. **Bug in magenta's initialisation** (`src/main_init.cpp:375`): people who start infected never
   realise their next few infections. This makes year 1 low and drags down 4-year averages.
2. **magenta blocks infection within uB (7.2 d) of an IB boost** (`src/person.cpp:489`). ms has no
   such rule. At EIR 20 this blocks 37% of otherwise eligible infections.
3. **magenta's realised EIR is lower than the input EIR, and lower than ms's**, by about 25-30% at
   EIR 0.55 and 7-10% at EIR 4-20. It starts about 10% low on day 1 and falls further over the first
   few months, so the vector/transmission side is not in equilibrium with the human side. This is
   what is left at low EIR once 1 and 2 are switched off. I have not traced its cause (see open
   items).

The +0.5 on ICA, ICA boosting while in D/T, the PfPR method and the maternal immunity method are all
real differences, but each moves clinical incidence by 3% or less.

## Results

N = 5e4, ft = 0.4, 4 years, averaged over 2 seeds. "eir" is the realised adult-equivalent annual EIR
(magenta: infectious bites / N / years / mean(psi), with mean(psi) = 0.766; ms: mean(`EIR_gamb`) x 365 / N).

| EIR | run | PfPR2-10 | clin. inc. (pppy) | realised eir |
|---|---|---|---|---|
| 0.55 | magenta as is | .036 | .096 | 0.39 |
| 0.55 | + init fix | .036 | .108 | 0.43 |
| 0.55 | + init fix, no IB gate | .043 | .114 | 0.44 |
| 0.55 | + all toggles below | .041 | .115 | 0.44 |
| 0.55 | **ms** | **.057** | **.152** | **0.61** |
| 4 | magenta as is | .201 | .347 | 3.49 |
| 4 | + init fix | .203 | .377 | 3.72 |
| 4 | + init fix, no IB gate | .210 | .403 | 3.79 |
| 4 | + all toggles below | .208 | .413 | 3.76 |
| 4 | **ms** | **.215** | **.432** | **4.16** |
| 20 | magenta as is | .454 | .515 | 18.3 |
| 20 | + init fix | .460 | .540 | 18.9 |
| 20 | + init fix, no IB gate | .458 | .614 | 19.0 |
| 20 | + all toggles below | .449 | .636 | 18.9 |
| 20 | **ms** | **.455** | **.659** | **20.4** |

"All toggles" = init fix + no IB gate + phi from pre-boost ICA + 0.5 + ICA/ID/IVA boosted only
for S/A/U. With these, magenta is within 4-5% of ms on clinical incidence at EIR 4 and 20. At EIR
0.55 clinical incidence per unit of realised EIR is the same in both models (0.115 / 0.44 = 0.26 and
0.152 / 0.61 = 0.25), so the remaining low-EIR gap is the realised EIR.

Single-toggle effects on clinical incidence, relative to "init fix" (seed-averaged):

| toggle | EIR 0.55 | EIR 4 | EIR 20 |
|---|---|---|---|
| no IB gate | +6% | +7% | +14% |
| phi from pre-boost ICA | -0.5% | +2% | +4% |
| pre-boost ICA + 0.5 | -4% | +1% | +2% |
| boost ICA/ID/IVA only in S/A/U | -2% | 0% | +1% |

In ms, removing the +0.5 changes clinical incidence by +1%, +1% and +2% at EIR 0.55, 4 and 20.
Computing phi from post-boost ICA, as magenta does, changes it by -3%, -1.5% and -1.5% relative to ms as released.

## Each candidate

### A. Initialisation bug: initially infected people cannot realise new infections (NEW, largest at 4 years)

`src/main_init.cpp:375` sets `m_number_of_realised_infections` to the person's starting MOI
(`floor(biting rate) + 1`) for everyone who starts in D/A/U/T. Their realisation vectors start
empty, or with one pending infection from `main_init.cpp:466-481`. Every later check uses
`m_infection_time_realisation_vector.size() > m_number_of_realised_infections`
(`person.cpp:319`, `:1398`, `:1626`), so the person's next MOI scheduled infections, and the
pending one from initialisation, are never realised or logged. This lasts until
`all_strain_clearance()` resets the counter (`person.cpp:778`), i.e. until they recover U->S,
are treated, or die. While this lasts they still get ICA/ID/IVA boosts from the bites.

- Size: at EIR 4 in year 1, 45.6k person-days had an infection scheduled but only 33.2k were
  realised (27% lost; with the fix, 1.3% lost over 4 years). Year-1 clinical incidence is 0.276 as is
  vs 0.356 with the fix. It settles by year 3-4. Over the 4-year window it lowers clinical incidence
  by 8% at EIR 4, 5% at EIR 20 and 11% at EIR 0.55.
- Direction: lowers incidence, prevalence and (through onward infectivity) realised EIR.
- Fix: set `m_number_of_realised_infections` to 0 at initialisation. The initial strains are in
  `m_active_strains`, not in the realisation vectors. Check what else reads that counter on
  reload (`main_saved_init.cpp:696`).

Related, minor: the T branch at `person.cpp:1520-1526` clears the realisation vectors without
resetting `m_number_of_realised_infections`. It looks harmless now, because nobody in T schedules
infections and T->P calls `all_strain_clearance()`, but it is the same pattern.

### B. IB gate on infection (NEW, large at high EIR)

`allocate_infection` only schedules an infection if
`m_IB_last_boost_time < t - uB || (int) m_IB_last_boost_time == t` (`person.cpp:489-490`). So a
successful bite within 7.2 days of an earlier IB-boosting bite (but not on the same day) cannot
infect. ms boosts IB with the same delay but never uses it to block infection
(`R/human_infection.R:30-36`, `:82`). The Griffin equilibrium has no such rule either. The same
condition is commented out in `allocate_bite` (`person.cpp:348-350`).

- Size: blocks 4.5%, 15% and 37% of successful bites on S/A/U people at EIR 0.55, 4 and 20.
  Clinical incidence rises 6%, 7% and 14% when it is removed. PfPR rises at low EIR (+20% at
  0.55) and barely changes at 20.
- Note: magenta boosts ICA/ID/IVA in `allocate_infection` *before* this gate, so blocked bites still
  boost clinical and severe immunity. ms only boosts these on infection.
- Fix: drop the gate (lines 489-490 and the matching brace) to match ms.

### C. The +0.5 on ICA in ms, and pre- vs post-boost ICA

- ms adds 0.5 to non-zero ICA (`R/human_infection.R:1005`; also IVA at `:1020`) and computes phi
  from the ICA **before** today's boost, because `boost_immunity` only queues the update
  (`:988-991`), so `calculate_clinical_infections` reads the old value.
- magenta boosts first (`person.cpp:447`) and then computes phi from the boosted ICA
  (`person.cpp:453`). So when a boost happens magenta uses ICA+1 where ms uses ICA+0.5, and when
  it does not (within uCA) magenta uses ICA where ms uses ICA+0.5.
- magenta's severe probability *is* computed before the IVA boost (`person.cpp:423`), so the two
  immunities are treated differently within magenta.
- Size: small. Up to ±4% on clinical incidence either way (tables above); the direction depends on
  EIR.
- Fix: compute phi before the ICA boost (as already done for theta), then decide separately whether
  to add the ms +0.5.

### D. Boosting while in D/T (candidate 2 in the brief)

Real difference (`person.cpp:381-388` commented out; boosts at `:426-450`; state check only at
`:478-483`), but small. Only 1-2.5% of ICA boosts go to people not in S/A/U, because D+T+P is only
about 1-3% of the population. Restricting boosts to S/A/U changes clinical incidence by 0-2%. The
bigger source of extra boosts is the IB gate (B).

Infections that resolve after the 12-day latent period while the person is in D/T, and so are not
counted (`person.cpp:1408-1419`), are rare: 13 per year at N=5e4, EIR 4. Pending infections wiped
by U->S recovery or death (`person.cpp:873-883`) are about 1% of infections.

### E. Clinical incidence counting

`log_daily_incidence` (`person.cpp:1620-1650`) counts one case per person per day when the next
pending infection, due tomorrow, resolves to D/T and the person is in S/A/U. Once A and B are
removed this matches ms's `calculate_clinical_infections` (`R/human_infection.R:608-640`) in
expectation. ms counts the expected value (`p_inc_clinical`) of one infection per bitten person per
day, from S/A/U only. magenta counts realised draws, which only adds noise. There is no asymptomatic
route to a case in either model. Not a source of the gap.

### F. PfPR2-10 by microscopy

ms (`R/render.R:28-60`) uses per-person q(age, ID) for A, plus D and Tr.
`severe_comparison.R` uses magenta's mean q over all 2-10 year olds, not over the A people. With
per-person q the PfPR is 0.5-1% lower (e.g. .2010 vs .2029 at EIR 4), so the comparison slightly
overstates magenta. Not a source of the gap. The script also compares a single end-of-run snapshot
for magenta against a 4-year average for ms, which only matters while magenta is still drifting
(see G).

### G. Burn-in and drift

- ms: year 1 is 5-10% above later years at EIR 0.55-4 (PfPR .064 -> .053 at 0.55), then stable
  from year 2 to year 8. Mean ICA in adults is flat, and in 2-5 year olds it falls by about 5% over
  8 years.
- magenta as is: clinical incidence rises over years 1-3 (EIR 4: .276, .360, .388, then about
  .38-.39), because of bug A. With A fixed it is flat from year 1 (EIR 4: .356, .384, .392, .392,
  ...), so 4 years is enough once A is fixed.
- magenta's realised EIR with all toggles at EIR 0.55: 0.50 over the first 18 days, 0.45 by 3
  months and 0.44 thereafter. At EIR 4: 3.55, then 3.44, then about 3.76-3.9. ms stays at or a bit
  above the input (0.59-0.66 at 0.55; about 4.1 at 4). See open items.

### H. Other differences checked

- **Maternal immunity**: magenta sets ICM at birth to PM × the mean ICA of 20-21 year olds ×
  the newborn's own biting rate (`person.cpp:857`; mothers summed at `:1546` and
  `main_init.cpp:363-368`). ms samples one age-20 mother in the same heterogeneity group
  (`R/mortality_processes.R:36-70`). Means are similar. It only affects infants (dCM = 68 d), and
  mean ICA in under-ones matches ms closely (e.g. 0.30 in ms at EIR 4). Minor.
- **Latency**: ms has no individual latent period. It uses EIR lagged by de = 12 d
  (`R/biting_process.R:145`), and the outcome and eligibility are decided at the bite. magenta draws
  the outcome and phi at the bite and realises it 12 days later. The only effect is the small losses
  in D, and is minor.
- **Prophylaxis**: magenta T (5 d) -> P (25 d, exponential), with full protection. ms Tr (5 d) ->
  S, with AL prophylaxis decaying as a Weibull (scale about 10.6 d). Treated people in magenta are
  protected about 15 days longer. Under 1% of person-time at EIR 20, so minor.
- **Parameters**: phi0/phi1/IC0/kC/uCA/dCA, b0/b1/IB0/kB/uB/dB, ID/q, PM/dCM, dur_D = dur_T = 5,
  dur_A = 195, dur_U = 110.3 vs 110, average age 21 y vs 7663 d, a0, rho and sigma² all match
  (`src/parameters.cpp:46-105`, ms `get_parameters()`). magenta truncates ages at 100 y.
  ID decays with `g_dB` rather than `g_dID` (`person.cpp:413`, `:1658`); both are 3650 d, so this
  has no effect now, but it would have one if either were changed.

## Ranked list

| rank | difference | real? | direction for magenta | size on clinical inc. (4 y) |
|---|---|---|---|---|
| 1 | Realised EIR below input / below ms (vector-human feedback) | yes, cause not traced | lower PfPR and incidence | ~-25% at EIR 0.55, ~-8% at 4-20 |
| 2 | IB gate blocks infection within uB (`person.cpp:489`) | yes, not in ms or Griffin | lower | -6% / -7% / -14% |
| 3 | Init bug, realised-infection counter (`main_init.cpp:375`) | bug | lower, mostly in years 1-2 | -11% / -8% / -5% |
| 4 | phi from post-boost ICA; no +0.5 | yes | small, sign depends on EIR | within ±4% |
| 5 | Boosts in D/T, and boosts from gate-blocked bites | yes | slightly lower | 0-2% (more is via 2) |
| 6 | PfPR mean-q in the script; snapshot vs average | method | overstates magenta PfPR | ~1% |
| 7 | Maternal immunity method, P vs drug prophylaxis, latency handling | yes | small | <1-2% |

## Suggested fixes (for OJ to decide)

1. Fix the init bug: `m_number_of_realised_infections = 0` at `main_init.cpp:375`. It is a clear
   bug and changes nothing once the transient has passed.
2. Remove the IB gate at `person.cpp:489-490`, so that every successful bite on an S/A/U person can
   infect, as in ms and the Griffin model. This changes magenta's behaviour at high EIR, and
   probably its fitted EIR-PfPR curves, so it is your call.
3. Compute phi before the ICA boost, as theta already is. Then decide whether to adopt ms's +0.5
   (for ICA and, separately, IVA; see the iva_plus_half runs in `severe_comparison.R`). These are
   small either way.
4. Optionally move the ICA/ID/IVA boosts after the gate and state check, so that only infections
   boost them, as in ms.
5. Investigate why magenta's realised EIR starts about 10% below the input and falls further. Do
   this after 1-2 are in, since lower human infectivity from A and B feeds into it. Suggested
   checks: the starting number of infectious mosquitoes against the equilibrium Iv; human
   infectivity (`reciprocal_infection_boolean`, the cA(q) / cU / gamma1 terms, and the
   delay_gam = 12.5 d lag). (`num_bites *= g_Q0` at `main_update.cpp:353` truncates an unsigned int, but that is under one bite a day, so negligible.)
   At low EIR this is the main remaining gap.
6. In `scripts/severe_comparison.R`, use per-person q for magenta PfPR. Ideally compare time
   averages, or the same window, in both models.

## Update: fixes 1-4 and 6 are in, and where the EIR shortfall comes from

Fixes 1-4 and 6 are in commit `7a920bb` (with an IVA left in the paper form, no +0.5). All tests
pass. `scripts/severe_comparison_results.csv` was regenerated (N = 1e5, 4 years):

| EIR | clinical, magenta before | clinical, magenta now | clinical, ms | PfPR2-10, magenta now | PfPR2-10, ms |
|---|---|---|---|---|---|
| 0.55 | 0.096 | 0.114 | 0.151 | 0.044 | 0.056 |
| 4 | 0.346 | 0.413 | 0.441 | 0.204 | 0.217 |
| 20 | 0.519 | 0.638 | 0.672 | 0.451 | 0.459 |

What is left is mostly the realised EIR (item 5). I found two causes in the transmission loop.
Neither is a plain bug: both come from magenta being individual-based in a way the equilibrium
(and ms) is not. So I have not changed either.

**5a. Gametocyte delay removes clinical episodes from onward transmission.**
`reciprocal_infection_boolean` (`src/person.cpp:54-68`) only counts strains acquired more than
`g_delay_gam` (12.5 d) ago.
- The acquisition day is not the infectious bite. It is the day the strain is realised, `dur_E`
  (12 d) after the bite, which is also the day the person moves to D, T or A. The strain is pushed
  with `parameters.g_current_time` at realisation (`src/person.cpp:1477-1480`).
- So a new infection only becomes infectious to mosquitoes 12.5 d after the person enters D or T
  (24.5 d after the bite). D and T each last about 5 d (a clinical case goes to one or the other),
  so the clinical episode itself contributes almost nothing.
  - An untreated case is usually A by day 12.5, so it transmits at cA, not cD.
  - A treated case is in P (cP = 0) by then, so cT never counts.
- ms lags the whole person's infectivity by 12.5 d. A person in D on day t contributes cD on day
  t + 12.5, whatever state they are in then. So the D and T days are moved later, not lost.
- The equilibrium solution also counts them in full (`FOIvij_eq` uses cD * D + cT * T,
  `R/equilibrium_init_create.R:247`). The lag makes no difference at equilibrium, which is why
  setting `g_delay_gam = 0` is a fair steady-state probe of the ms behaviour.

To match ms, keep the 12.5 d lag but apply it to the person's infectivity history: the
infectivity to use on day t is the person's state (D, T, A or U) on day t - 12.5. Strain age
should then only pick which strains a mosquito gets.

**5b. Mosquitoes bite every 3 days, but sporozoites appear at day 10.**
- A mosquito can only be infected at a bite. It then bites on a fixed 3-day cycle
  (`g_mosquito_next_biting_day_vector` = 3, `src/parameters.cpp:68`).
- The oocyst ruptures 10 d after the infecting bite (`src/mosquito.cpp:67`), so the first
  infectious bite is always on day 12. It is always 2 days late, never 0 or 1.
- Death is floor(Exp(mean 1/mu0)) + 1 days (`src/mosquito.cpp:254-263`), so P(alive on day k) =
  exp(-mu0 k).
- Expected infectious bites per mosquito infection:
  - magenta: exp(-12 mu0) / (1 - exp(-3 mu0)) = 0.627;
  - the ODE and equilibrium: exp(-10 mu0) / (3 mu0) = 0.675.
- The ODE number is what you get in the discrete model if the first infectious bite is on day
  11.45. In effect magenta is about half a day late, not 2 days. That half day still costs 7% of
  the infectious bites on each pass of the loop.
- The realised EIR starts at the input on day 0, because infectious mosquitoes are set directly
  from Iv. It then drifts down, most at low EIR.

**Probes** (scratch toggles, N = 5e4, 4 years, 2 seeds, all of fixes 1-4 on):
- `eip`: 15.7% of infecting bites get a rupture 3 d earlier (first infectious bite on day 9, not
  12). This makes the expected infectious bites match the ODE's 0.675 without changing the
  biting cycle. The draw is made once per bite: drawing per oocyst left rupture times out of order
  and stalled mosquitoes.
- `gam0`: `g_delay_gam = 0`.

The ms column is the mean of years 1-4 from 8-year runs (N = 5e4), and clinical incidence for ms is
from the N = 1e5 comparison.

| EIR | variant | adult-equiv realised EIR | ratio to input | PfPR2-10 | clinical |
|---|---|---|---|---|---|
| 0.55 | now | 0.436 | 0.79 | 0.041 | 0.115 |
| 0.55 | eip | 0.507 | 0.92 | 0.050 | 0.132 |
| 0.55 | gam0 | 0.481 | 0.87 | 0.047 | 0.126 |
| 0.55 | eip + gam0 | 0.564 | 1.03 | 0.055 | 0.146 |
| 0.55 | ms | 0.607 | 1.10 | 0.057 | 0.151 |
| 4 | now | 3.76 | 0.94 | 0.208 | 0.413 |
| 4 | eip | 4.13 | 1.03 | 0.220 | 0.446 |
| 4 | gam0 | 3.93 | 0.98 | 0.213 | 0.425 |
| 4 | eip + gam0 | 4.29 | 1.07 | 0.220 | 0.457 |
| 4 | ms | 4.16 | 1.04 | 0.217 | 0.441 |
| 20 | now | 18.9 | 0.94 | 0.449 | 0.635 |
| 20 | eip | 20.4 | 1.02 | 0.464 | 0.669 |
| 20 | gam0 | 19.2 | 0.96 | 0.450 | 0.644 |
| 20 | eip + gam0 | 20.8 | 1.04 | 0.470 | 0.677 |
| 20 | ms | 20.4 | 1.02 | 0.459 | 0.672 |

The biting-cycle timing (5b) is the main cause. Correcting it alone brings EIR 4 and 20 into line
with ms on realised EIR, PfPR and clinical incidence. EIR 0.55 is still about 15% below ms on EIR
and clinical incidence. Adding the gametocyte change (5a) closes most of that, but takes EIR 4 and
20 slightly past ms.

ms is not flat either. At an input of 0.55 its realised EIR is 0.66 in year 1, then 0.60, 0.58
and 0.59, settling about 6% above the input (0.585 by year 8). At 4 and 20 it settles 1-2% above.

Ways to fix 5b, all of them behaviour changes:
1. Draw the extrinsic incubation period per bite as in the probe, so the expected infectious bites
   match the ODE. This is a small change in `allocate_gametocytes`, but p depends on mu0 and the
   cycle length, so it would need recomputing for interventions that change either.
2. Recompute the equilibrium mosquito density (`mv0`) for the discrete 3-day cycle, so the input
   EIR is what is realised. This keeps the mosquito biology as it is.
3. Make the gonotrophic cycle stochastic (for example geometric with mean 3 d), which matches the
   ODE directly. This is the biggest change to the mosquito code and to interventions that use
   the biting vector.

**5c. Minor, initialisation of exposed mosquitoes** (`src/main_init.cpp:600-622`):
- Their remaining EIP is drawn from an exponential truncated at 10 d, not from the distribution
  the equilibrium implies.
- The multi-strain branch sets the rupture time without `g_current_time`. That is harmless at
  t = 0 but wrong otherwise.
Both only affect the first weeks.

Also note ms itself sits above the input EIR (0.61 at an input of 0.55), so ms is not exactly at
its own equilibrium either. Some of the remaining gap at EIR 0.55 is ms being high, not magenta
being low.

## How the numbers were made

Scratch copy of this branch, with environment toggles `MAG_FIX_INIT` (counter = 0 at init),
`MAG_NO_IBGATE`, `MAG_PREBOOST_PHI`, `MAG_PLUS_HALF` and `MAG_SAU_BOOST`, plus counters on scheduled
and realised infections, ICA boosts by state, and pending infections wiped by recovery or treatment.
These were run with `pipeline(EIR, N = 5e4, years, ft = 0.4, num_loci = 2)` for 4 years (seeds 11
and 12) and for 1, 2, 3, 5, 6 and 8 years with the same seed, to get per-year values by
differencing. ms 3.0.0 was run for 8 years at N = 5e4 (seeds 11 and 12) as released, without the
+0.5 on ICA, and with phi from post-boost ICA, using the same drug and treatment set-up as
`severe_comparison.R`. None of this code is on the branch.
