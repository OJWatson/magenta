# Severe disease, hospitalisations and deaths: magenta vs malariasimulation
#
# Runs magenta and malariasimulation (3.0.0, with postie 1.1.0 for the
# hospitalisation and death assumptions) from equilibrium in matched perennial
# settings (same EIR, same ft, no drug failure) at EIRs that give a PfPR2-10
# (microscopy) of about 5%, 20% and 45% in magenta's equilibrium solution.
#
# malariasimulation is run twice per EIR: as released, where severe_immunity()
# adds 0.5 to any non-zero IVA, and with that removed, which is the form in
# Griffin et al. 2015 Proc R Soc B and the one magenta uses.
#
# Severe incidence is the expected value in both models (magenta's Severe
# loggers, malariasimulation's p_inc_severe renderer), so it carries less noise
# than counting draws.
#
# Usage (from the package root, with this branch of magenta installed):
#   Rscript scripts/severe_comparison.R
# N, YEARS and CORES can be set as environment variables. Writes
# scripts/severe_comparison_results.csv

library(magenta)
library(malariasimulation)

N <- as.numeric(Sys.getenv("N", 1e5))
YEARS <- as.numeric(Sys.getenv("YEARS", 4))
CORES <- as.numeric(Sys.getenv("CORES", 9))
FT <- 0.4
EIRS <- c(0.55, 4, 20)
out_csv <- file.path("scripts", "severe_comparison_results.csv")

# ages in days, magenta's under 5 logger is age <= 5 * 365
age_max <- 100 * 365
age_05 <- 5 * 365

# ---------------------------------------------------------------------------
# magenta
# ---------------------------------------------------------------------------
run_magenta <- function(eir, seed) {

  out <- pipeline(EIR = eir, N = N, years = YEARS, ft = FT, num_loci = 2,
                  seed = seed,
                  housekeeping_list = magenta:::housekeeping_list_create(
                    quiet = TRUE, quiet_test = TRUE, cluster = TRUE))
  lg <- out$Loggers
  rates <- severe_rates(lg)

  # PfPR2-10 by microscopy from the final population, as in pipeline_utils
  mpl <- model_param_list_create()
  two_ten <- lg$Ages < 3650 & lg$Ages > 720
  q <- magenta:::q_fun(mpl$d1, lg$ID, mpl$ID0, mpl$kD,
                       sapply(lg$Ages, magenta:::fd, mpl$fD0, mpl$aD, mpl$gammaD))
  pfpr <- sum(lg$InfectionStates %in% c(1, 4) & two_ten) / sum(two_ten) +
    sum(lg$InfectionStates == 2 & two_ten) / sum(two_ten) * mean(q[two_ten])

  data.frame(model = "magenta", iva_plus_half = FALSE, EIR = eir,
             pfpr_2_10 = pfpr,
             clinical_all = rates$clinical[1], clinical_05 = rates$clinical[2],
             severe_all = rates$severe[1], severe_05 = rates$severe[2])
}

# ---------------------------------------------------------------------------
# malariasimulation
# ---------------------------------------------------------------------------
severe_immunity_paper <- function(age, acquired_immunity, maternal_immunity,
                                  parameters) {
  fv <- 1 - (1 - parameters$fv0) / (1 + (age / parameters$av)^parameters$gammav)
  parameters$theta0 * (parameters$theta1 + (1 - parameters$theta1) /
    (1 + fv * ((acquired_immunity + maternal_immunity) / parameters$iv0)^parameters$kv))
}

run_malsim <- function(eir, seed, plus_half) {

  # runs in its own forked process, so this only changes this run
  if (!plus_half) {
    utils::assignInNamespace("severe_immunity", severe_immunity_paper,
                             "malariasimulation")
  }
  set.seed(seed)

  p <- get_parameters(list(
    human_population = N,
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365 - 1,
    clinical_incidence_rendering_min_ages = c(0, 0),
    clinical_incidence_rendering_max_ages = c(age_05, age_max),
    severe_incidence_rendering_min_ages = c(0, 0),
    severe_incidence_rendering_max_ages = c(age_05, age_max)
  ))
  # no drug failure, so ft treated cases all clear as in magenta
  al <- AL_params
  al[1] <- 1
  p <- set_drugs(p, list(al))
  p <- set_clinical_treatment(p, 1, 1, FT)
  p <- set_equilibrium(p, eir)

  o <- run_simulation(YEARS * 365, p)
  s <- function(x) sum(o[[x]])
  a05 <- paste0("_0_", age_05)
  aall <- paste0("_0_", age_max)

  data.frame(model = "malariasimulation", iva_plus_half = plus_half, EIR = eir,
             pfpr_2_10 = s("n_detect_lm_730_3649") / s("n_age_730_3649"),
             clinical_all = 365 * s(paste0("p_inc_clinical", aall)) / s(paste0("n_age", aall)),
             clinical_05 = 365 * s(paste0("p_inc_clinical", a05)) / s(paste0("n_age", a05)),
             severe_all = 365 * s(paste0("p_inc_severe", aall)) / s(paste0("n_age", aall)),
             severe_05 = 365 * s(paste0("p_inc_severe", a05)) / s(paste0("n_age", a05)))
}

# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
jobs <- expand.grid(model = c("magenta", "malsim_paper", "malsim_plus_half"),
                    eir = EIRS, stringsAsFactors = FALSE)

res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]
  switch(j$model,
         magenta = run_magenta(j$eir, seed = 100L + i),
         malsim_paper = run_malsim(j$eir, seed = 100L + i, plus_half = FALSE),
         malsim_plus_half = run_malsim(j$eir, seed = 100L + i, plus_half = TRUE))
}, mc.cores = CORES, mc.preschedule = FALSE)

failed <- vapply(res, inherits, logical(1), "try-error")
if (any(failed)) {
  print(res[failed])
  stop("some runs failed")
}
res <- do.call(rbind, res)

# hospitalisations and deaths per 1000 per year with the postie assumptions
for (ag in c("all", "05")) {
  so <- severe_outcomes(res[[paste0("severe_", ag)]], ft = FT)
  res[[paste0("hosp_per_1000_", ag)]] <- 1000 * so$severe_hospital
  res[[paste0("deaths_per_1000_", ag)]] <- 1000 * so$mortality
}

# check against postie itself for the malariasimulation rows
ms <- res$model == "malariasimulation"
pst <- postie:::severe_incidence_mortality(
  data.frame(severe = res$severe_all[ms], ft = FT, ft_sev = 0.8),
  treatment_scaler = 0.42, hosp_sev_cfr = 0.065, community_sev_cfr = 0.6)
stopifnot(isTRUE(all.equal(1000 * pst$mortality, res$deaths_per_1000_all[ms])))

res$N <- N
res$years <- YEARS
res$ft <- FT
res[] <- lapply(res, function(x) if (is.numeric(x)) signif(x, 4) else x)
res <- res[order(res$EIR, res$model, res$iva_plus_half), ]
write.csv(res, out_csv, row.names = FALSE)
print(res)
