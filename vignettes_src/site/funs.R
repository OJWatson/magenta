## magenta and malariasimulation runs of one site file unit, reduced to yearly
## PfPR2-10 (microscopy) and clinical incidence (all ages and under 5)

## The two Kenya units, cut from the site file (see run_site_comparison.R)
units <- c("Busia urban", "Nyandarua rural")
unit_sites <- "kenya_units.rds"
years <- 2000:2024
burnin <- 20

read_unit <- function(unit) readRDS(unit_sites)[[unit]]

## magenta ----

## Per save: microscopy PfPR2-10 (each person's detection probability, as
## malariasimulation), and the clinical loggers (daily mean counts)
magenta_summary <- function(sim.out, mpl) {
  l <- sim.out$Loggers
  age <- l$Ages
  st <- l$InfectionStates
  q <- magenta:::q_fun(mpl$d1, l$ID, mpl$ID0, mpl$kD,
                       magenta:::fd(age, mpl$fD0, mpl$aD, mpl$gammaD))
  micro <- (st %in% c(1, 4)) + (st == 2) * q
  two_ten <- age >= 2 * 365 & age < 10 * 365
  data.frame(days = l$Log_Counter, n = length(age), n_05 = sum(age < 5 * 365),
             pfpr_2_10 = mean(micro[two_ten]),
             clinical = l$Incidence, clinical_05 = l$Incidence_05)
}

magenta_save <- function(res, i, sim.out, mpl, num_loci) {
  res[[i]] <- magenta_summary(sim.out, mpl)
  res
}

## AL with no resistance, non-malarial fevers off, saves every 73 days
run_magenta <- function(unit, seed, N) {
  t0 <- Sys.time()
  site <- read_unit(unit)
  dl <- magenta::drug_list_create(drugs = list(magenta::drug_create_al()))
  a <- magenta::site_to_pipeline(site, years = years, burnin = burnin, drug_list = dl,
                                 update_length = 73)
  out <- do.call(magenta::pipeline, c(a, list(
    N = N, seed = seed, update_save_func = magenta_save,
    nmf_list = magenta:::nmf_list_create(nmf_flag = FALSE),
    housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE, cluster = TRUE))))
  mpl <- do.call(magenta::model_param_list_create,
                 a[names(a) %in% names(formals(magenta::model_param_list_create))])
  # pipeline() returns the simulation itself as the last element
  out[[length(out)]] <- magenta_summary(out[[length(out)]], mpl)
  iv <- do.call(rbind, out)
  iv$year <- attr(a, "year")[floor((seq_len(nrow(iv)) - 1) * 73 / 365) + 1]
  iv <- iv[iv$year %in% years, ]
  py <- function(n) tapply(n * iv$days, iv$year, sum) / 365
  cases <- function(x) tapply(x * iv$days, iv$year, sum)
  data.frame(model = "magenta", unit = unit, seed = seed, year = years,
             version = as.character(utils::packageVersion("magenta")),
             pfpr_2_10 = as.numeric(tapply(iv$pfpr_2_10, iv$year, mean)),
             clinical = as.numeric(cases(iv$clinical) / py(iv$n)),
             clinical_05 = as.numeric(cases(iv$clinical_05) / py(iv$n_05)),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

## malariasimulation, as foresite runs it ----
## demography "site" uses the site's life table, as foresite does, and
## "exponential" magenta's exponential population (mean age 21 years), so that
## all-ages incidence can be compared like for like

run_malsim <- function(unit, seed, N, demography = c("site", "exponential")) {
  demography <- match.arg(demography)
  t0 <- Sys.time()
  set.seed(seed)
  site <- read_unit(unit)
  itn <- site$interventions$itn
  site$interventions$itn$implementation$itn_input_dist <- site::site_usage_to_model_distribution(
    usage = itn$use$itn_use,
    usage_year = itn$use$year,
    usage_day_of_year = itn$use$usage_day_of_year,
    distribution_year = itn$implementation$year,
    distribution_day_of_year = itn$implementation$distribution_day_of_year,
    distribution_lower = itn$implementation$distribution_lower,
    distribution_upper = itn$implementation$distribution_upper,
    net_loss_function = netz::net_loss_map,
    half_life = itn$retention_half_life
  )
  p <- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    overrides = list(human_population = N),
    parasite = "falciparum",
    start_year = years[1] - burnin,
    end_year = max(years),
    eir = site$eir$eir[site$eir$sp == "pf"],
    prevalence = c(2, 10) * 365
  )
  if (demography == "exponential") {
    p$custom_demography <- FALSE
    p$average_age <- 21 * 365
    p <- malariasimulation::set_equilibrium(p, init_EIR = site$eir$eir[site$eir$sp == "pf"])
  }
  s <- malariasimulation::run_simulation(timesteps = p$timesteps, parameters = p)
  s$year <- years[1] - burnin + floor((s$timestep - 1) / 365)
  s <- s[s$year %in% years, ]
  # site_parameters() renders incidence in 0-5, 5-15 and 15+
  inc <- grep("^n_inc_clinical_", names(s), value = TRUE)
  age <- sub("n_inc_clinical_", "n_age_", inc)
  py <- function(n) tapply(n, s$year, sum) / 365
  data.frame(model = if (demography == "site") "malariasimulation" else
               "malariasimulation (exponential ages)", unit = unit, seed = seed, year = years,
             version = as.character(utils::packageVersion("malariasimulation")),
             pfpr_2_10 = as.numeric(tapply(s$n_detect_lm_730_3649 / s$n_age_730_3649,
                                           s$year, mean)),
             clinical = as.numeric(tapply(rowSums(s[inc]), s$year, sum) / py(rowSums(s[age]))),
             clinical_05 = as.numeric(tapply(s$n_inc_clinical_0_1824, s$year, sum) /
                                        py(s$n_age_0_1824)),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

## The three runs: magenta, malariasimulation with the site's life table, and
## malariasimulation with magenta's exponential ages
models <- c("magenta", "malariasimulation", "malariasimulation_exp")

run_model <- function(model, unit, seed, N) {
  switch(model,
         magenta = run_magenta(unit, seed, N),
         malariasimulation = run_malsim(unit, seed, N),
         malariasimulation_exp = run_malsim(unit, seed, N, "exponential"))
}
