suppressPackageStartupMessages({
  library(magenta)
  library(malariasimulation)
})

FT <- 0.4
CHUNK <- 73L          # days per output period (5 per year)
age_max <- 100 * 365
age_05 <- 5 * 365

# pipeline(update_save = TRUE) calls simulation_R(pl, seed = seed) for every
# chunk, which re-runs set.seed(seed) with the same seed each chunk. Wrap it so
# each call within a run gets a different seed (counter is per forked process,
# reset in run_magenta).
.seed_env <- new.env()
.seed_env$n <- 0L
.orig_sim_R <- magenta:::simulation_R
.wrapped_sim_R <- function(param_list, seed) {
  .seed_env$n <- .seed_env$n + 1L
  .orig_sim_R(param_list, seed = seed + 1000003L * .seed_env$n)
}
utils::assignInNamespace("simulation_R", .wrapped_sim_R, "magenta")

mpl_default <- model_param_list_create()

summ_loggers <- function(lg, mpl = mpl_default) {
  ages <- lg$Ages
  st <- lg$InfectionStates
  n <- length(ages)
  two_ten <- ages < 3650 & ages > 720
  q <- magenta:::q_fun(mpl$d1, lg$ID, mpl$ID0, mpl$kD,
                       sapply(ages, magenta:::fd, mpl$fD0, mpl$aD, mpl$gammaD))
  pfpr <- (sum(st %in% c(1, 4) & two_ten) + sum(q[st == 2 & two_ten])) / sum(two_ten)
  days <- lg$Log_Counter
  data.frame(days = days,
             realised_eir = 365 * lg$Treatments$daily_infectious_bite_counters / n / days /
               mean(1 - mpl$rho * exp(-ages / mpl$a0)),
             pfpr_2_10 = pfpr,
             clinical_all = 365 * lg$Incidence / n,
             clinical_05 = 365 * lg$Incidence_05 / sum(ages <= age_05))
}

run_magenta <- function(eir, seed, N, years) {
  .seed_env$n <- 0L
  save_fn <- function(res, i, sim.out, mpl, num_loci) {
    res[[i]] <- summ_loggers(sim.out$Loggers)
    res
  }
  out <- pipeline(EIR = eir, N = N, years = years, ft = FT, num_loci = 2,
                  seed = seed, update_save = TRUE, update_length = CHUNK,
                  update_save_func = save_fn,
                  housekeeping_list = magenta:::housekeeping_list_create(
                    quiet = TRUE, quiet_test = TRUE, cluster = TRUE))
  k <- length(out)
  out[[k]] <- summ_loggers(out[[k]]$Loggers)
  df <- do.call(rbind, out[seq_len(k)])
  df$period <- seq_len(nrow(df))
  df$t_end <- cumsum(df$days)
  df$model <- "magenta"; df$EIR <- eir; df$seed <- seed
  df
}

run_malsim <- function(eir, seed, N, years) {
  set.seed(seed)
  p <- get_parameters(list(
    human_population = N,
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365 - 1,
    clinical_incidence_rendering_min_ages = c(0, 0),
    clinical_incidence_rendering_max_ages = c(age_05, age_max)
  ))
  al <- AL_params
  al[1] <- 1
  p <- set_drugs(p, list(al))
  p <- set_clinical_treatment(p, 1, 1, FT)
  p <- set_equilibrium(p, eir)
  o <- run_simulation(years * 365, p)
  per <- (o$timestep - 1) %/% CHUNK + 1
  a05 <- paste0("_0_", age_05); aall <- paste0("_0_", age_max)
  ag <- function(x) tapply(x, per, sum)
  df <- data.frame(
    days = as.numeric(ag(rep(1, nrow(o)))),
    realised_eir = 365 * as.numeric(ag(o$EIR_gamb)) / N / as.numeric(ag(rep(1, nrow(o)))),
    pfpr_2_10 = as.numeric(ag(o$n_detect_lm_730_3649 / o$n_age_730_3649)) /
      as.numeric(ag(rep(1, nrow(o)))),
    clinical_all = 365 * as.numeric(ag(o[[paste0("p_inc_clinical", aall)]])) /
      as.numeric(ag(o[[paste0("n_age", aall)]])),
    clinical_05 = 365 * as.numeric(ag(o[[paste0("p_inc_clinical", a05)]])) /
      as.numeric(ag(o[[paste0("n_age", a05)]])))
  df$period <- seq_len(nrow(df))
  df$t_end <- cumsum(df$days)
  df$model <- "malariasimulation"; df$EIR <- eir; df$seed <- seed
  df
}

equilibrium_outputs <- function(eir) {
  eq <- equilibrium_init_create(
    age_vector = c(seq(0, 2, 0.25), 3:10, seq(12.5, 20, 2.5), seq(25, 100, 5)),
    het_brackets = 5, EIR = eir, ft = FT,
    model_param_list = model_param_list_create(),
    country = NULL, admin = NULL, quiet = TRUE)
  by_age <- function(x) rowSums(x[, , 1])
  two_ten <- eq$age >= 730 - 1 & eq$age < 3650 - 1
  u5 <- eq$age < age_05 - 1
  detected <- by_age(eq$init_D) + by_age(eq$init_T) + rowSums(eq$init_A[, , 1] * eq$p_det_eq)
  clin <- rowSums(eq$phi_eq * eq$FOI * (eq$init_S + eq$init_A + eq$init_U)[, , 1])
  data.frame(EIR = eir, realised_eir = eir,
             pfpr_2_10 = sum(detected[two_ten]) / sum(eq$den[two_ten]),
             clinical_all = 365 * sum(clin) / sum(eq$den),
             clinical_05 = 365 * sum(clin[u5]) / sum(eq$den[u5]))
}
