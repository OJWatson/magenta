# Forensic output collection for magenta vs malariasimulation vs deterministic.
# Sourced after ../funs.R (which loads both packages and wraps magenta's
# simulation_R so each chunk gets a different seed).

BANDS_LO <- c(0, 1, 2, 5, 10, 15, 20, 40) * 365
BAND_NAMES <- c("0-1", "1-2", "2-5", "5-10", "10-15", "15-20", "20-40", "40+")
band_of <- function(age_days) {
  factor(BAND_NAMES[findInterval(age_days, BANDS_LO)], levels = BAND_NAMES)
}
# discrete-day mean sojourn of rexpint1(d) + 1, i.e. floor(Exp(mean d)) + 1
disc_mean <- function(d) 1 / (1 - exp(-1 / d))

# Per-band means of a set of per-person columns, plus n
band_means <- function(band, df) {
  out <- lapply(df, function(x) as.numeric(tapply(x, band, mean)))
  out <- as.data.frame(out)
  out$n <- as.numeric(table(band))
  out$band <- levels(band)
  out
}

# ---------------------------------------------------------------- magenta ---

# Summarise a magenta snapshot. lg = Loggers (or NULL at t=0), st = full state
# from the getter (or NULL), period-level counters come from lg.
magenta_snapshot <- function(lg, st, mpl, t_now) {
  if (!is.null(lg)) {
    ages <- lg$Ages; is <- lg$InfectionStates
    imm <- data.frame(IB = lg$IB, ICA = lg$ICA, ICM = lg$ICM, ID = lg$ID,
                      IVA = lg$IVA, IVM = lg$IVM)
  } else {
    pl <- st$population_List
    ages <- pl$Ages; is <- pl$Infection_States
    imm <- data.frame(IB = pl$IB, ICA = pl$ICA, ICM = pl$ICM, ID = pl$ID,
                      IVA = pl$IVA, IVM = pl$IVM)
  }
  N <- length(ages)
  band <- band_of(ages)
  fd <- magenta:::fd(ages, mpl$fD0, mpl$aD, mpl$gammaD)
  q <- mpl$d1 + (1 - mpl$d1) / (1 + fd * (imm$ID / mpl$ID0)^mpl$kD)
  psi <- 1 - mpl$rho * exp(-ages / mpl$a0)
  per <- data.frame(S = is == 0, D = is == 1, A = is == 2, U = is == 3,
                    T = is == 4, P = is == 5)
  per$prev_micro <- per$D + per$T + per$A * q
  per$prev_any <- per$D + per$A + per$U + per$T
  # clinical incidence per person-day by Little's law with magenta's discrete
  # sojourns (validated against the logged totals)
  per$clin_little <- per$D / disc_mean(1 / mpl$rD) + per$T / disc_mean(1 / mpl$rT)
  per <- cbind(per, imm)
  # probability scales from immunity (person level, then averaged by band)
  per$b_inf <- mpl$b0 * (mpl$b1 + (1 - mpl$b1) / (1 + (imm$IB / mpl$IB0)^mpl$kB))
  per$phi_clin <- mpl$phi0 * (mpl$phi1 + (1 - mpl$phi1) /
                    (1 + ((imm$ICA + imm$ICM) / mpl$IC0)^mpl$kC))
  fv_ <- 1 - (1 - mpl$fV0) / (1 + (ages / mpl$aV)^mpl$gammaV)
  per$theta_sev <- mpl$theta0 * (mpl$theta1 + (1 - mpl$theta1) /
                    (1 + fv_ * ((ifelse(imm$IVA > 0, imm$IVA + 0.5, imm$IVA) + imm$IVM) / mpl$IV0)^mpl$kV))
  per$q_det <- q
  out <- list(N = N, omega = mean(psi), t = t_now)
  if (!is.null(st)) {
    pl <- st$population_List
    ev <- st$populations_event_and_strains_List
    zeta <- pl$Zetas
    pi_ <- zeta * psi
    tt <- t_now - 1  # last simulated day
    lag_s <- floor(mpl$delayGam); lag_l <- ceiling(mpl$delayGam)
    w <- mpl$delayGam - lag_s
    hd <- ev$State_history_days; hs <- ev$State_history_states
    lagged <- function(day) {
      vapply(seq_len(N), function(i) {
        d <- hd[[i]]
        if (length(d) == 0) return(as.integer(is[i]))
        k <- max(1L, sum(d <= day))
        as.integer(hs[[i]][k])
      }, integer(1))
    }
    s12 <- lagged(tt - lag_s); s13 <- lagged(tt - lag_l)
    cA <- mpl$cU + (mpl$cD - mpl$cU) * q^mpl$gamma1
    cfun <- function(s) ifelse(s == 1, mpl$cD, ifelse(s == 2, cA,
                        ifelse(s == 3, mpl$cU, ifelse(s == 4, mpl$cT, 0))))
    acq <- ev$Strain_day_of_acquisition_vectors
    cacq <- ev$Cleared_strain_day_of_acquisition_vectors
    cday <- ev$Cleared_strain_days
    gam <- vapply(seq_len(N), function(i) {
      any(unlist(acq[[i]]) <= tt - lag_s) ||
        any(unlist(cacq[[i]]) <= tt - lag_s & unlist(cday[[i]]) > tt - lag_l)
    }, logical(1))
    p_lag <- gam * ((1 - w) * cfun(s12) + w * cfun(s13))
    p_now <- cfun(is)
    per$pi <- pi_
    per$p_lag <- p_lag
    per$pi_p_lag <- pi_ * p_lag
    per$pi_p_now <- pi_ * p_now
    per$gam <- gam
    out$kappa_lag <- sum(pi_ * p_lag) / sum(pi_)
    out$kappa_now <- sum(pi_ * p_now) / sum(pi_)
    out$mean_zeta <- mean(zeta)
    out$sum_pi <- sum(pi_)
    # mosquitoes
    sc <- st$scourge_List
    on <- !sc$Mosquito_Off_Season
    ms <- sc$Mosquito_Infection_States[on]
    out$mos_n <- sum(on); out$mos_total <- length(on)
    out$mos_S_snap <- mean(ms == 0); out$mos_E_snap <- mean(ms == 1)
    out$mos_I_snap <- mean(ms == 2)
    out$mos_resid_life <- mean(sc$Mosquito_Day_of_death[on] - t_now)
    out$mos_next_meal <- mean(sc$Mosquito_Day_of_next_blood_meal[on] - t_now)
    out$g_mean_mv <- st$parameters_List$g_mean_mv
  }
  out$by_band <- band_means(band, per)
  out
}

magenta_period <- function(lg, mpl) {
  n <- length(lg$Ages); days <- lg$Log_Counter
  tr <- lg$Treatments
  data.frame(days = days,
             bites = tr$daily_bite_counters, inf_bites = tr$daily_infectious_bite_counters,
             inc_all = lg$Incidence / n,                       # per person-day
             inc_05 = lg$Incidence_05 / sum(lg$Ages <= 5 * 365),
             sev_all = lg$Severe / n, sev_05 = lg$Severe_05 / sum(lg$Ages <= 5 * 365),
             S = lg$S, D = lg$D, A = lg$A, U = lg$U, T = lg$T, P = lg$P,
             Mos_S = lg$Mosquitoes$Mos_S, Mos_E = lg$Mosquitoes$Mos_E,
             Mos_I = lg$Mosquitoes$Mos_I)
}

magenta_get_state <- function(ptr, seed) {
  pl2 <- magenta:::param_list_simulation_get_create(statePtr = ptr)
  magenta:::simulation_R(pl2, seed = seed)
}

run_magenta_forensic <- function(eir, seed, N, years, chunk = CHUNK) {
  mpl <- model_param_list_create()
  hk <- magenta:::housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)
  # t = 0: same seed counter as the main run, so the same initial state
  .seed_env$n <- 0L
  st0 <- pipeline(EIR = eir, N = N, years = years, ft = FT, num_loci = 2, seed = seed,
                  set_up_only = TRUE, full_save = TRUE, housekeeping_list = hk)
  snaps <- list(magenta_snapshot(NULL, st0, mpl, 0L))
  rm(st0); gc()
  periods <- list()
  save_fn <- function(res, i, sim.out, mpl, num_loci) {
    st <- magenta_get_state(sim.out$Ptr, seed)
    t_now <- st$parameters_List$g_current_time
    res[[i]] <- list(period = magenta_period(sim.out$Loggers, mpl),
                     snap = magenta_snapshot(sim.out$Loggers, st, mpl, t_now))
    res
  }
  .seed_env$n <- 0L
  out <- pipeline(EIR = eir, N = N, years = years, ft = FT, num_loci = 2,
                  seed = seed, update_save = TRUE, update_length = chunk,
                  update_save_func = save_fn, housekeeping_list = hk)
  k <- length(out)
  lgk <- out[[k]]$Loggers
  out[[k]] <- list(period = magenta_period(lgk, mpl),
                   snap = magenta_snapshot(lgk, NULL, mpl, chunk * k))
  list(model = "magenta", EIR = eir, seed = seed, N = N,
       snaps = c(snaps, lapply(out[seq_len(k)], `[[`, "snap")),
       periods = do.call(rbind, lapply(out[seq_len(k)], `[[`, "period")))
}

# ------------------------------------------------------- malariasimulation ---

malsim_params <- function(eir, N, p25 = FALSE) {
  mn <- BANDS_LO; mx <- c(BANDS_LO[-1] - 1, age_max)
  r <- list(human_population = N)
  # align ms to magenta's parameter values (ms rounds several published values)
  mp <- mpl_default
  r$average_age <- 1 / mp$eta; r$du <- 1 / mp$rU; r$cd <- mp$cD; r$cu <- mp$cU; r$ct <- mp$cT
  r$ad <- mp$aD; r$b0 <- mp$b0; r$ib0 <- mp$IB0; r$kb <- mp$kB; r$ub <- mp$uB
  r$phi0 <- mp$phi0; r$phi1 <- mp$phi1; r$uc <- mp$uCA
  for (v in c("prevalence", "clinical_incidence", "severe_incidence", "age_group",
              "ib", "id", "ica", "iva", "icm", "ivm")) {
    r[[paste0(v, "_rendering_min_ages")]] <- mn
    r[[paste0(v, "_rendering_max_ages")]] <- mx
  }
  p <- get_parameters(r)
  al <- AL_params; al[1] <- 1
  # p25: prophylaxis matched to magenta (exponential, mean 30.5 d from treatment)
  if (p25) { al[3] <- 1; al[4] <- 30.5 }
  p <- set_drugs(p, list(al))
  p <- set_clinical_treatment(p, 1, 1, FT)
  set_equilibrium(p, eir)
}

malsim_snapshot <- function(state, p) {
  v <- state$individual$variables
  t <- state$timesteps
  N <- p$human_population
  is <- integer(N)
  cats <- c(S = 0L, D = 1L, A = 2L, U = 3L, Tr = 4L)
  for (k in names(cats)) is[v$state[[k]]] <- cats[[k]]
  ages <- t - v$birth
  band <- band_of(ages)
  psi <- 1 - p$rho * exp(-ages / p$a0)
  pi_ <- v$zeta * psi
  fd <- 1 - (1 - p$fd0) / (1 + (ages / p$ad)^p$gammad)
  q <- p$d1 + (1 - p$d1) / (1 + fd * (v$id / p$id0)^p$kd)
  per <- data.frame(S = is == 0, D = is == 1, A = is == 2, U = is == 3,
                    T = is == 4, P = FALSE)
  per$prev_micro <- per$D + per$T + per$A * q
  per$prev_any <- per$D + per$A + per$U + per$T
  per$clin_little <- NA_real_
  per$IB <- v$ib; per$ICA <- v$ica; per$ICM <- v$icm; per$ID <- v$id
  per$IVA <- v$iva; per$IVM <- v$ivm
  per$b_inf <- p$b0 * (p$b1 + (1 - p$b1) / (1 + (v$ib / p$ib0)^p$kb))
  per$phi_clin <- p$phi0 * (p$phi1 + (1 - p$phi1) / (1 + ((v$ica + v$icm) / p$ic0)^p$kc))
  fv_ <- 1 - (1 - p$fv0) / (1 + (ages / p$av)^p$gammav)
  per$theta_sev <- p$theta0 * (p$theta1 + (1 - p$theta1) /
                    (1 + fv_ * ((ifelse(v$iva > 0, v$iva + 0.5, v$iva) + v$ivm) / p$iv0)^p$kv))
  per$q_det <- q
  per$pi <- pi_; per$p_lag <- NA_real_; per$pi_p_lag <- NA_real_
  per$pi_p_now <- pi_ * v$infectivity
  per$gam <- NA
  list(N = N, omega = mean(psi), t = t, kappa_now = sum(pi_ * v$infectivity) / sum(pi_),
       mean_zeta = mean(v$zeta), sum_pi = sum(pi_),
       by_band = band_means(band, per))
}

malsim_period <- function(d, p) {
  a <- p$Q0 * p$blood_meal_rates
  inc <- function(pre) {
    cols <- grep(paste0("^", pre, "_"), names(d), value = TRUE)
    sapply(cols, function(cc) sum(d[[cc]]))
  }
  nage <- inc("n_age")
  out <- data.frame(days = nrow(d),
    EIR_adult = 365 * mean(d$EIR_gamb) / p$human_population,
    FOIM = mean(d$FOIM_gamb), kappa_lag = mean(d$FOIM_gamb) / a, a = a,
    mu = mean(d$mu_gamb), infectivity_mean = mean(d$infectivity),
    total_M = mean(d$total_M_gamb), Sm = mean(d$Sm_gamb_count),
    Pm = mean(d$Pm_gamb_count), Im = mean(d$Im_gamb_count),
    S = mean(d$S_count), D = mean(d$D_count), A = mean(d$A_count),
    U = mean(d$U_count), T = mean(d$Tr_count),
    n_infections = sum(d$n_infections))
  # by-band period sums (person-days) for rendered age outputs
  bb <- data.frame(band = BAND_NAMES, person_days = as.numeric(nage),
    inc_clin = as.numeric(inc("n_inc_clinical")), inc_sev = as.numeric(inc("n_inc_severe")),
    detect_lm = as.numeric(inc("n_detect_lm")), detect_pcr = as.numeric(inc("n_detect_pcr")))
  for (im in c("ib", "id", "ica", "iva", "icm", "ivm")) {
    cols <- grep(paste0("^", im, "_mean_"), names(d), value = TRUE)
    bb[[toupper(im)]] <- sapply(cols, function(cc) mean(d[[cc]]))
  }
  list(period = out, by_band = bb)
}

run_malsim_forensic <- function(eir, seed, N, years, chunk = CHUNK, p25 = FALSE) {
  set.seed(seed)
  p <- malsim_params(eir, N, p25)
  eq <- list(total_M = p$total_M, init_foim = p$init_foim, init_EIR = p$init_EIR,
             a = p$Q0 * p$blood_meal_rates, mum = p$mum,
             init_ib = p$init_ib, init_ica = p$init_ica, init_id = p$init_id,
             init_icm = p$init_icm, init_iva = p$init_iva, init_ivm = p$init_ivm)
  r <- run_resumable_simulation(1, p)
  snaps <- list(malsim_snapshot(r$state, p))
  per1 <- malsim_period(r$data, p)
  periods <- list(per1$period); bands <- list(cbind(per1$by_band, period = 0L))
  ends <- c(seq(chunk, years * 365, by = chunk))
  st <- r$state
  for (i in seq_along(ends)) {
    r <- run_resumable_simulation(ends[i], p, initial_state = st, restore_random_state = TRUE)
    st <- r$state
    d <- r$data[r$data$timestep > (if (i == 1) 1 else ends[i - 1]), ]
    pp <- malsim_period(d, p)
    periods[[i + 1]] <- pp$period
    bands[[i + 1]] <- cbind(pp$by_band, period = i)
    snaps[[i + 1]] <- malsim_snapshot(st, p)
  }
  list(model = if (p25) "malariasimulation_p25" else "malariasimulation", EIR = eir, seed = seed, N = N, eq = eq,
       snaps = snaps, periods = do.call(rbind, periods), period_bands = do.call(rbind, bands))
}

# ------------------------------------------------------------------ tasks ---

run_forensic <- function(model, eir, rep, seed, N, years, out_dir = "forensic/results",
                         chunk = CHUNK) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(out_dir, sprintf("%s_eir%s_rep%02d.rds", model, eir, rep))
  st <- system.time(res <- switch(model,
    magenta = run_magenta_forensic(eir, seed, N, years, chunk),
    malariasimulation = run_malsim_forensic(eir, seed, N, years, chunk),
    malariasimulation_p25 = run_malsim_forensic(eir, seed, N, years, chunk, p25 = TRUE)))
  res$rep <- rep; res$elapsed <- st[["elapsed"]]; res$host <- Sys.info()[["nodename"]]
  saveRDS(res, f)
  list(file = f, elapsed = st[["elapsed"]])
}
