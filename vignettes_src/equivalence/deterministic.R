# Deterministic Griffin equilibrium by age band and mosquito quantities (fine ages)
suppressMessages(library(magenta))
FT <- 0.4
BANDS_LO <- c(0, 1, 2, 5, 10, 15, 20, 40) * 365
BAND_NAMES <- c("0-1", "1-2", "2-5", "5-10", "10-15", "15-20", "20-40", "40+")
det_one <- function(eir, ages = c(seq(0, 10, 0.1), seq(10.5, 100, 0.5)), label = "fine") {
  mpl <- model_param_list_create()
  eq <- equilibrium_init_create(age_vector = ages, het_brackets = 5, EIR = eir, ft = FT,
    model_param_list = mpl, country = NULL, admin = NULL, quiet = TRUE)
  s3 <- function(x) apply(x, c(1, 2), sum)          # age x het
  w <- eq$den %o% eq$het_wt                             # population weights age x het
  S <- s3(eq$init_S); D <- s3(eq$init_D); A <- s3(eq$init_A); U <- s3(eq$init_U)
  T <- s3(eq$init_T); P <- s3(eq$init_P)
  # states are proportions of the whole population, so per-person = state / w
  band <- factor(BAND_NAMES[findInterval(eq$age, BANDS_LO)], levels = BAND_NAMES)
  ag <- function(x) as.numeric(tapply(rowSums(x), band, sum))
  den_b <- as.numeric(tapply(eq$den, band, sum))
  psi <- 1 - eq$rho * exp(-eq$age / eq$a0)
  pi_w <- (psi %o% eq$rel_foi) * w          # biting weight (rel_foi = het zeta)
  infl <- eq$cD * D + eq$cT * T + eq$cU * U + eq$cA_eq * A   # each is weight-scaled
  # immunity means: weighted by population
  imm <- function(x) as.numeric(tapply(rowSums(s3(x) * w) , band, sum)) / den_b
  # note init_I* arrays hold per-person immunity in each age x het cell (num_int 1 used)
  imm1 <- function(x) { y <- x[, , 1]; as.numeric(tapply(rowSums(y * w), band, sum)) / den_b }
  clin <- eq$phi_eq * eq$FOI * (S + A + U)
  bites <- eq$EIR_eq * w
  bb <- data.frame(band = BAND_NAMES, den = den_b,
    S = ag(S) / den_b, D = ag(D) / den_b, A = ag(A) / den_b, U = ag(U) / den_b,
    T = ag(T) / den_b, P = ag(P) / den_b,
    prev_micro = (ag(D) + ag(T) + ag(A * eq$p_det_eq)) / den_b,
    prev_any = (ag(D) + ag(T) + ag(A) + ag(U)) / den_b,
    clin_per_day = ag(clin) / den_b,
    IB = imm1(eq$init_IB), ICA = imm1(eq$init_ICA), ICM = imm1(eq$init_ICM),
    ID = imm1(eq$init_ID), IVA = imm1(eq$init_IVA),
    inf_bites_per_day = ag(bites) / den_b,
    pi_share = as.numeric(tapply(rowSums(pi_w), band, sum)) / sum(pi_w),
    kappa_contrib = as.numeric(tapply(rowSums(psi * infl %*% diag(eq$rel_foi)), band, sum)) / sum(pi_w))
  IVM_c <- outer(exp(-eq$age / mpl$dVM), eq$IVM_init_eq)
  cm <- function(y) as.numeric(tapply(rowSums(y * w), band, sum)) / den_b
  IBc <- eq$init_IB[, , 1]; ICc <- eq$init_ICA[, , 1] + eq$init_ICM[, , 1]
  IVAc <- eq$init_IVA[, , 1]
  bb$IVM <- cm(IVM_c)
  bb$b_inf <- cm(mpl$b0 * (mpl$b1 + (1 - mpl$b1) / (1 + (IBc / mpl$IB0)^mpl$kB)))
  bb$phi_clin <- cm(eq$phi_eq)
  fv_ <- 1 - (1 - mpl$fV0) / (1 + (eq$age / mpl$aV)^mpl$gammaV)
  thetac <- mpl$theta0 * (mpl$theta1 + (1 - mpl$theta1) /
      (1 + fv_ * ((ifelse(IVAc > 0, IVAc + 0.5, IVAc) + IVM_c) / mpl$IV0)^mpl$kV))
  bb$theta_sev <- cm(thetac)
  sevc <- thetac * eq$FOI * (S + A + U)     # expected severe cases per day, as clin
  u5_ <- eq$age < 5 * 365
  bb$q_det <- cm(eq$p_det_eq)
  glob <- data.frame(EIR = eir, ages = label, omega = eq$omega, mv0 = eq$mv0,
    Sv = eq$init_Sv, Ev = eq$init_Ev, Iv = eq$init_Iv, FOIv = eq$FOIv_eq,
    kappa = eq$FOIv_eq / eq$av0,
    kappa_direct = sum(psi * infl %*% diag(eq$rel_foi)) / sum(pi_w),
    av0 = eq$av0, mu0 = eq$mu0, EIR_check = 365 * eq$mv0 * eq$av0 * eq$init_Iv / eq$omega,
    pfpr_2_10 = with(bb, NA), clin_all = 365 * sum(clin) / sum(eq$den),
    sev_all = 365 * sum(sevc) / sum(eq$den), sev_u5 = 365 * sum(rowSums(sevc)[u5_]) / sum(eq$den[u5_]),
    S = sum(S) / sum(w), D = sum(D) / sum(w), A = sum(A) / sum(w), U = sum(U) / sum(w),
    T = sum(T) / sum(w), P = sum(P) / sum(w))
  list(glob = glob, bands = cbind(EIR = eir, ages = label, bb))
}
res <- list()
for (e in c(0.5, 1, 2, 4, 8, 16, 32, 64)) res[[length(res) + 1]] <- det_one(e)
g <- do.call(rbind, lapply(res, `[[`, "glob")); b <- do.call(rbind, lapply(res, `[[`, "bands"))
# PfPR2-10 microscopy and clinical incidence under 5 from fine ages (same as forensic quantities)
ages <- c(seq(0, 10, 0.1), seq(10.5, 100, 0.5))
g$pfpr_2_10 <- NA_real_; g$clin_u5 <- NA_real_
for (k in seq_along(res)) {
  e <- res[[k]]; eir <- e$glob$EIR
  mpl <- model_param_list_create()
  eq <- equilibrium_init_create(age_vector = ages, het_brackets = 5, EIR = eir, ft = FT,
    model_param_list = mpl, country = NULL, admin = NULL, quiet = TRUE)
  s3 <- function(x) apply(x, c(1, 2), sum)
  S <- rowSums(s3(eq$init_S)); D <- rowSums(s3(eq$init_D)); A <- rowSums(s3(eq$init_A))
  U <- rowSums(s3(eq$init_U)); T <- rowSums(s3(eq$init_T))
  clin <- rowSums(eq$phi_eq * eq$FOI * (s3(eq$init_S) + s3(eq$init_A) + s3(eq$init_U)))
  a <- eq$age; i2 <- a >= 2 * 365 & a < 10 * 365; u5 <- a < 5 * 365
  g$pfpr_2_10[k] <- sum(D[i2] + T[i2] + A[i2]) / sum(eq$den[i2])
  g$clin_u5[k] <- 365 * sum(clin[u5]) / sum(eq$den[u5])
}
write.csv(g, "deterministic_global.csv", row.names = FALSE)
write.csv(b, "deterministic_by_band.csv", row.names = FALSE)
saveRDS(list(global = g, bands = b), "deterministic.rds")
print(g[, c("EIR","mv0","kappa","pfpr_2_10","clin_all","clin_u5")], digits = 4)
