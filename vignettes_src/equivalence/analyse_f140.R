# Final magenta 1.4.0 vs malariasimulation 3.0.0 equivalence sweep (analysis of saved rds only).
# Run from final_140/:  Rscript analysis/analyse_f140.R [W0_years] [W1_years]
suppressMessages(library(ggplot2))
args <- commandArgs(TRUE)
Y0 <- if (length(args) >= 1) as.numeric(args[1]) else 4
Y1 <- if (length(args) >= 2) as.numeric(args[2]) else 10
out <- "analysis"; CH <- 73
BANDS <- c("0-1", "1-2", "2-5", "5-10", "10-15", "15-20", "20-40", "40+")
COL <- c(magenta = "#c51b7d", ms_p25 = "#2166ac", ms = "#4d9221", deterministic = "grey40")
MODS <- c("magenta", "ms_p25", "ms")
MLAB <- c(magenta = "magenta", ms_p25 = "malariasimulation, 25-day prophylaxis (primary)",
          ms = "malariasimulation, as released", deterministic = "deterministic equilibrium")
EIRS <- c(0.5, 1, 2, 4, 8, 16, 32, 64)
BC <- c("S", "D", "A", "U", "T", "P", "prev_micro", "IB", "ICA", "ICM", "ID", "IVA", "IVM",
        "b_inf", "phi_clin", "theta_sev", "q_det")

det <- read.csv("deterministic_global.csv"); det <- det[det$ages == "fine", ]
detb <- read.csv("deterministic_by_band.csv"); detb <- detb[detb$ages == "fine", ]
detb$ICtot <- detb$ICA + detb$ICM; detb$IVtot <- detb$IVA + detb$IVM
mu <- 0.132; a <- 0.92 / 3; fv <- 1 / 3; Q0 <- 0.92; tau <- 10
s <- exp(-mu); p_d <- 1 - s
z_cont <- function(k) a * k / (a * k + mu) * exp(-mu * tau)
z_disc <- function(k) { q <- fv * Q0 * k; q / (q + p_d) * (s^10 + s^11) / 2 }
det$zk <- det$Iv / z_cont(det$kappa)

# ------------------------------------------------------------------ read ---
files <- list.files("results", pattern = "\\.rds$", full.names = TRUE)
cat(length(files), "files\n")
glob <- list(); bnd <- list(); rt <- list()
for (f in files) {
  r <- readRDS(f)
  mod <- switch(r$model, magenta = "magenta", malariasimulation = "ms", malariasimulation_p25 = "ms_p25")
  K <- 50; N <- r$N
  sn <- r$snaps[2:(K + 1)]; pr <- r$periods
  if (mod != "magenta") pr <- pr[2:(K + 1), ]
  stopifnot(nrow(pr) == K)
  pf <- sapply(sn, function(x) { b <- x$by_band; i <- b$band %in% c("2-5", "5-10")
    sum(b$prev_micro[i] * b$n[i]) / sum(b$n[i]) })
  omega <- sapply(sn, `[[`, "omega")
  if (mod == "magenta") {
    d <- data.frame(EIR_adult = 365 * pr$inf_bites / (pr$days * N) / omega,
      m = r$snaps[[1]]$mos_n / N,
      kappa = sapply(sn, function(x) if (is.null(x$kappa_lag)) NA_real_ else x$kappa_lag),
      z = pr$inf_bites / pr$bites,
      clin_all = pr$inc_all * 365, clin_u5 = pr$inc_05 * 365,
      sev_all = pr$sev_all * 365, sev_u5 = pr$sev_05 * 365,
      S = pr$S, D = pr$D, A = pr$A, U = pr$U, T = pr$T, P = pr$P)
  } else {
    pb <- r$period_bands; pb <- pb[pb$period >= 1, ]
    u5 <- pb$band %in% BANDS[1:3]
    agg <- function(v, i) tapply(pb[[v]][i], pb$period[i], sum)
    al <- rep(TRUE, nrow(pb))
    d <- data.frame(EIR_adult = pr$EIR_adult, m = pr$total_M / N * det$omega[1],
      kappa = pr$kappa_lag, z = pr$Im / pr$total_M,
      clin_all = 365 * agg("inc_clin", al) / agg("person_days", al),
      clin_u5 = 365 * agg("inc_clin", u5) / agg("person_days", u5),
      sev_all = 365 * agg("inc_sev", al) / agg("person_days", al),
      sev_u5 = 365 * agg("inc_sev", u5) / agg("person_days", u5),
      S = pr$S / N, D = pr$D / N, A = pr$A / N, U = pr$U / N, T = pr$T / N, P = NA_real_)
  }
  d$pfpr_2_10 <- pf; d$chunk <- 1:K; d$model <- mod; d$EIR <- r$EIR; d$rep <- r$rep
  d$zk <- if (mod == "magenta") d$z / z_disc(d$kappa) else d$z / z_cont(d$kappa)
  glob[[length(glob) + 1]] <- d
  b <- do.call(rbind, lapply(0:K, function(k) {
    x <- r$snaps[[k + 1]]$by_band; cbind(x[c("band", "n", BC)], chunk = k) }))
  b$model <- mod; b$EIR <- r$EIR; b$rep <- r$rep
  bnd[[length(bnd) + 1]] <- b
  rt[[length(rt) + 1]] <- data.frame(model = mod, EIR = r$EIR, rep = r$rep, elapsed_min = r$elapsed / 60, host = r$host)
  rm(r, sn, pr); invisible(gc(verbose = FALSE))
}
glob <- do.call(rbind, glob); bnd <- do.call(rbind, bnd)
bnd$ICtot <- bnd$ICA + bnd$ICM; bnd$IVtot <- bnd$IVA + bnd$IVM
glob$year <- (glob$chunk - 0.5) * CH / 365
bnd$year <- bnd$chunk * CH / 365          # snapshot at end of chunk
saveRDS(list(glob = glob, bnd = bnd), file.path(out, "f140_extracted.rds"))
rt <- do.call(rbind, rt)
rts <- do.call(rbind, lapply(split(rt, list(rt$model, rt$EIR), drop = TRUE), function(x)
  data.frame(model = x$model[1], EIR = x$EIR[1], n = nrow(x), mean_min = mean(x$elapsed_min),
             se_min = sd(x$elapsed_min) / sqrt(nrow(x)), max_min = max(x$elapsed_min))))
write.csv(rts[order(rts$model, rts$EIR), ], file.path(out, "f140_runtime.csv"), row.names = FALSE)

# ------------------------------------------------------------- helpers ---
msd <- function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
inwin <- function(df) df$year > Y0 & df$year < Y1
pctdiff <- function(m1, s1, m2, s2) {
  R <- m1 / m2; se <- abs(R) * sqrt((s1 / m1)^2 + (s2 / m2)^2)
  c(pct = 100 * (R - 1), lo = 100 * (R - 1 - 1.96 * se), hi = 100 * (R - 1 + 1.96 * se))
}
# one comparison row from a data frame with columns model, v (rep-level window means) and a deterministic value
cmp_row <- function(x, dv) {
  g <- lapply(MODS, function(mo) msd(x$v[x$model == mo])); names(g) <- MODS
  r <- data.frame(det = dv)
  for (mo in MODS) { r[[paste0("mean_", mo)]] <- g[[mo]][["mean"]]; r[[paste0("se_", mo)]] <- g[[mo]][["se"]] }
  for (cm in list(c("magenta", "ms_p25"), c("magenta", "ms"), c("ms_p25", "ms"))) {
    nm <- paste(cm, collapse = "_vs_")
    A <- g[[cm[1]]]; B <- g[[cm[2]]]
    p <- if (all(is.finite(c(A, B))) && B[["mean"]] != 0 && A[["mean"]] != 0) pctdiff(A[["mean"]], A[["se"]], B[["mean"]], B[["se"]]) else c(pct = NA, lo = NA, hi = NA)
    r[[paste0("pct_", nm)]] <- p[["pct"]]; r[[paste0("lo_", nm)]] <- p[["lo"]]; r[[paste0("hi_", nm)]] <- p[["hi"]]
    r[[paste0("signif_", nm)]] <- p[["lo"]] > 0 | p[["hi"]] < 0
  }
  for (mo in MODS) { r[[paste0("pct_", mo, "_vs_det")]] <- 100 * (g[[mo]][["mean"]] / dv - 1)
                     r[[paste0("pctse_", mo, "_vs_det")]] <- 100 * g[[mo]][["se"]] / dv }
  r
}

# ------------------------------------------------ burn-in + drift ---
bi_cols <- c(EIR_realised = "EIR_adult", pfpr_2_10 = "pfpr_2_10", kappa = "kappa", clin_all = "clin_all", sev_all = "sev_all")
dref <- function(e, o) switch(o, EIR_realised = e, pfpr_2_10 = det$pfpr_2_10[det$EIR == e], kappa = det$kappa[det$EIR == e],
                              clin_all = det$clin_all[det$EIR == e], sev_all = det$sev_all[det$EIR == e])
bi <- do.call(rbind, lapply(names(bi_cols), function(nm) {
  x <- aggregate(glob[[bi_cols[[nm]]]], glob[c("model", "EIR", "chunk", "year")], function(z) c(mean(z, na.rm = TRUE), sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))))
  data.frame(x[1:4], outcome = nm, mean = x$x[, 1], se = x$x[, 2])
}))
for (v in c("ICA", "IVA")) {                       # slow immunity variables, adults 20-40
  g <- bnd[bnd$band == "20-40", ]
  x <- aggregate(g[[v]], g[c("model", "EIR", "chunk", "year")], function(z) c(mean(z), sd(z) / sqrt(length(z))))
  bi <- rbind(bi, data.frame(x[1:4], outcome = paste0(v, "_20-40"), mean = x$x[, 1], se = x$x[, 2]))
}
bi$ref <- mapply(function(e, o) if (grepl("_20-40", o)) detb[[sub("_20-40", "", o)]][detb$EIR == e & detb$band == "20-40"] else dref(e, o), bi$EIR, bi$outcome)
bi$ratio <- bi$mean / bi$ref; bi$ratio_se <- bi$se / bi$ref
write.csv(bi, file.path(out, "f140_burnin.csv"), row.names = FALSE)
W <- glob[inwin(glob), ]
drift <- do.call(rbind, lapply(split(W, list(W$model, W$EIR), drop = TRUE), function(x)
  do.call(rbind, lapply(names(bi_cols), function(nm) {
    v <- bi_cols[[nm]]; mid <- (Y0 + Y1) / 2
    sl <- sapply(split(x, x$rep), function(y) coef(lm(y[[v]] ~ y$year))[2] / mean(y[[v]]) * 100)
    data.frame(model = x$model[1], EIR = x$EIR[1], outcome = nm, slope_pct_per_yr = mean(sl), slope_se = sd(sl) / sqrt(length(sl)),
               z = mean(sl) / (sd(sl) / sqrt(length(sl))))
  }))))
write.csv(drift, file.path(out, "f140_burnin_drift.csv"), row.names = FALSE)

# ------------------------------------------------ headline ---
OUT <- c(EIR_adult = "realised_EIR", m = "m", kappa = "kappa", z = "sporozoite_rate", zk = "z_given_kappa",
         pfpr_2_10 = "PfPR2-10", clin_all = "clinical_all", clin_u5 = "clinical_u5",
         sev_all = "severe_all", sev_u5 = "severe_u5",
         S = "S", D = "D", A = "A", U = "U", T = "T", P = "P")
pr <- aggregate(W[names(OUT)], W[c("model", "EIR", "rep")], function(z) mean(z, na.rm = TRUE))
detv <- function(e, o) { d <- det[det$EIR == e, ]
  switch(o, EIR_adult = e, m = d$mv0, kappa = d$kappa, z = d$Iv, zk = d$zk, pfpr_2_10 = d$pfpr_2_10,
         clin_all = d$clin_all, clin_u5 = d$clin_u5, sev_all = d$sev_all, sev_u5 = d$sev_u5,
         S = d$S, D = d$D, A = d$A, U = d$U, T = d$T, P = d$P) }
rows <- function(vars) do.call(rbind, lapply(EIRS, function(e) do.call(rbind, lapply(vars, function(o) {
  cbind(data.frame(EIR = e, outcome = OUT[[o]]), cmp_row(data.frame(model = pr$model[pr$EIR == e], v = pr[[o]][pr$EIR == e]), detv(e, o))) }))))
h1 <- rows(names(OUT)[1:10]); h2 <- rows(names(OUT)[11:16])
write.csv(h1, file.path(out, "f140_headline.csv"), row.names = FALSE)

# ------------------------------------------------ bands: immunity, states ---
WB <- bnd[bnd$chunk >= 1 & inwin(bnd), ]
IMMV <- c("IB", "ICA", "ICM", "ICtot", "ID", "IVA", "IVM", "IVtot", "b_inf", "phi_clin", "theta_sev", "q_det")
STV <- c("S", "D", "A", "U", "T", "P", "prev_micro")
bp <- aggregate(WB[c(IMMV, STV)], WB[c("model", "EIR", "rep", "band")], mean)
bandrows <- function(vars) do.call(rbind, lapply(EIRS, function(e) do.call(rbind, lapply(BANDS, function(b) do.call(rbind, lapply(vars, function(o) {
  i <- bp$EIR == e & bp$band == b
  cbind(data.frame(EIR = e, band = b, variable = o), cmp_row(data.frame(model = bp$model[i], v = bp[[o]][i]), detb[[o]][detb$EIR == e & detb$band == b])) }))))))
hi <- bandrows(IMMV); hs <- bandrows(STV)
hs$band <- factor(hs$band, BANDS); hi$band <- factor(hi$band, BANDS)
# maternal immunity is ~1e-12 or smaller from age 10; relative differences there are numerical noise
mat0 <- hi$variable %in% c("ICM", "IVM") & hi$band %in% BANDS[5:8]
hi[mat0, grep("^(pct|lo|hi|signif|pctse)_", names(hi))] <- NA
write.csv(hi, file.path(out, "f140_immunity_bands.csv"), row.names = FALSE)
write.csv(hs, file.path(out, "f140_states_bands.csv"), row.names = FALSE)

# ------------------------------------------------ compact rds ---
tidy <- function(df, keys) do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  r <- df[i, ]; k <- r[keys]
  rbind(cbind(k, model = c(MODS, "deterministic"),
    mean = c(r$mean_magenta, r$mean_ms_p25, r$mean_ms, r$det),
    lo = c(r$mean_magenta - 1.96 * r$se_magenta, r$mean_ms_p25 - 1.96 * r$se_ms_p25, r$mean_ms - 1.96 * r$se_ms, r$det),
    hi = c(r$mean_magenta + 1.96 * r$se_magenta, r$mean_ms_p25 + 1.96 * r$se_ms_p25, r$mean_ms + 1.96 * r$se_ms, r$det)))
}))
recode <- function(d) { d$model <- c(magenta = "magenta", ms_p25 = "malariasimulation_p25", ms = "malariasimulation", deterministic = "deterministic")[d$model]; rownames(d) <- NULL; d }
pcols <- c("pct_magenta_vs_ms_p25", "lo_magenta_vs_ms_p25", "hi_magenta_vs_ms_p25", "signif_magenta_vs_ms_p25",
           "pct_magenta_vs_ms", "lo_magenta_vs_ms", "hi_magenta_vs_ms", "signif_magenta_vs_ms")
ai <- audit_meta <- NULL
meta <- list(magenta_sha = "1d424e6813a559739f8a4dee288aac3ea42d6065", magenta_version = "1.4.0",
  malariasimulation = "v3.0.0", N = 1e5, reps = 10, years = 10, chunk_days = CH, window_years = c(Y0, Y1),
  drug = "clinical treatment coverage 0.4; drug efficacy 1; ms_p25 prophylaxis Weibull shape 1 scale 30.5 d; ms default AL shape 11.3 scale 10.6",
  models = MLAB, ci = "mean +/- 1.96 replicate SE", date = as.character(Sys.Date()))
es <- list(headline = recode(tidy(h1, c("EIR", "outcome"))), immunity = recode(tidy(hi, c("EIR", "band", "variable"))),
           states = recode(tidy(hs, c("EIR", "band", "variable"))),
           pct_headline = h1[c("EIR", "outcome", pcols)], pct_immunity = hi[c("EIR", "band", "variable", pcols)], meta = meta)
es$immunity$band <- as.character(es$immunity$band); es$states$band <- as.character(es$states$band)
saveRDS(es, file.path(out, "equivalence_summary.rds"), compress = "xz")

# ---------------------------------------------------------------- figures ---
th <- theme_bw(base_size = 9) + theme(strip.background = element_rect(fill = "grey95"), legend.position = "bottom", panel.grid.minor = element_blank())
olab <- c(realised_EIR = "Realised EIR", m = "Mosquito density m", kappa = "Infectivity kappa", PfPR2_10 = "PfPR2-10 (microscopy)",
          clinical_all = "Clinical incidence, all ages", clinical_u5 = "Clinical incidence, under 5",
          severe_all = "Severe incidence, all ages", severe_u5 = "Severe incidence, under 5")
okeys <- c("realised_EIR", "m", "kappa", "PfPR2-10", "clinical_all", "clinical_u5", "severe_all", "severe_u5")
olabv <- setNames(unname(olab), okeys)
xs <- scale_x_continuous(trans = "log2", breaks = EIRS, labels = EIRS)
mcol <- scale_colour_manual(values = COL, labels = MLAB, name = NULL)
pdg <- position_dodge(width = 0.08)
# (a)
da <- es$headline[es$headline$outcome %in% okeys, ]; da$outcome <- factor(olabv[da$outcome], levels = olab)
dd_ <- da[da$model == "deterministic", ]; dm <- da[da$model != "deterministic", ]
dm$model <- factor(dm$model, c("magenta", "malariasimulation_p25", "malariasimulation"))
cmap <- c(magenta = COL[["magenta"]], malariasimulation_p25 = COL[["ms_p25"]], malariasimulation = COL[["ms"]], deterministic = "grey40")
lmap <- c(magenta = MLAB[["magenta"]], malariasimulation_p25 = MLAB[["ms_p25"]], malariasimulation = MLAB[["ms"]], deterministic = MLAB[["deterministic"]])
pa <- ggplot(dm, aes(EIR, mean, colour = model)) +
  geom_line(data = dd_, aes(colour = model), linetype = "dashed", linewidth = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = pdg) + geom_line(position = pdg, linewidth = 0.4) +
  geom_point(position = pdg, size = 1.2) + facet_wrap(~outcome, scales = "free_y", ncol = 4) + xs +
  scale_colour_manual(values = cmap, labels = lmap, name = NULL) +
  labs(x = "Nominal EIR (per adult per year)", y = "Mean over window (+/- 1.96 SE)") + th + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_a_outcomes.png"), pa, width = 10, height = 6, dpi = 200)
# (b)
SER <- c(ms_p25 = "vs malariasimulation, 25-day prophylaxis (primary)", ms = "vs malariasimulation, as released")
sercol <- scale_colour_manual(values = COL[c("ms_p25", "ms")], labels = SER, name = NULL)
longcmp <- function(df, keys) rbind(
  cbind(df[keys], series = "ms_p25", pct = df$pct_magenta_vs_ms_p25, lo = df$lo_magenta_vs_ms_p25, hi = df$hi_magenta_vs_ms_p25, sig = df$signif_magenta_vs_ms_p25),
  cbind(df[keys], series = "ms", pct = df$pct_magenta_vs_ms, lo = df$lo_magenta_vs_ms, hi = df$hi_magenta_vs_ms, sig = df$signif_magenta_vs_ms))
db <- longcmp(h1[h1$outcome %in% okeys, ], c("EIR", "outcome")); db$outcome <- factor(olabv[db$outcome], levels = olab)
pb <- ggplot(db, aes(EIR, pct, colour = series)) + geom_hline(yintercept = 0, colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = pdg) + geom_line(position = pdg, linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 1.6, position = pdg) + sercol +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "95% CI excludes 0") + facet_wrap(~outcome, scales = "free_y", ncol = 4) + xs +
  labs(x = "Nominal EIR (per adult per year)", y = "magenta vs malariasimulation (% difference)") + th + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_b_pctdiff.png"), pb, width = 10, height = 6, dpi = 200)
# (c) immunity by band, EIR 2, 8, 32
immlab <- c(IB = "IB", ICA = "ICA", ICM = "ICM", ID = "ID", IVA = "IVA", IVM = "IVM")
plot_bands <- function(vars, eirs, fn, h, ylab) {
  d <- es$immunity[es$immunity$variable %in% vars & es$immunity$EIR %in% eirs, ]
  d <- d[!(d$variable %in% c("ICM", "IVM") & d$band %in% BANDS[5:8]), ]
  d$band <- factor(d$band, BANDS); d$variable <- factor(d$variable, vars); d$EIRl <- factor(paste("EIR", d$EIR), paste("EIR", eirs))
  dm <- d[d$model != "deterministic", ]; dd <- d[d$model == "deterministic", ]
  dm$model <- factor(dm$model, c("magenta", "malariasimulation_p25", "malariasimulation"))
  p <- ggplot(dm, aes(band, mean, colour = model, group = model)) +
    geom_line(data = dd, linetype = "dashed", linewidth = 0.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = position_dodge(width = 0.4)) +
    geom_line(linewidth = 0.4, position = position_dodge(width = 0.4)) + geom_point(size = 1, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = cmap, labels = lmap, name = NULL) + facet_grid(variable ~ EIRl, scales = "free_y") +
    labs(x = "Age band (years)", y = ylab) + th + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(colour = guide_legend(nrow = 2))
  ggsave(file.path(out, fn), p, width = 9, height = h, dpi = 200)
}
plot_bands(names(immlab), c(2, 8, 32), "fig_c_immunity_bands.png", 11, "Mean immunity (window average, +/- 1.96 SE)")
plot_bands(c("b_inf", "phi_clin", "theta_sev", "q_det"), c(2, 8, 32), "fig_f_prob_scales.png", 8, "Mean probability (window average)")
# (d)
vv <- c("IB", "ICA", "ICM", "ICtot", "ID", "IVA", "IVM", "IVtot", "b_inf", "phi_clin", "theta_sev", "q_det")
dd <- longcmp(hi[hi$EIR == 8 & hi$variable %in% vv, ], c("band", "variable")); dd$variable <- factor(dd$variable, vv)
pd <- ggplot(dd, aes(band, pct, colour = series, group = series)) + geom_hline(yintercept = 0, colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = position_dodge(width = 0.4)) +
  geom_line(linewidth = 0.4, position = position_dodge(width = 0.4)) + geom_point(size = 1.2, position = position_dodge(width = 0.4)) +
  sercol + facet_wrap(~variable, ncol = 4, scales = "free_y") +
  labs(x = "Age band (years)", y = "magenta vs malariasimulation (% difference), EIR 8") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_d_immunity_pctdiff_eir8.png"), pd, width = 10, height = 7, dpi = 200)
# (e) burn-in
dc <- bi[bi$model %in% MODS & bi$outcome %in% c("EIR_realised", "pfpr_2_10", "kappa", "clin_all", "sev_all", "ICA_20-40", "IVA_20-40"), ]
dc$EIRl <- factor(paste("EIR", dc$EIR), paste("EIR", EIRS)); dc$outcome <- factor(dc$outcome, unique(dc$outcome))
pe <- ggplot(dc, aes(year, ratio, colour = model, fill = model)) +
  annotate("rect", xmin = Y0, xmax = Y1, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed") +
  geom_ribbon(aes(ymin = ratio - 1.96 * ratio_se, ymax = ratio + 1.96 * ratio_se), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.35) +
  scale_colour_manual(values = COL[MODS], labels = MLAB[MODS], name = NULL) + scale_fill_manual(values = COL[MODS], guide = "none") +
  scale_x_continuous(breaks = c(0, 5, 10)) + facet_grid(outcome ~ EIRl, scales = "free_y") +
  labs(x = "Year", y = "Ratio to deterministic (mean of 10 reps; grey = window)") + th +
  theme(axis.title.y = element_text(size = 7), panel.spacing.x = unit(0.08, "in"), axis.text = element_text(size = 6)) + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_e_burnin.png"), pe, width = 12, height = 9, dpi = 200)

op <- options(width = 250)
cat("\nWindow", Y0, "-", Y1, "\n")
print(h1[h1$outcome %in% c("realised_EIR", "kappa", "PfPR2-10", "clinical_all", "clinical_u5", "severe_all", "severe_u5"),
  c("EIR", "outcome", "mean_magenta", "mean_ms_p25", "mean_ms", "det", "pct_magenta_vs_ms_p25", "signif_magenta_vs_ms_p25", "pct_magenta_vs_ms", "signif_magenta_vs_ms", "pct_magenta_vs_det", "pct_ms_p25_vs_det", "pct_ms_vs_det")], digits = 3, row.names = FALSE)
cat("\nDrift |z|>3\n"); print(drift[abs(drift$z) > 3, ], digits = 3, row.names = FALSE)
options(op)
