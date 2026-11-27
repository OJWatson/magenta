context("site")

itn_mu_fv <- function(...) {
  mpl <- model_param_list_create(...)
  eq <- equilibrium_init_create(
    age_vector = c(seq(0, 10, 0.5), seq(11, 100, 2)), het_brackets = 5,
    ft = 0.3, EIR = 20, country = "Kenya", admin = "Kisumu",
    model_param_list = mpl, quiet = TRUE
  )
  mu_fv_create(eq, ft = c(0.3, 0.4, 0.4, 0.5), itn_cov = c(0.2, 0.5, 0.5, 0.6),
               irs_cov = c(0, 0.1, 0.1, 0), years = 4, full = TRUE)
}

test_that("default itn decay clock gives the same mu and fv as before", {
  # reference made with magenta at 546ff8f, before the clock option
  ref <- readRDS(test_path("fixtures", "mu_fv_breakpoint.rds"))
  out <- suppressMessages(itn_mu_fv())
  expect_equal(out$mu, ref$mu, tolerance = 1e-10)
  expect_equal(out$fv, ref$fv, tolerance = 1e-10)

  # nets are new at every yearly breakpoint (0, 365, 730, 1095)
  loss <- log(2) / (2.64 * 365)
  t <- seq_along(out$ITN_decay)
  expect_equal(out$ITN_decay, exp(-(t - pmin(365 * floor(t / 365), 1095)) * loss))
})

test_that("campaign and constant itn decay clocks", {
  loss <- log(2) / (2.64 * 365)

  camp <- suppressMessages(itn_mu_fv(itn_decay_clock = "campaign", itn_campaign_t0 = 100))
  t <- seq_along(camp$ITN_decay)
  expect_equal(camp$ITN_decay, exp(-((t - 100) %% (3 * 365)) * loss))

  const <- suppressMessages(itn_mu_fv(itn_decay_clock = "constant"))
  k <- loss * 3 * 365
  expect_equal(unique(round(const$ITN_decay, 12)), round((1 - exp(-k)) / k, 12))

  # older nets kill fewer mosquitoes
  ref <- readRDS(test_path("fixtures", "mu_fv_breakpoint.rds"))
  expect_lt(mean(const$mu), mean(ref$mu))
})

test_that("itn and irs half-life arguments take effect", {
  def <- model_param_list_create()
  expect_equal(def$itn_loss, log(2) / (2.64 * 365))
  expect_equal(def$irs_loss, log(2) / (0.5 * 365))
  expect_equal(def$itn_clock_mode, 0)

  mpl <- model_param_list_create(itn_half_life = 900, irs_half_life = 100,
                                 itn_decay_clock = "constant")
  expect_equal(mpl$itn_loss, log(2) / 900)
  expect_equal(mpl$irs_loss, log(2) / 100)
  expect_equal(mpl$itn_clock_mode, 2)
})

test_that("site_ft", {
  dl <- drug_list_create(drugs = list(drug_create_al()))
  eff <- drug_create_al()$lpf[1]
  expect_equal(site_ft(c(0, 0.5), dl), c(0, 0.5 / eff))
  expect_equal(site_ft(0.99, dl), 1)

  # weighted by partner drug ratios
  dl2 <- drug_list_create(drugs = list(drug_create_al(), drug_create_dhappq()),
                          number_of_drugs = 2, partner_drug_ratios = c(0.25, 0.75))
  eff2 <- 0.25 * eff + 0.75 * drug_create_dhappq()$lpf[1]
  expect_equal(site_ft(0.4, dl2), 0.4 / eff2)
})

test_that("site_seasonality", {
  site <- readRDS(test_path("fixtures", "site_busia_urban.rds"))
  s <- site_seasonality(site)
  sp <- site$seasonality$seasonality_parameters
  expect_named(s, c("a0", "a1", "a2", "a3", "b1", "b2", "b3", "theta_c"))
  expect_equal(unname(s[1:7]), unlist(unname(sp[c("g0", "g1", "g2", "g3", "h1", "h2", "h3")])))
  # the normalised curve has mean 1 (no flooring for this unit)
  w <- 2 * pi * (1:365) / 365
  r <- s[["a0"]] + s[["a1"]] * cos(w) + s[["a2"]] * cos(2 * w) + s[["a3"]] * cos(3 * w) +
    s[["b1"]] * sin(w) + s[["b2"]] * sin(2 * w) + s[["b3"]] * sin(3 * w)
  expect_equal(mean(pmax(r / s[["theta_c"]], 0.001)), 1)

  # equilibrium_init_create uses it in place of the admin lookup
  eq <- equilibrium_init_create(age_vector = c(0, 1, 5, 10, 20, 40, 60), het_brackets = 3,
                                ft = 0, EIR = 10, model_param_list = model_param_list_create(),
                                seasonality = s, quiet = TRUE)
  expect_equal(eq$ssa0, s[["a0"]])
  expect_equal(eq$ssb3, s[["b3"]])
  expect_equal(eq$theta_c, s[["theta_c"]])
})

test_that("site_to_pipeline", {
  skip_if_not_installed("site")
  site <- readRDS(test_path("fixtures", "site_busia_urban.rds"))
  int <- site$interventions
  dl <- drug_list_create(drugs = list(drug_create_al()))
  a <- site_to_pipeline(site, burnin = 3, drug_list = dl)

  # every argument is one pipeline() or model_param_list_create() takes
  ok <- c(names(formals(pipeline)), names(formals(model_param_list_create)))
  expect_true(all(names(a) %in% ok))

  expect_equal(a$years, 8)
  expect_equal(attr(a, "year"), 1997:2004)
  expect_equal(a$EIR, site$eir$eir)
  expect_equal(a$seasonality, site_seasonality(site))

  # untreated burn-in, then the site's yearly values
  expect_equal(a$ft[1:3], c(0, 0, 0))
  expect_equal(a$itn_cov[1:3], c(0, 0, 0))
  expect_equal(a$ft[4:8], site_ft(int$treatment$implementation$tx_cov, dl))
  expect_equal(a$itn_cov[4:8], int$itn$use$itn_use)
  expect_equal(a$irs_cov[4:8], 0.75 * int$irs$implementation$irs_cov)
  expect_equal(attr(a, "tx_cov")[4:8], int$treatment$implementation$tx_cov)

  # burn-in at the first year's values
  b <- site_to_pipeline(site, burnin = 2, burnin_interventions = "first", drug_list = dl)
  expect_equal(b$itn_cov[1:3], rep(int$itn$use$itn_use[1], 3))

  # nets from site::net_efficacy at the use weighted resistance
  res <- site$vectors$pyrethroid_resistance
  r <- round(stats::weighted.mean(res$pyrethroid_resistance, int$itn$use$itn_use), 2)
  ne <- site::net_efficacy
  ne <- ne[ne$net_type == "pyrethroid_only" & round(ne$pyrethroid_resistance, 2) == r, ]
  expect_equal(c(a$d_ITN0, a$r_ITN0, a$r_ITN1), c(ne$dn0, ne$rn0, ne$rnm))
  expect_equal(a$itn_half_life, ne$gamman * log(2))
  expect_equal(a$ITN_interval, int$itn$retention_half_life)
  expect_equal(a$itn_decay_clock, "constant")

  # species weighted vector parameters
  v <- site$vectors$vector_species
  expect_equal(a$Q0, sum(v$prop * v$Q0))
  expect_equal(a$bites_Bed, sum(v$prop * v$phi_bednets))
  expect_equal(a$bites_Indoors, sum(v$prop * v$phi_indoors))
  expect_equal(a$mu0, sum(v$prop * v$mum))

  # and the parameters are accepted
  mpl <- do.call(model_param_list_create, a[names(a) %in% names(formals(model_param_list_create))])
  expect_equal(mpl$itn_clock_mode, 2)
  expect_equal(mpl$Q0, a$Q0)
})

test_that("core parameters and seasonality reach the individual model on a fresh start", {
  hk <- housekeeping_list_create(quiet = TRUE, cluster = TRUE)
  site <- readRDS(test_path("fixtures", "site_busia_urban.rds"))
  s <- site_seasonality(site)

  # defaults are the individual model's own
  def <- pipeline(EIR = 5, N = 200, years = 0.05, full_save = TRUE, seed = 1L,
                  housekeeping_list = hk)$parameters_List
  expect_equal(def$g_Q0, 0.92)
  expect_equal(def$g_theta, rep(1, 365))

  fs <- pipeline(EIR = 5, N = 200, years = 0.05, full_save = TRUE, seed = 1L,
                 housekeeping_list = hk, seasonality = s,
                 Q0 = 0.85, mu0 = 0.12, b0 = 0.55, sigma2 = 1.5, PM = 0.7)$parameters_List
  expect_equal(fs$g_Q0, 0.85)
  expect_equal(fs$g_mu0, 0.12)
  expect_equal(fs$g_b0, 0.55)
  expect_equal(fs$g_PM, 0.7)
  expect_equal(fs$g_zeta_sdlog, sqrt(1.5))
  expect_equal(fs$g_zeta_meanlog, -0.75)
  expect_equal(fs$g_theta, seasonal_theta(s))
  expect_equal(mean(fs$g_theta), 1)

  # the admin lookup is unchanged
  eq <- equilibrium_init_create(age_vector = c(0, 1, 5, 10, 20, 40, 60), het_brackets = 3,
                                ft = 0, EIR = 10, model_param_list = model_param_list_create(),
                                country = "Kenya", admin = "Kisumu", quiet = TRUE)
  expect_equal(equilibrium_ss_create(eq)$theta, seasonal_profile("Kisumu", "Kenya"))
})
