#------------------------------------------------
#' Model Parameter List Creation
#'
#' \code{model_param_list_create} creates list of model parameters to be used
#' within \code{equilibrium_init_create}
#'
#' @param eta Death rate for exponential population distribution, i.e. 1/Mean
#'   Population Age. Default = 0.0001305
#' @param rho Age-dependent biting parameter. Default = 0.85
#' @param a0 Age-dependent biting parameter. Default = 2920
#' @param sigma2 Variance of the log heterogeneity in biting rates. Default = 1.67
#' @param max_age Maximum age in days. Default = 100*365
#' @param rA Rate of leaving asymptomatic infection. Default = 0.00512821
#' @param rT Rate of leaving treatment. Default = 0.2
#' @param rD Rate of leaving clinical disease. Default = 0.2
#' @param rU Rate of recovering from subpatent infection. Default = 0.00906627
#' @param rP Rate of leaving prophylaxis. Default = 0.05
#' @param dE Latent period of human infection. Default = 12
#' @param delayGam Lag from parasites to infectious gametocytes. Default = 12.5
#' @param cD Untreated disease contribution to infectiousness. Default = 0.0676909
#' @param cT Treated disease contribution to infectiousness. Default =   0.322 * cD
#' @param cU Subpatent disease contribution to infectiousness. Default = 0.006203
#' @param gamma1 Parameter for infectiousness of state A. Default = 1.82425
#' @param d1 Minimum probability due to maximum immunity. Default = 0.160527
#' @param dID Inverse of decay rate. Default = 3650
#' @param ID0 Scale parameter. Default = 1.577533
#' @param kD Shape parameter. Default = 0.476614
#' @param uD Duration in which immunity is not boosted. Default = 9.44512
#' @param aD Scale parameter relating age to immunity. Default = 8001.99
#' @param fD0 Time-scale at which immunity changes with age. Default = 0.007055
#' @param gammaD Shape parameter relating age to immunity. Default = 4.8183
#' @param alphaA PCR detection probability parameters state A. Default = 0.757
#' @param alphaU PCR detection probability parameters state U. Default = 0.186
#' @param b0 Maximum probability due to no immunity. Default = 0.590076
#' @param b1 Maximum relative reduction due to immunity. Default = 0.5
#' @param dB Inverse of decay rate. Default = 3650
#' @param IB0 Scale parameter. Default = 43.8787
#' @param kB Shape parameter. Default = 2.15506
#' @param uB Duration in which immunity is not boosted. Default = 7.19919
#' @param phi0 Maximum probability due to no immunity. Default = 0.791666
#' @param phi1 Maximum relative reduction due to immunity. Default = 0.000737
#' @param dCA Inverse of decay rate. Default = 10950
#' @param IC0 Scale parameter. Default = 18.02366
#' @param kC Shape parameter. Default = 2.36949
#' @param uCA Duration in which immunity is not boosted. Default = 6.06349
#' @param PM New-born immunity relative to mother. Default = 0.774368
#' @param dCM Inverse of decay rate of maternal immunity. Default = 67.6952
#' @param delayMos Extrinsic incubation period. Default = 10
#' @param tau1 Duration of host seeking, assumed to be constant between species. Default = 0.69
#' @param tau2 Duration of mosquito resting after feed. Default = 2.31
#' @param mu0 Daily mortality of adult mosquitos. Default = 0.132
#' @param Q0 Anthrophagy probability. Default = 0.92
#' @param chi Endophily probability. Default = 0.86
#' @param bites_Bed Percentage of bites indoors and in bed. Default = 0.89
#' @param bites_Indoors Percentage of bites indoors . Default = 0.97
#' @param muEL Per capita daily mortality rate of early stage larvae (low density). Default = 0.0338
#' @param muLL Per capita daily mortality rate of late stage larvae (low density). Default = 0.0348
#' @param muPL Per capita daily mortality rate of pupae. Default = 0.249
#' @param dEL Development time of early stage larvae. Default = 6.64
#' @param dLL Development time of late stage larvae. Default = 3.72
#' @param dPL Development time of pupae. Default = 0.643
#' @param gammaL Relative effect of density dependence on late instars relative
#'   to early instars. Default = 13.25
#' @param km Seasonal carrying capacity. Default = 11
#' @param cm Seasonal birth rate. Default = 0.05
#' @param betaL Number of eggs laid per day per mosquito. Default = 21.2
#' @param num_int Number of possible interventions. Default = 4
#' @param itn_cov The proportion of people that use an ITN. Default = 0
#' @param irs_cov The proportion of people living in houses that have been sprayed. Default = 0
#' @param ITN_IRS_on Time of ITN and IRS to be activated. Default = -1, i.e. never.
#' @param DY Duration of year (days). Default = 365
#' @param d_ITN0 Probability of dying with an encounter with ITN (max). Default = 0.41
#' @param r_ITN0 Probability of repeating behaviour with ITN (max). Default = 0.56
#' @param r_ITN1 Probability of repeating behaviour with ITN (min). Default = 0.24
#' @param r_IRS0 Probability of repeating behaviour with IRS (min). Default = 0.6
#' @param d_IRS0 Probability of dying with an encounter with IRS (max). Default = 1
#' @param irs_half_life IRS half life. Default =   0.5 * DY
#' @param itn_half_life ITN half life. Default =   2.64 * DY
#' @param IRS_interval How long before IRS is repeated, i.e. when IRS decay = 1. Default =   1 * DY
#' @param ITN_interval How long before ITN is repeated, i.e. when IRS decay = 1.  Default =   3 * DY
#'
#'
#' @return A \code{magenta_core_params}: a list with the fields the
#'   simulation and the equilibrium solution read, see \code{\link{cpp_contract}}.
#' @export
#' @examples
#' mpl <- model_param_list_create(itn_cov = 0.3)
#' mpl$cov
#' mpl$itn_loss
#' validate_model_param_list(mpl)

model_param_list_create <- function(
  # age, heterogeneity in exposure,
  eta = 0.0001305,
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100 * 365,
  #  rate of leaving infection states
  rA = 0.00512821,
  rT = 0.2,
  rD = 0.2,
  rU = 0.00906627,
  rP = 1 / 25,
  #  human latent period and time lag from asexual parasites to
  dE = 12,
  delayGam = 12.5,
  # human infectiousness to mosquitoes
  cD = 0.0676909,
  cT = 0.322 * cD,
  cU = 0.006203,
  gamma1 = 1.82425,
  #  Immunity reducing probability of detection
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  # Immunity reducing probability of infection
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  # Immunity reducing probability of clinical disease
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  # entomological parameters
  delayMos = 10,
  tau1 = 0.69,
  tau2 = 2.31,
  mu0 = 0.132,
  Q0 = 0.92,
  chi = 0.86,
  bites_Bed = 0.89,
  bites_Indoors = 0.97,
  # larval parameters daily density dependent mortality rate of egg
  muEL = 0.0338,
  muLL = 0.0348,
  muPL = 0.249,
  dEL = 6.64,
  dLL = 3.72,
  dPL = 0.643,
  gammaL = 13.25,
  km = 11,
  cm = 0.05,
  betaL = 21.2,
  # intervention parameters
  num_int = 4,
  itn_cov = 0,
  irs_cov = 0,
  ITN_IRS_on = -1,
  DY = 365,
  d_ITN0 = 0.41,
  r_ITN0 = 0.56,
  r_ITN1 = 0.24,
  r_IRS0 = 0.6,
  d_IRS0 = 1,
  irs_half_life = 0.5 * DY,
  itn_half_life = 2.64 * DY,
  IRS_interval = 1 * DY,
  ITN_interval = 3 * DY
) {
  mpl <- as.list(environment())

  # derived mosquito quantities
  mpl$fv0 <- 1 / (tau1 + tau2)
  mpl$av0 <- Q0 * mpl$fv0 # daily feeeding rate on humans
  mpl$Surv0 <- exp(-mu0 * delayMos) # probability of surviving incubation period
  mpl$p10 <- exp(-mu0 * tau1) # probability of surviving one feeding cycle
  mpl$p2 <- exp(-mu0 * tau2) # probability of surviving one resting cycle

  # derived larval quantities {White et al. 2011 Parasites and Vectors}
  mpl$eov <- betaL / mu0 * (exp(mu0 / mpl$fv0) - 1)
  mpl$b_lambda <- (gammaL * muLL / muEL - dEL / dLL + (gammaL - 1) * muLL * dEL)
  mpl$lambda <- -0.5 * mpl$b_lambda +
    sqrt(0.25 * mpl$b_lambda^2 +
      gammaL * betaL * muLL * dEL / (2 * muEL * mu0 * dLL * (1 + dPL * muPL)))

  # derived intervention quantities
  mpl$pop_split <- rep(1 / num_int, num_int)
  # {No intervention} {ITN only} {IRS only} {Both ITN and IRS}
  mpl$cov <- c(
    (1 - itn_cov) * (1 - irs_cov), itn_cov * (1 - irs_cov),
    (1 - itn_cov) * irs_cov, itn_cov * irs_cov
  )
  mpl$irs_loss <- log(2) / irs_half_life
  mpl$itn_loss <- log(2) / itn_half_life

  # order the fields as the contract lists them
  mpl <- mpl[cpp_contract$core$fields]
  validate_model_param_list(new_param_list(mpl, "magenta_core_params"))
}

#' @rdname model_param_list_create
#' @param x List to validate.
#' @export
validate_model_param_list <- function(x) {
  what <- "model_param_list"
  check_fields(x, cpp_contract$core, what)
  vectors <- c("pop_split", "cov")
  for (field in setdiff(cpp_contract$core$fields, vectors)) {
    check_number(x[[field]], paste0(what, "$", field))
  }
  for (field in c("itn_cov", "irs_cov")) {
    check_number(x[[field]], paste0(what, "$", field), min = 0, max = 1)
  }
  check_number(x$pop_split, paste0(what, "$pop_split"), min = 0, max = 1, n = x$num_int)
  check_number(x$cov, paste0(what, "$cov"), min = 0, max = 1, n = 4)
  invisible(x)
}
