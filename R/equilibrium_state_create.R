#' Equilibrium list creation
#'
#' \code{equilibrium_ss_create} creates a starting state from
#'   the output of \code{equilibrium_init_create}
#'
#' @param eqInit Equilibrium initialisation state as
#'   output from \code{equilibrium_init_create}
#'
#' @keywords internal
equilibrium_ss_create <- function(eqInit) {
  ## create equilibrium state for return
  equilibrium_state <- list(
    "age_brackets" = eqInit$age_bounds,
    "het_brackets" = eqInit$het_bounds,
    "Smat" = eqInit$init_S[, , 1],
    "Dmat" = eqInit$init_D[, , 1],
    "Amat" = eqInit$init_A[, , 1],
    "Umat" = eqInit$init_U[, , 1],
    "Tmat" = eqInit$init_T[, , 1],
    "Pmat" = eqInit$init_P[, , 1],
    "IBmat" = eqInit$init_IB[, , 1],
    "ICAmat" = eqInit$init_ICA[, , 1],
    "ICMmat" = eqInit$init_ICM[, , 1],
    "IDmat" = eqInit$init_ID[, , 1],
    "FOI" = eqInit$FOI,
    "phi" = eqInit$phi_eq,
    "Sv" = eqInit$init_Sv * eqInit$mv0,
    "Ev" = eqInit$init_Ev * eqInit$mv0,
    "Iv" = eqInit$init_Iv * eqInit$mv0,
    "ICM_Init" = eqInit$ICM_init_eq,
    "theta" = seasonal_profile(eqInit$admin, eqInit$country)
  )

  return(equilibrium_state)
}


#' Spatial Parameter List Creation
#'
#' \code{spl_create} creates the user required spatial
#' list to be passed to the initialisation function if needed.
#'
#' @param spatial_type Numeric for spatial type. Default = 0
#' @param human_importation_rate_vector Vector of human importations into this
#'   region. Default = NULL
#' @param mosquito_imporation_rate_vector Vector of mosquito infection
#'   importations into this region. Default = NULL
#' @param cotransmission_freq_vector Vector of size of cotransmissions. This is
#'   preallocated here and regenerated every simulation update length.
#'   Default = rep(1, 10000)
#' @param oocyst_freq_vector Vector of number of ooocysts. This is
#'   preallocated here and regenerated every simulation update length.
#'   Default = rep(1, 10000)
#' @param num_human_infs Number of human infections from last day. Required for
#'   metapopulation. Default = 0
#' @param num_mos_infs Number of mosquito infections from last day. Required for
#'   metapopulation. Default = 0
#' @param plaf Population Level Frequency. Default = rep(0.5, 24)
#' @param island_imports_plaf_linked_flag Boolean. Whether imported barcodes for
#'  island model are drawn dependent on other sites. Currently, if TRUE, imported
#'  barcodes will either be all 0 or will be 1 at all loci where plaf > 0,
#'  if \code{rbinom(1, 1, plaf[plaf>0][1])} is TRUE. I.e. the first loci greater
#'  than 0 is used to determine if the incoming barcode is 1 at all sites where
#'  plaf is greater than 0.
#' @keywords internal
spl_create <- function(spatial_type = 0,
                       human_importation_rate_vector = NULL,
                       mosquito_imporation_rate_vector = NULL,
                       cotransmission_freq_vector = rep(1, 10000),
                       oocyst_freq_vector = rep(1, 10000),
                       num_human_infs = 0,
                       num_mos_infs = 0,
                       plaf = rep(0.5, 24),
                       island_imports_plaf_linked_flag = FALSE) {
  # if non spatial
  if (spatial_type == 0) {
    spatial_list <- list(
      "spatial_type" = spatial_type,
      "cotransmission_freq_vector" = cotransmission_freq_vector,
      "oocyst_freq_vector" = oocyst_freq_vector,
      "imported_cotransmissions_events" = 0,
      "imported_oocyst_events" = 0,
      "plaf" = plaf,
      "island_imports_plaf_linked_flag" = island_imports_plaf_linked_flag
    )

    # if island spatial
  } else if (spatial_type == 1) {
    # first work out how frequent importation is:
    hum_imp_rate <- sum(human_importation_rate_vector, na.rm = TRUE)
    mos_imp_rate <- sum(mosquito_imporation_rate_vector, na.rm = TRUE)

    spatial_list <- list(
      "spatial_type" = spatial_type,
      "cotransmission_freq_vector" = cotransmission_freq_vector,
      "oocyst_freq_vector" = oocyst_freq_vector,
      "imported_cotransmissions_events" = hum_imp_rate,
      "imported_oocyst_events" = mos_imp_rate,
      "plaf" = plaf,
      "island_imports_plaf_linked_flag" = island_imports_plaf_linked_flag
    )

    # metapopulation spatial lists were never finished
  } else {
    stop("metapopulation spatial lists are not implemented")
  }

  spatial_list
}


#' Odin generator function
#'
#' \code{generate_default_model} creates the user required generator for
#' the odin model.
#'
#' @param ft Treatment seeking numeric
#' @param age Vector of ages
#' @param generator generator for model.
#' @param dat Output of \code{equilibrium_init_create}
#' @param dde Use dde in soliving eventual odin. Default = FALSE
#'
#' @keywords internal
generate_default_model <- function(ft, age, dat, generator, dde = TRUE) {
  generator(
    init_S = dat$S,
    init_T = dat$T,
    init_D = dat$D,
    init_A = dat$A,
    init_U = dat$U,
    init_P = dat$P,
    init_ICA = dat$ICA,
    init_ICM = dat$ICM,
    init_ID = dat$ID,
    init_IB = dat$IB,
    init_Sv = dat$Sv,
    init_Ev = dat$Ev,
    init_Iv = dat$Iv,
    init_PL = dat$PL,
    init_LL = dat$LL,
    init_EL = dat$EL,
    na = dat$na,
    nh = dat$nh,
    age_rate = dat$age_rate,
    foi_age = dat$foi_age,
    het_wt = dat$het_wt,
    rel_foi = dat$rel_foi,
    omega = dat$omega,
    pi = pi,
    x_I = dat$x_I,
    mv0 = dat$mv0,
    mu0 = dat$mu0,
    ssa0 = dat$ssa0,
    ssa1 = dat$ssa1,
    ssa2 = dat$ssa2,
    ssa3 = dat$ssa3,
    ssb1 = dat$ssb1,
    ssb2 = dat$ssb2,
    ssb3 = dat$ssb3,
    theta_c = dat$theta_c,
    den = dat$den,
    age59 = dat$age59,
    age05 = dat$age05,
    age2years = dat$age2years,
    age10years = dat$age10years,
    two_to_10_length = dat$age10years - dat$age2years,
    age = age * 365,
    ft = ft,
    use_dde = dde,
    eta = dat$eta,
    rA = dat$rA,
    rT = dat$rT,
    rD = dat$rD,
    rU = dat$rU,
    rP = dat$rP,
    dE = dat$dE,
    delayGam = dat$delayGam,
    delayMos = dat$delayM,
    cD = dat$cD,
    cT = dat$cT,
    cU = dat$cU,
    gamma1 = dat$gamma1,
    d1 = dat$d1,
    dID = dat$dID,
    ID0 = dat$ID0,
    kD = dat$kD,
    uD = dat$uD,
    aD = dat$aD,
    fD0 = dat$fD0,
    gammaD = dat$gammaD,
    b0 = dat$b0,
    b1 = dat$b1,
    dB = dat$dB,
    IB0 = dat$IB0,
    kB = dat$kB,
    uB = dat$uB,
    phi0 = dat$phi0,
    phi1 = dat$phi1,
    dCA = dat$dCA,
    IC0 = dat$IC0,
    kC = dat$kC,
    uCA = dat$uCA,
    dCM = dat$dCM,
    tau1 = dat$tau1,
    tau2 = dat$tau2,
    Q0 = dat$Q0,
    chi = dat$chi,
    bites_Bed = dat$bites_Bed,
    bites_Indoors = dat$bites_Indoors,
    p10 = dat$p10,
    p2 = dat$p2,
    muEL = dat$muEL,
    muLL = dat$muLL,
    muPL = dat$muPL,
    dEL = dat$dEL,
    dLL = dat$dLL,
    dPL = dat$dPL,
    gammaL = dat$gammaL,
    beta_larval0 = dat$betaL,
    num_int = dat$num_int,
    itn_cov = dat$itn_cov,
    irs_cov = dat$irs_cov,
    int_times = dat$int_times,
    ITN_IRS_on = dat$ITN_IRS_on,
    d_ITN0 = dat$d_ITN0,
    r_ITN0 = dat$r_ITN0,
    r_ITN1 = dat$r_ITN1,
    r_IRS0 = dat$r_IRS0,
    d_IRS0 = dat$d_IRS0,
    IRS_interval = dat$IRS_interval,
    ITN_interval = dat$ITN_interval,
    irs_loss = dat$irs_loss,
    itn_loss = dat$itn_loss
  )
}
