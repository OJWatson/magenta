# Mosquito death and feeding rates under ITN/IRS, in closed form.
#
# The deterministic model in inst/extdata/odin_itn_irs.R computes its
# intervention block (coverage categories, net and spray decay, the repeat and
# successful feeding probabilities, zbar and wbar) from the time t and the
# parameters alone; no state variable feeds into mu or fv. pipeline() only
# needs those two rates, so this evaluates that block directly at the same
# daily time points mu_fv_create() returns, in the same order of operations as
# the generated C, which makes the result bit-identical to the odin output
# (tests/testthat/test-mu-fv.R checks this against the model).

#' Mosquito death and feeding rates from intervention coverage
#'
#' Evaluates the intervention block of the deterministic model (the one in
#' \code{inst/extdata/odin_itn_irs.R}) directly: the ITN and IRS decay, the
#' per-category feeding and repellency probabilities, and from them the
#' mosquito feeding rate \code{fv} and death rate \code{mu}, for every day of
#' the simulation. Neither rate depends on the model's state, so this gives
#' the same numbers as running the model and is what \code{\link{pipeline}}
#' uses.
#'
#' @param params The model parameter list from
#'   \code{\link{model_param_list_create}}, which holds the entomological and
#'   intervention parameters, or the equilibrium list from
#'   \code{\link{equilibrium_init_create}}, which carries it.
#' @param ft Treatment seeking, one value per year (only its length matters
#'   here; it is kept so the arguments match \code{\link{mu_fv_create}}).
#' @param itn_cov,irs_cov Bed net and spraying coverage, one value per year.
#'   A shorter vector is padded by repeating its last value.
#' @param years Length of the simulation in years.
#' @param int_times Days at which the yearly coverage values start. Default
#'   \code{NULL} uses \code{c(-25, 0, 365, 730, ...)}, i.e. the year-one
#'   coverage from before the simulation starts and a change on each
#'   anniversary.
#' @param itn_irs_on Day from which interventions act; \code{0} as in the
#'   pipeline, where every simulated day is on or after it.
#' @return Data frame with a row per day (\code{round(years * 365)} rows) and
#'   columns \code{mu} and \code{fv}.
#' @export
#' @examples
#' # nets from year two, spraying from year three, over four years
#' rates <- mu_fv_from_coverage(
#'   model_param_list_create(), ft = 0.4,
#'   itn_cov = c(0, 0.3, 0.6, 0.6), irs_cov = c(0, 0, 0.2, 0.2), years = 4
#' )
#' rates[c(1, 366, 731), ]
mu_fv_from_coverage <- function(params, ft, itn_cov, irs_cov, years,
                                int_times = NULL, itn_irs_on = 0) {
  n_years <- ceiling(years)
  first_then <- function(x, what) {
    x <- pad_annual(x, n_years, what)
    c(x[1], x)
  }
  ft <- first_then(ft, "ft")
  itn_vector <- first_then(itn_cov, "itn_cov")
  irs_vector <- first_then(irs_cov, "irs_cov")
  if (is.null(int_times)) {
    int_times <- c(-25, seq(0, by = 365, length.out = n_years))
  }
  if (length(int_times) != length(itn_vector)) {
    stop("int_times must have one entry per yearly value plus one", call. = FALSE)
  }
  p <- params
  num_int <- p$num_int

  tt <- seq(1, round(years * 365), 1)
  if (any(tt < int_times[1])) {
    stop("int_times must start on or before day 1", call. = FALSE)
  }
  # interpolate(t_vector, x, "constant"): the value at the last breakpoint on
  # or before t, held beyond the last one
  idx <- findInterval(tt, int_times)
  eff_itn <- itn_vector[idx]
  eff_irs <- irs_vector[idx]
  int_on <- int_times[idx]
  before <- tt < itn_irs_on
  eff_on <- ifelse(before, itn_irs_on, int_on)

  # net and spray decay since the start of the current coverage period
  itn_decay <- ifelse(before, 0, exp(-((tt - eff_on) %% p$ITN_interval) * p$itn_loss))
  irs_decay <- ifelse(before, 0, exp(-((tt - eff_on) %% p$IRS_interval) * p$irs_loss))
  d_itn <- ifelse(before, 0, p$d_ITN0 * itn_decay)
  r_itn <- ifelse(before, 0, p$r_ITN1 + (p$r_ITN0 - p$r_ITN1) * itn_decay)
  s_itn <- ifelse(before, 1, 1 - d_itn - r_itn)
  r_irs <- ifelse(before, 0, p$r_IRS0 * irs_decay)
  d_irs <- ifelse(before, 0, p$chi * p$d_IRS0 * irs_decay)
  s_irs <- ifelse(before, 1, 1 - d_irs)

  bed <- p$bites_Bed
  indoors <- p$bites_Indoors
  # coverage of, probability of feeding in, and of being repelled by, each
  # intervention category: none, ITN only, IRS only, both
  cov <- cbind(
    (1 - eff_itn) * (1 - eff_irs), eff_itn * (1 - eff_irs),
    (1 - eff_itn) * eff_irs, eff_itn * eff_irs
  )
  w <- cbind(
    1,
    1 - bed + bed * s_itn,
    1 - indoors + indoors * (1 - r_irs) * s_irs,
    1 - indoors + bed * (1 - r_irs) * s_itn * s_irs + (indoors - bed) * (1 - r_irs) * s_irs
  )
  z <- cbind(
    0,
    bed * r_itn,
    indoors * r_irs,
    bed * (r_irs + (1 - r_irs) * r_itn) + (indoors - bed) * r_irs
  )
  # summed category by category, as the C sum does
  zh <- 0
  wh <- 0
  for (k in seq_len(num_int)) {
    zh <- zh + cov[, k] * z[, k]
    wh <- wh + cov[, k] * w[, k]
  }
  zh <- ifelse(before, 0, zh)
  wh <- ifelse(before, 1, wh)

  zbar <- p$Q0 * zh
  fv <- 1 / (p$tau1 / (1 - zbar) + p$tau2)
  wbar <- 1 - p$Q0 + p$Q0 * wh
  p1 <- wbar * p$p10 / (1 - zbar * p$p10)
  mu <- -fv * log(p1 * p$p2)
  data.frame(mu = mu, fv = fv)
}
