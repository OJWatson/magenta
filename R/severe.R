#------------------------------------------------
#' Severe disease rates from simulation loggers
#'
#' \code{severe_rates} turns the \code{Loggers} returned by a simulation into
#' annual rates per person of clinical and severe disease, for all ages and
#' for children aged 0 to 5 years.
#'
#' The severe output is the expected number of severe cases, i.e. the sum of
#' the probability that each new infection is severe (Griffin et al. 2015
#' Proc R Soc B 282:20142657), rather than a random draw. It is averaged over the logging
#' period in the same way as \code{Incidence}. The denominators are the
#' population at the end of the logging period.
#'
#' @param loggers The \code{Loggers} list from a simulation, e.g.
#'   \code{pipeline(...)$Loggers}.
#'
#' @return A data frame with columns \code{age_group} ("all" and "0-5"),
#'   \code{clinical} and \code{severe}, both per person per year.
#' @export
severe_rates <- function(loggers) {

  if (!all(c("Severe", "Severe_05", "Incidence", "Incidence_05", "Ages") %in%
           names(loggers))) {
    stop("loggers must contain Severe, Severe_05, Incidence, Incidence_05 and Ages")
  }

  # same age cut off as the C++ loggers, i.e. <= 5 years in days
  n <- length(loggers$Ages)
  n_05 <- sum(loggers$Ages <= 5 * 365)

  data.frame(
    age_group = c("all", "0-5"),
    clinical = 365 * c(loggers$Incidence / n, loggers$Incidence_05 / n_05),
    severe = 365 * c(loggers$Severe / n, loggers$Severe_05 / n_05)
  )
}

#------------------------------------------------
#' Hospitalisations and deaths from severe disease
#'
#' \code{severe_outcomes} converts a rate of severe disease into rates of
#' hospitalised and community severe cases and deaths, following
#' \code{severe_incidence_mortality} in the postie package (1.1.0), so that
#' the numbers line up with malariasimulation outputs processed with postie.
#'
#' The model severe rate is first scaled by \code{severe_scaler} and then by
#' the effect of treatment, \code{(ft * (1 - treatment_scaler) + (1 - ft)) /
#' (baseline_treatment * (1 - treatment_scaler) + (1 - baseline_treatment))}.
#' A fraction \code{ft_sev} of severe cases are hospitalised, with case
#' fatality \code{hosp_sev_cfr}; the rest die with \code{community_sev_cfr}.
#'
#' The case fatality ratios (6.5\% in hospital from Reyburn et al., 60\% in
#' the community from Lubell et al.) and the hospitalised share of 0.8 are
#' from Griffin et al. 2016 Lancet Infect Dis 16:465-72, supplement section
#' 3.4.4. The \code{treatment_scaler} of 0.42 is a postie assumption based on
#' Mousa et al. 2020 PLoS Med, and the 10\% \code{baseline_treatment} is
#' hard-coded in postie (the historical coverage in the Griffin et al. 2016
#' supplement).
#'
#' @param severe Numeric rate of severe disease from the model, e.g. the
#'   \code{severe} column of \code{\link{severe_rates}}. The outputs are in the
#'   same units.
#' @param ft Proportion of clinical cases treated in the simulation.
#' @param treatment_scaler Reduction in severe disease from treating a
#'   clinical case. Default = 0.42
#' @param ft_sev Proportion of severe cases that are hospitalised.
#'   Default = 0.8
#' @param hosp_sev_cfr Case fatality ratio of hospitalised severe cases.
#'   Default = 0.065
#' @param community_sev_cfr Case fatality ratio of severe cases in the
#'   community. Default = 0.6
#' @param baseline_treatment Treatment coverage in the data the severe model
#'   was fitted to. Default = 0.1 (fixed in postie)
#' @param severe_scaler Multiplier applied to the model severe rate.
#'   Default = 1/0.8 (fixed in postie)
#'
#' @return A data frame with columns \code{severe_model} (the input),
#'   \code{severe}, \code{severe_hospital}, \code{severe_community},
#'   \code{mortality_hospital}, \code{mortality_community} and
#'   \code{mortality}.
#' @export
#' @examples
#' # severe rate of 0.01 per person per year with 40% of cases treated
#' severe_outcomes(0.01, ft = 0.4)
severe_outcomes <- function(severe, ft,
                            treatment_scaler = 0.42,
                            ft_sev = 0.8,
                            hosp_sev_cfr = 0.065,
                            community_sev_cfr = 0.6,
                            baseline_treatment = 0.1,
                            severe_scaler = 1 / 0.8) {

  if (any(treatment_scaler > 1 | treatment_scaler < 0)) {
    stop("treatment_scaler must be between 0 and 1")
  }
  if (any(ft > 1 | ft < 0) || any(ft_sev > 1 | ft_sev < 0)) {
    stop("ft and ft_sev must be between 0 and 1")
  }

  ts <- 1 - treatment_scaler
  baseline_scaler <- baseline_treatment * ts + (1 - baseline_treatment)

  sev <- severe * severe_scaler * (ft * ts + (1 - ft)) / baseline_scaler
  hosp <- sev * ft_sev
  comm <- sev * (1 - ft_sev)

  data.frame(
    severe_model = severe,
    severe = sev,
    severe_hospital = hosp,
    severe_community = comm,
    mortality_hospital = hosp * hosp_sev_cfr,
    mortality_community = comm * community_sev_cfr,
    mortality = hosp * hosp_sev_cfr + comm * community_sev_cfr
  )
}
