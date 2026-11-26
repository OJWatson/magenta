#------------------------------------------------
#' Seasonal parameters from a site file
#'
#' Turns the rainfall Fourier coefficients of a malariasimulation site file
#' (\code{site$seasonality$seasonality_parameters}, one unit) into magenta's
#' seasonal parameters. magenta's seasonal curve is
#' \code{max((a0 + sum(a_i cos + b_i sin)) / theta_c, 0.001)}, so
#' \code{a0 = g0}, \code{a_i = g_i}, \code{b_i = h_i}, and \code{theta_c} is the
#' yearly mean of the curve floored at 0.001, which is malariasimulation's
#' normalising constant.
#'
#' @param site Site file for one unit, e.g. from \code{site::subset_site()}
#'
#' @return Named numeric vector (a0, a1, a2, a3, b1, b2, b3, theta_c) for the
#'   \code{seasonality} argument of \code{\link{pipeline}}
#' @export
site_seasonality <- function(site) {
  sp <- site$seasonality$seasonality_parameters
  if (is.null(sp) || nrow(sp) != 1) {
    stop("site needs exactly one row of seasonality_parameters (one unit)")
  }
  g <- unlist(sp[1, c("g0", "g1", "g2", "g3", "h1", "h2", "h3")])
  w <- 2 * pi * (1:365) / 365
  r <- g[["g0"]] + g[["g1"]] * cos(w) + g[["g2"]] * cos(2 * w) + g[["g3"]] * cos(3 * w) +
    g[["h1"]] * sin(w) + g[["h2"]] * sin(2 * w) + g[["h3"]] * sin(3 * w)
  c(a0 = g[["g0"]], a1 = g[["g1"]], a2 = g[["g2"]], a3 = g[["g3"]],
    b1 = g[["h1"]], b2 = g[["h2"]], b3 = g[["h3"]], theta_c = mean(pmax(r, 0.001)))
}

#------------------------------------------------
#' Treatment seeking from a site file's effective treatment coverage
#'
#' Site files give \code{tx_cov}, the effective treatment coverage, and
#' malariasimulation treats that fraction of clinical cases with a fully
#' effective drug. magenta instead treats a fraction \code{ft} with the drugs in
#' \code{drug_list}, which can fail. \code{site_ft()} sets \code{ft} so that the
#' drugs' baseline (no resistance) efficacy gives the site's effective coverage:
#' \code{ft = tx_cov / efficacy}, capped at 1, where efficacy is the mean of the
#' drugs' wild type \code{lpf} weighted by \code{partner_drug_ratios} (as
#' \code{\link{pipeline}} uses for the deterministic model). Resistance then
#' lowers the effective coverage below the site's.
#'
#' The site file also splits \code{tx_cov} into ACT and non-ACT
#' (\code{prop_act}), but both are given at full efficacy in malariasimulation,
#' so the split does not change the cleared fraction. It is not used here: all
#' treatment is with \code{drug_list}.
#'
#' @param tx_cov Numeric vector of effective treatment coverage
#' @param drug_list Drug list, see \code{\link{drug_list_create}}
#'
#' @return Numeric vector of ft
#' @export
site_ft <- function(tx_cov, drug_list = drug_list_create()) {
  lpfs <- vapply(drug_list$drugs, function(x) x$lpf[1], numeric(1))
  efficacy <- stats::weighted.mean(lpfs, drug_list$partner_drug_ratios)
  pmin(tx_cov / efficacy, 1)
}

#------------------------------------------------
#' Pipeline arguments from a malariasimulation site file
#'
#' Turns one unit of an Imperial malariasimulation site file (e.g. a foresite
#' calibrated, scaled site, reduced to one unit with \code{site::subset_site()})
#' into the arguments for \code{\link{pipeline}}, so magenta runs on the site's
#' calibrated EIR without its own fitting.
#'
#' \itemize{
#'   \item EIR: the site's \code{pf} EIR (adult EIR, as magenta's \code{EIR}).
#'   \item Seasonality: \code{\link{site_seasonality}}.
#'   \item Burn-in: \code{burnin} years before the first site year. With
#'     \code{burnin_interventions = "none"} the burn-in has no treatment, ITNs or
#'     IRS (so the starting equilibrium is at ft = 0); with \code{"first"} it
#'     holds the first year's values.
#'   \item Yearly interventions: \code{itn_cov} is the site's ITN use,
#'     \code{irs_cov} the site's IRS coverage times \code{irs_scale} (0.75, the
#'     adjustment malariasimulation runs apply), and \code{ft} from
#'     \code{\link{site_ft}} with \code{drug_list}. Years missing from the site
#'     are 0.
#'   \item Nets: \code{d_ITN0}, \code{r_ITN0} and \code{r_ITN1} are
#'     \code{dn0}, \code{rn0} and \code{rnm} from \code{site::net_efficacy}
#'     (\code{net_type}) at the unit's pyrethroid resistance, averaged over the
#'     site years weighted by ITN use (magenta has one net efficacy for the whole
#'     run). The ITN half-life is \code{gamman * log(2)} (malariasimulation's
#'     decay is \code{exp(-t / gamman)}), and nets are aged with the "constant"
#'     clock over a life of the site's retention half-life (see
#'     \code{itn_decay_clock} in \code{\link{model_param_list_create}}). This is
#'     the mean efficacy of nets in use when nets are topped up every year to
#'     hold use steady, as the site's yearly distributions do. Without a
#'     retention half-life the net life is 3 years.
#'   \item Vectors: \code{Q0}, \code{bites_Bed}, \code{bites_Indoors} and
#'     \code{mu0} are the species means (weighted by \code{prop}) of \code{Q0},
#'     \code{phi_bednets}, \code{phi_indoors} and \code{mum}.
#'   \item \code{eta}: given (default 1 / (22 years)); magenta's population is
#'     exponential so the site's demography is not used.
#' }
#'
#' Drugs, resistance, non-malarial fevers, barcodes and the other
#' \code{pipeline} arguments are left to the user: add them to the returned
#' list. Use the same \code{drug_list} here and in \code{pipeline}.
#'
#' @param site Site file for one unit, e.g. from \code{site::subset_site()}
#' @param years Calendar years to run after the burn-in. Default = all years
#'   with treatment coverage in the site file
#' @param burnin Years of burn-in before \code{years}. Default = 20
#' @param burnin_interventions "none" (default) or "first", see Details
#' @param drug_list Drug list used for \code{ft}, see \code{\link{site_ft}}
#' @param irs_scale Multiplier on the site's IRS coverage. Default = 0.75
#' @param net_type Net type in \code{site::net_efficacy}. Default =
#'   "pyrethroid_only"
#' @param eta Human death rate per day. Default = 1 / (22 * 365)
#' @param update_length Days between updates in \code{pipeline}, which must
#'   divide 365. Default = 365
#'
#' @return Named list of arguments for \code{pipeline}, with attributes
#'   \code{year} (calendar year of each simulated year, burn-in included) and
#'   \code{tx_cov} (the site's treatment coverage in each simulated year, e.g.
#'   for \code{\link{severe_outcomes}})
#' @export
#' @examples
#' \dontrun{
#' s <- readRDS("calibrated_scaled_site.rds")
#' unit <- site::subset_site(s, data.frame(name_1 = "Busia", urban_rural = "urban"))
#' args <- site_to_pipeline(unit)
#' out <- do.call(pipeline, c(args, list(N = 10000)))
#' }
site_to_pipeline <- function(site,
                             years = NULL,
                             burnin = 20,
                             burnin_interventions = c("none", "first"),
                             drug_list = drug_list_create(),
                             irs_scale = 0.75,
                             net_type = "pyrethroid_only",
                             eta = 1 / (22 * 365),
                             update_length = 365) {

  if (!requireNamespace("site", quietly = TRUE)) {
    stop("site_to_pipeline needs the site package (mrc-ide/site)")
  }
  burnin_interventions <- match.arg(burnin_interventions)
  int <- site$interventions

  # EIR
  eir <- site$eir[site$eir$sp == "pf", ]
  if (nrow(eir) != 1) stop("site needs exactly one pf EIR (one unit)")

  # yearly inputs
  tx <- int$treatment$implementation
  if (is.null(years)) years <- sort(unique(tx$year))
  by_year <- function(d, col, fun = max) {
    v <- vapply(years, function(y) {
      x <- d[[col]][d$year == y]
      if (length(x) == 0) 0 else fun(x)
    }, numeric(1))
    v
  }
  tx_cov <- by_year(tx, "tx_cov")
  itn_use <- by_year(int$itn$use, "itn_use")
  irs_cov <- by_year(int$irs$implementation, "irs_cov")

  pad <- function(x) {
    c(rep(if (burnin_interventions == "none") 0 else x[1], burnin), x)
  }
  tx_cov <- pad(tx_cov)

  # nets at the unit's resistance, weighted by use
  res <- site$vectors$pyrethroid_resistance
  res <- res[res$year %in% years, ]
  w <- itn_use[match(res$year, years)]
  if (sum(w) == 0) w <- rep(1, length(w))
  resistance <- round(stats::weighted.mean(res$pyrethroid_resistance, w), 2)
  ne <- site::net_efficacy
  ne <- ne[ne$net_type == net_type &
             round(ne$pyrethroid_resistance, 2) == resistance, ]
  if (nrow(ne) != 1) stop("no net_efficacy row for ", net_type, " at resistance ", resistance)

  retention <- int$itn$retention_half_life

  # species weighted vector parameters
  v <- site$vectors$vector_species
  if (abs(sum(v$prop) - 1) > 1e-6) stop("vector_species prop must sum to 1 (one unit)")
  wsum <- function(col) sum(v$prop * v[[col]])

  args <- list(
    EIR = eir$eir,
    ft = site_ft(tx_cov, drug_list),
    itn_cov = pad(itn_use),
    irs_cov = irs_scale * pad(irs_cov),
    years = burnin + length(years),
    update_save = TRUE,
    update_length = update_length,
    seasonality = site_seasonality(site),
    drug_list = drug_list,
    d_ITN0 = ne$dn0,
    r_ITN0 = ne$rn0,
    r_ITN1 = ne$rnm,
    itn_half_life = ne$gamman * log(2),
    ITN_interval = if (is.null(retention)) 3 * 365 else retention,
    itn_decay_clock = "constant",
    Q0 = wsum("Q0"),
    bites_Bed = wsum("phi_bednets"),
    bites_Indoors = wsum("phi_indoors"),
    mu0 = wsum("mum"),
    eta = eta
  )
  attr(args, "year") <- c(years[1] - rev(seq_len(burnin)), years)
  attr(args, "tx_cov") <- tx_cov
  args
}
