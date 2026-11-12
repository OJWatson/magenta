#' Creates age brackets
#'
#' Makes age brackets given, number of brackets, max age, and spacing
#'
#' @param max_age Maximum age in years
#' @param num_age_brackets Number of age brackets
#' @param geometric_age_brackets Boolean for geometric brackets
#' @keywords internal
age_brackets <- function(max_age = 100,
                         num_age_brackets = 20,
                         geometric_age_brackets = TRUE) {
  if (geometric_age_brackets) {
    ## Create the geometric age brackets
    ratio <- (max_age / 0.1)^(1 / num_age_brackets)
    age_vector <- 0.1 * ratio**(1:num_age_brackets)
    age_vector[1] <- 0
  } else {
    age_vector <- seq(0, max_age, num_age_brackets)
  }
  return(age_vector)
}

#' Spatial matrix grab
#'
#' Grabs spatial incidence and mosquitoFOI matrix from database
#'
#' @param country Country string
#' @param admin Admin string
#' @param year_range Year range
#'
#' @export
#' @examples
#' spl <- spl_grab("Democratic Republic of the Congo", "Kinshasa", 2010:2012)
#' dim(spl$incidence)

spl_grab <- function(country, admin, year_range) {
  importations <- magenta::importations

  # grab importations
  imp_nms <- colnames(importations$`2000`$incidence)
  ads <- strsplit(imp_nms, "|", fixed = TRUE) %>%
    lapply(function(x) x[2]) %>%
    unlist()
  cc <- strsplit(imp_nms, "|", fixed = TRUE) %>%
    lapply(function(x) x[1]) %>%
    unlist()
  imp_match <- which(cc == unique(cc)[match_clean(country, unique(cc))] &
    ads == unique(ads)[match_clean(admin, unique(ads))])

  message(
    "Requested: ", admin, ", ", country,
    "\nReturned: ", ads[imp_match], ", ", cc[imp_match]
  )

  incidence <- lapply(importations, function(x) x$incidence[imp_match, ]) %>%
    unlist() %>%
    matrix(ncol = length(imp_nms), byrow = TRUE)

  mosquitoFOI <- lapply(importations, function(x) x$mosquitoFOI[imp_match, ]) %>%
    unlist() %>%
    matrix(ncol = length(imp_nms), byrow = TRUE)

  # grab year bits
  end_year <- tail(year_range, 1)
  start_year <- year_range[1]
  years <- 2000:2017

  # alter the end
  if (end_year > 2017) {
    extra <- length(tail(years, 1):end_year) - 1
    incidence <- rbind(incidence, t(replicate(extra, incidence[18, ])))
    mosquitoFOI <- rbind(mosquitoFOI, t(replicate(extra, mosquitoFOI[18, ])))
  } else if (end_year <= 2017 && end_year >= 2000) {
    incidence <- incidence[1:(end_year - 1999), , drop = FALSE]
    mosquitoFOI <- mosquitoFOI[1:(end_year - 1999), , drop = FALSE]
  } else {
    incidence <- NULL
    mosquitoFOI <- NULL
  }

  # alter the start
  if (start_year < 2000) {
    incidence <- rbind(
      matrix(0, ncol = 593, nrow = length(start_year:1999)),
      incidence
    )
    mosquitoFOI <- rbind(
      matrix(0, ncol = 593, nrow = length(start_year:1999)),
      mosquitoFOI
    )
  } else {
    incidence <- incidence[which(years == start_year):dim(incidence)[1], , drop = FALSE]
    mosquitoFOI <- mosquitoFOI[which(years == start_year):dim(mosquitoFOI)[1], , drop = FALSE]
  }


  return(list(
    "incidence" = incidence,
    "mosquitoFOI" = mosquitoFOI
  ))
}


#' Intervention grab
#'
#' @details Grabs ITN, IRS  ft from database
#'
#' N.B. The admin region used for ft can be checked
#'
#' @param country Country string
#' @param admin Admin string
#' @param year_range Year range
#' @param final_itn_cov Final ITN coverage
#' @param final_irs_cov Final IRS coverage
#' @param final_ft Final treatment coverage
#'
#' @export
#' @examples
#' ints <- intervention_grab(
#'   "Democratic Republic of the Congo", "Kinshasa", year_range = 2010:2015
#' )
#' ints$itn_cov

intervention_grab <- function(country, admin, year_range,
                              final_itn_cov = NULL,
                              final_irs_cov = NULL,
                              final_ft = NULL) {
  itn_2000_2015 <- magenta::itn_2000_2015
  irs_2000_2015 <- magenta::irs_2000_2015
  admin_units_seasonal <- magenta::admin_units_seasonal

  # incorporate hitoric interventions
  ad_mat <- itn_2000_2015$admin[match_clean(admin, itn_2000_2015$admin)]
  c_mat <- itn_2000_2015$country[match_clean(country, itn_2000_2015$country)]

  itn_rows <- itn_2000_2015$admin == ad_mat & itn_2000_2015$country == c_mat
  irs_rows <- irs_2000_2015$admin == ad_mat & irs_2000_2015$country == c_mat
  itn_cov <- itn_2000_2015$value[which(itn_rows)]
  irs_cov <- irs_2000_2015$value[which(irs_rows)]

  # admatch for ft
  adm <- admin_match(admin, country, quiet = TRUE)
  ft <- admin_units_seasonal$ft[adm]
  ft <- rep(ft, 16)

  # check it's the right length
  if (length(irs_cov) != 16) {
    irs_cov <- irs_2000_2015$value[which(irs_2000_2015$country == country)]
  }
  if (length(irs_cov) != 16) {
    irs_cov <- rep(0, 16)
  }

  # grab year bits
  end_year <- tail(year_range, 1)
  start_year <- year_range[1]
  years <- 2000:2015

  # check the final coverages
  if (is.null(final_itn_cov)) final_itn_cov <- tail(itn_cov, 1)
  if (is.null(final_irs_cov)) final_irs_cov <- tail(irs_cov, 1)
  if (is.null(final_ft)) final_ft <- tail(ft, 1)

  # alter the end
  if (end_year > 2015) {
    itn_cov <- c(itn_cov, final_itn_cov)
    irs_cov <- c(irs_cov, final_irs_cov)
    ft <- c(ft, final_ft)

    itn_cov <- approx(x = c(years, end_year), y = itn_cov, xout = 2000:end_year)$y
    irs_cov <- approx(x = c(years, end_year), y = irs_cov, xout = 2000:end_year)$y
    ft <- approx(x = c(years, end_year), y = ft, xout = 2000:end_year)$y
  } else if (end_year <= 2015 && end_year > 2000) {
    itn_cov <- itn_cov[1:(end_year - 1999)]
    irs_cov <- irs_cov[1:(end_year - 1999)]
    ft <- ft[1:(end_year - 1999)]
  } else {
    itn_cov <- rep(0, length(year_range))
    irs_cov <- rep(0, length(year_range))
    ft <- rep(0, length(year_range))
    return(list(
      "itn_cov" = itn_cov,
      "irs_cov" = irs_cov,
      "ft" = ft
    ))
  }

  if (year_range[1] < 2000) {
    itn_cov <- c(rep(0, length(start_year:1999)), itn_cov)
    irs_cov <- c(rep(0, length(start_year:1999)), irs_cov)
    ft <- c(rep(0, length(start_year:1999)), ft)
  } else {
    itn_cov <- itn_cov[which(years == start_year):length(itn_cov)]
    irs_cov <- irs_cov[which(years == start_year):length(irs_cov)]
    ft <- ft[which(years == start_year):length(ft)]
  }


  return(list(
    "itn_cov" = itn_cov,
    "irs_cov" = irs_cov,
    "ft" = ft
  ))
}

#' Mosquito death and feeding rates for a run
#'
#' The daily mosquito death rate \code{mu} and feeding rate \code{fv} under
#' the run's bed net and spraying coverage, which the simulation reads day by
#' day. Both come from the intervention block of the deterministic model and
#' depend on nothing but the coverage and the parameters, so they are
#' evaluated directly by \code{\link{mu_fv_from_coverage}}; with
#' \code{full = TRUE} the deterministic model itself is run instead and its
#' whole output returned (see \code{\link{run_deterministic_model}}).
#'
#' @param eqInit Output of \code{\link{equilibrium_init_create}}.
#' @param ft Treatment seeking, one value per year.
#' @param itn_cov,irs_cov Bed net and spraying coverage, one value per year.
#' @param years Length of the simulation in years.
#' @param int_times Days at which the yearly values start; see
#'   \code{\link{mu_fv_from_coverage}}.
#' @param full Run the deterministic model and return all of its output.
#'   Default = FALSE
#' @param odin_model Unused; kept for compatibility.
#' @return With \code{full = FALSE}, a data frame with a row per day and
#'   columns \code{mu} and \code{fv}; otherwise the list of model outputs
#'   from \code{\link{run_deterministic_model}}.
#' @keywords internal
mu_fv_create <- function(eqInit,
                         ft,
                         itn_cov,
                         irs_cov,
                         years,
                         int_times = NULL,
                         full = FALSE,
                         odin_model = NULL) {
  if (all(irs_cov == 0) && all(itn_cov == 0)) {
    out <- data.frame(
      "mu" = rep(0.132, length(seq_len(round(years * 365)) + 1)),
      "fv" = rep(0.333, length(seq_len(round(years * 365)) + 1))
    )
  } else if (full) {
    out <- run_odin_model(eqInit, ft, itn_cov, irs_cov, years, int_times)
  } else {
    out <- mu_fv_from_coverage(eqInit, ft, itn_cov, irs_cov, years, int_times)
  }
  out
}

# Run the deterministic model over the simulation and return every output
# (transformed variables, including t, mu and fv), as the pipeline did
# before the closed form.
run_odin_model <- function(eqInit, ft, itn_cov, irs_cov, years, int_times = NULL) {
  message("\nRunning Deterministic Model for ", years, " years")

  # the yearly vectors come padded from build_schedule(); the first value
  # is repeated so that the model starts from the year-one coverage
  n_years <- ceiling(years)
  first_then <- function(x) {
    x <- pad_annual(x, n_years, deparse(substitute(x)))
    c(x[1], x)
  }
  eqInit$ft_vector <- first_then(ft)
  eqInit$itn_vector <- first_then(itn_cov)
  eqInit$irs_vector <- first_then(irs_cov)

  # set up the intervetnion times
  if (is.null(int_times)) {
    int_times <- c(-25, seq(0, by = 365, length.out = n_years))
  }
  eqInit$t_vector <- int_times
  eqInit$ITN_IRS_on <- 0

  # split pop accordingly
  eqInit$pop_split <- rep(1 / eqInit$num_int, eqInit$num_int)

  extensions <- which(unlist(lapply(eqInit, function(x) {
    (length(dim(x)) == 3)
  })))
  mat <- matrix(0, eqInit$na, eqInit$nh)
  for (e in names(extensions)) {
    if (grepl("_I", e)) {
      eqInit[[e]] <- array(eqInit[[e]][, , 1], c(eqInit$na, eqInit$nh, eqInit$num_int))
    } else {
      eqInit[[e]] <- vapply(rep(1 / eqInit$num_int, eqInit$num_int),
        FUN = function(x) {
          x * eqInit[[e]][, , 1]
        }, mat
      )
    }
  }

  # build model
  state <- eqInit[names(eqInit) %in% odin_itn_irs$private_fields$user]

  # weird catach for when init_ICM values fall below 5e-32 that must catch
  state$init_ICM
  model <- odin_itn_irs$new(user = state, use_dde = TRUE)

  # create model and simulate
  tt <- seq(1, round(years * 365), 1)
  mod_run <- model$run(t = tt, n_history = 1000, step_size_max = 10)
  model$transform_variables(mod_run)
}

#' Run the deterministic model
#'
#' Runs the deterministic transmission model that \code{magenta} carries
#' (\code{inst/extdata/odin_itn_irs.R}, an early copy of the Imperial
#' malaria model with ITN and IRS) from its equilibrium for a given setting,
#' and returns every model output through time. The simulation itself only
#' takes its mosquito death and feeding rates from this model (see
#' \code{\link{mu_fv_create}}); the full run is for checking the
#' individual-based simulation against the deterministic one, for example
#' that the prevalence and incidence a large population settles at for an
#' \code{EIR} match.
#'
#' @inheritParams equilibrium_init_create
#' @param years Length of the run in years.
#' @param itn_cov,irs_cov Bed net and spraying coverage, one value per year.
#'   Default = 0
#' @param int_times Days at which the yearly values start; see
#'   \code{\link{mu_fv_from_coverage}}.
#' @return List of the model's outputs, one vector (or array) per output with
#'   a value per day: \code{t}, the compartment totals (\code{Sout},
#'   \code{Tout}, ...), prevalence and incidence, and the mosquito rates
#'   \code{mu} and \code{fv}, among others.
#' @export
#' @examples
#' \donttest{
#' # the deterministic model and the simulation at the same EIR, no resistance
#' det <- run_deterministic_model(EIR = 10, ft = 0.4, years = 2, itn_cov = 0.3)
#' # daily prevalence and clinical incidence from the deterministic model
#' det$prev[c(1, 365, 730)]
#' det$inc[c(1, 365, 730)]
#'
#' # the simulation's incidence over the same two years, per 365-day period
#' sim <- pipeline(
#'   EIR = 10, ft = 0.4, years = 2, N = 10000, itn_cov = 0.3,
#'   update_save = TRUE, update_length = 365, output = save_loggers()
#' )
#' sapply(sim$output, function(x) x$Loggers$Incidence)
#' }
run_deterministic_model <- function(EIR, ft, years,
                                    itn_cov = 0, irs_cov = 0,
                                    age_vector = age_brackets(100, 40, TRUE),
                                    het_brackets = 5,
                                    model_param_list = model_param_list_create(),
                                    country = NULL, admin = NULL,
                                    int_times = NULL, quiet = FALSE) {
  eq_init <- equilibrium_init_create(
    age_vector = age_vector, het_brackets = het_brackets, ft = ft[1],
    EIR = EIR, model_param_list = model_param_list, country = country,
    admin = admin, quiet = quiet
  )
  run_odin_model(eq_init, ft, itn_cov, irs_cov, years, int_times)
}


drug_list_update <- function(drug_list, year, tf) {
  if (drug_list$resistance_flag) {
    if (!drug_list$mft_flag) {
      # if sequential cycling
      if (drug_list$sequential_cycling > 0) {
        if (!is.na(tf) && tf > drug_list$sequential_cycling &&
          year > drug_list$next_temporal_cycle) {
          drug_list$drug_choice <- drug_list$drug_choice + 1
          drug_list$next_temporal_cycle <- year + drug_list$sequential_update
          if (drug_list$drug_choice == (drug_list$number_of_drugs)) {
            drug_list$drug_choice <- 0
          }

          message("Sequential cycling: switched to drug ", drug_list$drug_choice)
        }
      }

      # if temporal cycling
      if (drug_list$temporal_cycling > 0) {
        if (drug_list$next_temporal_cycle == year) {
          drug_list$drug_choice <- drug_list$drug_choice + 1
          if (drug_list$drug_choice == (drug_list$number_of_drugs)) {
            drug_list$drug_choice <- 0
          }
          drug_list$next_temporal_cycle <- drug_list$next_temporal_cycle +
            drug_list$temporal_cycling
        }
      }
    }
  }

  return(drug_list)
}

# convert allele frequencies
pop_alf <- function(nums, nl, weighted = FALSE) {
  if (is.raw(unlist(nums))) {
    if (weighted) {
      res <- colSums(do.call(rbind, lapply(nums, function(x) {
        colMeans(matrix(as.numeric(x), ncol = nl))
      }))) / length(nums)
    } else {
      res <- matrix(as.numeric(do.call(rbind, nums)), ncol = nl)
      res <- colSums(res) / nrow(res)
    }
  } else {
    res <- rep(NA, nl)
  }
  return(unlist(res))
}

# convert to lineages frequencies
lineages <- function(nums, nl) {
  if (is.raw(unlist(nums))) {
    res <- table(
      factor(
        apply(matrix(as.numeric(do.call(rbind, nums)), ncol = nl), 1, bitsToInt),
        levels = 0:((2^nl) - 1)
      )
    )
  } else {
    res <- rep(NA, 2^nl)
  }
  return(unlist(res))
}
