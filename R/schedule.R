#' Build the per-period simulation schedule
#'
#' \code{build_schedule} turns the yearly inputs to \code{\link{pipeline}} into
#' one data frame with a row per update period. The main loop reads a row at a
#' time and never has to recycle, pad or index the yearly vectors itself.
#'
#' @section Periods:
#' The simulation covers \code{round(years * 365)} days split into
#' \code{round(years * 365 / update_length)} periods of \code{update_length}
#' days. The last period is shortened if it would run past the final day, and
#' any days beyond the last full period are not simulated, which is why a one
#' year simulation with monthly (30 day) updates runs for 360 days. A period
#' belongs to the year in which its first day falls. Setting
#' \code{update_length} to the total number of days gives a single period.
#'
#' @section Padding rule:
#' Every yearly input (\code{ft}, \code{itn_cov}, \code{irs_cov},
#' \code{resistance_flag}, \code{mutation_flag}, the rows of \code{plaf} and of
#' the two spatial matrices) is extended to \code{ceiling(years)} entries by
#' repeating its last value. A single value therefore applies to every year,
#' and a vector shorter than the simulation holds its final value to the end.
#' An input longer than \code{ceiling(years)} is an error, so that a mismatch
#' between the inputs and the length of the simulation is caught rather than
#' silently truncated.
#'
#' @param years Length of the simulation in years. A vector of calendar years
#'   (e.g. \code{2000:2015}) is taken to mean its length.
#' @param update_length Length of each period in days.
#' @param ft Treatment seeking rate, one value or one per year.
#' @param itn_cov ITN coverage, one value or one per year.
#' @param irs_cov IRS coverage, one value or one per year.
#' @param resistance_flag Whether resistance is simulated, one value or one
#'   per year.
#' @param mutation_flag Whether mutation is simulated, one value or one per
#'   year.
#' @param plaf Population level allele frequencies: a vector (one per locus)
#'   or a matrix with one row per year.
#' @param spatial_incidence Human importation rates: \code{NULL}, a vector with
#'   one value per year, or a matrix with one row per year (as returned by
#'   \code{\link{spl_grab}}).
#' @param spatial_mosquitoFOI Mosquito importation rates, in the same form as
#'   \code{spatial_incidence}.
#' @param drug_choice The drug in use at the start. The column is updated by
#'   the pipeline as drug cycling switches drugs, so that the schedule returned
#'   with a simulation records the drug used in each period.
#' @return A data frame of class \code{magenta_schedule} with the columns
#'   \code{period}, \code{day_start}, \code{day_end}, \code{year}, \code{ft},
#'   \code{itn_cov}, \code{irs_cov}, \code{resistance_flag},
#'   \code{mutation_flag}, \code{drug_choice} and the list columns \code{plaf},
#'   \code{importation_human} and \code{importation_mosquito}. The padded
#'   yearly vectors are kept in the \code{annual} attribute for the
#'   deterministic mosquito model.
#' @export
#' @examples
#' sched <- build_schedule(
#'   years = 2, update_length = 365, ft = c(0.2, 0.6), itn_cov = 0.3,
#'   plaf = rep(0.5, 4)
#' )
#' sched[, c("period", "day_start", "day_end", "year", "ft", "itn_cov")]
build_schedule <- function(years, update_length, ft = 0.4, itn_cov = 0,
                           irs_cov = 0, resistance_flag = FALSE,
                           mutation_flag = FALSE, plaf = rep(0.5, 24),
                           spatial_incidence = NULL, spatial_mosquitoFOI = NULL,
                           drug_choice = 0) {
  if (length(years) > 1) {
    years <- length(years)
  }
  check_number(years, "years", min = 1 / 365)
  check_number(update_length, "update_length", min = 1, integer = TRUE)

  n_years <- ceiling(years)
  n_days <- round(years * 365)
  n_periods <- max(1L, as.integer(round(n_days / update_length)))

  day_start <- (seq_len(n_periods) - 1L) * update_length + 1L
  day_end <- pmin(seq_len(n_periods) * update_length, n_days)
  year <- pmin(as.integer(floor(day_start / 365)) + 1L, n_years)

  annual <- list(
    ft = check_number(pad_annual(ft, n_years, "ft"), "ft", min = 0, max = 1, n = NULL),
    itn_cov = check_number(pad_annual(itn_cov, n_years, "itn_cov"), "itn_cov",
      min = 0, max = 1, n = NULL
    ),
    irs_cov = check_number(pad_annual(irs_cov, n_years, "irs_cov"), "irs_cov",
      min = 0, max = 1, n = NULL
    ),
    resistance_flag = check_flag(
      pad_annual(resistance_flag, n_years, "resistance_flag"), "resistance_flag",
      n = NULL
    ),
    mutation_flag = check_flag(
      pad_annual(mutation_flag, n_years, "mutation_flag"), "mutation_flag",
      n = NULL
    ),
    plaf = pad_annual_rows(plaf, n_years, "plaf"),
    spatial_incidence = pad_annual_rows(spatial_incidence, n_years, "spatial_incidence"),
    spatial_mosquitoFOI = pad_annual_rows(
      spatial_mosquitoFOI, n_years, "spatial_mosquitoFOI"
    )
  )
  check_number(annual$plaf, "plaf", min = 0, max = 1, n = NULL)
  # the importation matrices from spl_grab() hold NAs for admin units without
  # data; spl_create() sums them with na.rm = TRUE
  if (!identical(dim(annual$spatial_incidence), dim(annual$spatial_mosquitoFOI))) {
    stop(
      "spatial_incidence and spatial_mosquitoFOI must have the same dimensions",
      call. = FALSE
    )
  }

  row_or_null <- function(m, i) {
    if (is.null(m)) {
      return(NULL)
    }
    m[i, ]
  }

  schedule <- data.frame(
    period = seq_len(n_periods),
    day_start = day_start,
    day_end = day_end,
    year = year,
    ft = annual$ft[year],
    itn_cov = annual$itn_cov[year],
    irs_cov = annual$irs_cov[year],
    resistance_flag = annual$resistance_flag[year],
    mutation_flag = annual$mutation_flag[year],
    drug_choice = rep(drug_choice, n_periods)
  )
  schedule$plaf <- lapply(year, function(y) annual$plaf[y, ])
  schedule$importation_human <- lapply(year, function(y) {
    row_or_null(annual$spatial_incidence, y)
  })
  schedule$importation_mosquito <- lapply(year, function(y) {
    row_or_null(annual$spatial_mosquitoFOI, y)
  })

  attr(schedule, "annual") <- annual
  class(schedule) <- c("magenta_schedule", "data.frame")
  schedule
}

# Extend a yearly vector to n_years entries by repeating its last value, or
# stop if it is longer than the simulation.
pad_annual <- function(x, n_years, what) {
  if (length(x) == 0) {
    stop(sprintf("%s must have at least one value", what), call. = FALSE)
  }
  if (length(x) > n_years) {
    stop(sprintf(
      "%s has %d values but the simulation covers %d year(s)",
      what, length(x), n_years
    ), call. = FALSE)
  }
  c(x, rep(x[length(x)], n_years - length(x)))
}

# The same rule for inputs with one row per year. A vector is one row that
# applies to every year, except for the spatial inputs, where a plain vector
# is one value per year as spl_grab() and the old pipeline() accepted.
pad_annual_rows <- function(x, n_years, what) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric", what), call. = FALSE)
  }
  if (!is.matrix(x)) {
    if (what == "plaf") {
      x <- matrix(x, nrow = 1)
    } else {
      x <- matrix(pad_annual(x, n_years, what), ncol = 1)
    }
  }
  if (nrow(x) > n_years) {
    stop(sprintf(
      "%s has %d rows but the simulation covers %d year(s)",
      what, nrow(x), n_years
    ), call. = FALSE)
  }
  extra <- n_years - nrow(x)
  if (extra > 0) {
    x <- rbind(x, matrix(rep(x[nrow(x), ], each = extra), nrow = extra))
  }
  unname(x)
}

#' @export
print.magenta_schedule <- function(x, ...) {
  days <- if (is.null(x$day_end)) NA else max(x$day_end)
  cat("<magenta_schedule> ", nrow(x), " period(s) over ", days, " days\n", sep = "")
  shown <- c(
    "period", "day_start", "day_end", "year", "ft", "itn_cov", "irs_cov",
    "resistance_flag", "mutation_flag", "drug_choice"
  )
  df <- x
  class(df) <- "data.frame"
  print(df[intersect(shown, names(df))], ...)
  invisible(x)
}
