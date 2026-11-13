#' Build the per-period simulation schedule
#'
#' \code{build_schedule} turns the yearly inputs to \code{\link{pipeline}} into
#' one data frame with a row per update period. The main loop reads a row at a
#' time and never has to recycle, pad or index the yearly vectors itself, and
#' the yearly value of an input is the value in the first row of that year.
#'
#' @section Periods:
#' The simulation covers \code{round(years * 365)} days split into
#' \code{round(years * 365 / update_length)} periods of \code{update_length}
#' days. The last period is shortened if it would run past the final day, and
#' any days beyond the last full period are not simulated, which is why a one
#' year simulation with monthly (30 day) updates runs for 360 days. A period
#' belongs to the year in which its first day falls. Every simulated year must
#' contain the start of a period, because the yearly inputs are switched when a
#' period starts a new year: an \code{update_length} that leaves a simulated
#' year with no period starting in it (one longer than a year, or between 244
#' and 364 days in some runs) is an error.
#'
#' @section Padding rule:
#' Every yearly input (\code{ft}, \code{itn_cov}, \code{irs_cov},
#' \code{resistance_flag} and \code{mutation_flag}) is extended to
#' \code{ceiling(years)} entries by repeating its last value. A single value
#' therefore applies to every year, and a vector shorter than the simulation
#' holds its final value to the end. An input longer than \code{ceiling(years)}
#' is an error, so that a mismatch between the inputs and the length of the
#' simulation is caught rather than silently truncated. The allele frequencies
#' and importation rates, which are a vector per year rather than a value, are
#' padded by the same rule and held beside the schedule (see
#' \code{\link{pipeline}}); they are not columns of it.
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
#' @param drug_choice The drug in use at the start. The column is updated by
#'   the pipeline as drug cycling switches drugs, so that the schedule returned
#'   with a simulation records the drug used in each period.
#' @return A data frame of class \code{magenta_schedule} with one row per
#'   period and the columns \code{period}, \code{day_start}, \code{day_end},
#'   \code{year}, \code{ft}, \code{itn_cov}, \code{irs_cov},
#'   \code{resistance_flag}, \code{mutation_flag} and \code{drug_choice}.
#' @export
#' @examples
#' build_schedule(years = 2, update_length = 365, ft = c(0.2, 0.6), itn_cov = 0.3)
#' # monthly periods: the yearly values repeat down the rows of each year
#' head(build_schedule(years = 2, update_length = 30, ft = c(0.2, 0.6)), 14)
build_schedule <- function(years, update_length, ft = 0.4, itn_cov = 0,
                           irs_cov = 0, resistance_flag = FALSE,
                           mutation_flag = FALSE, drug_choice = 0) {
  if (length(years) > 1) {
    years <- length(years)
  }
  check_number(years, "years", min = 1 / 365)
  check_number(update_length, "update_length", min = 1, integer = TRUE)

  # whole periods only, as the package has always counted them
  n_years <- ceiling(years)
  n_days <- round(years * 365)
  n_periods <- max(1L, as.integer(round(n_days / update_length)))

  day_start <- (seq_len(n_periods) - 1L) * update_length + 1L
  day_end <- pmin(seq_len(n_periods) * update_length, n_days)
  year <- pmin(as.integer(floor(day_start / 365)) + 1L, n_years)

  # the yearly inputs are switched when a period starts a new year, so a year
  # that is simulated but has no period starting in it would never get its
  # values; say so rather than run with the previous year's
  last_simulated_year <- (max(day_end) - 1L) %/% 365L + 1L
  if (last_simulated_year > max(year)) {
    stop(sprintf(paste0(
      "update_length = %d leaves no period starting in year %d, which is still ",
      "simulated; the yearly inputs switch when a period starts a year, so use ",
      "update_length = 365 or one short enough that every year starts a period"
    ), update_length, max(year) + 1L), call. = FALSE)
  }

  # each yearly input padded to the simulation length, then read by year
  ft <- check_number(pad_annual(ft, n_years, "ft"), "ft", min = 0, max = 1, n = NULL)
  itn_cov <- check_number(pad_annual(itn_cov, n_years, "itn_cov"), "itn_cov",
    min = 0, max = 1, n = NULL
  )
  irs_cov <- check_number(pad_annual(irs_cov, n_years, "irs_cov"), "irs_cov",
    min = 0, max = 1, n = NULL
  )
  resistance_flag <- check_flag(
    pad_annual(resistance_flag, n_years, "resistance_flag"), "resistance_flag",
    n = NULL
  )
  mutation_flag <- check_flag(
    pad_annual(mutation_flag, n_years, "mutation_flag"), "mutation_flag",
    n = NULL
  )

  schedule <- data.frame(
    period = seq_len(n_periods),
    day_start = day_start,
    day_end = day_end,
    year = year,
    ft = ft[year],
    itn_cov = itn_cov[year],
    irs_cov = irs_cov[year],
    resistance_flag = resistance_flag[year],
    mutation_flag = mutation_flag[year],
    drug_choice = rep(drug_choice, n_periods)
  )
  class(schedule) <- c("magenta_schedule", "data.frame")
  schedule
}

#' The yearly inputs that are a vector per year
#'
#' The allele frequencies and the two importation rate inputs to
#' \code{\link{pipeline}} are one vector per year rather than one value, so
#' they do not fit the schedule's columns. This pads each to a matrix with
#' \code{ceiling(years)} rows by the schedule's rule (last row repeated, more
#' rows than years an error) and checks it, so that the simulation can read
#' row \code{year} when a period starts a new year.
#'
#' @param years Length of the simulation in years, or a vector of calendar
#'   years.
#' @param plaf Population level allele frequencies: a vector (one per locus,
#'   applying to every year) or a matrix with one row per year.
#' @param spatial_incidence Human importation rates: \code{NULL}, a vector with
#'   one value per year, or a matrix with one row per year (as returned by
#'   \code{\link{spl_grab}}).
#' @param spatial_mosquitoFOI Mosquito importation rates, in the same form as
#'   \code{spatial_incidence}.
#' @return List with the padded \code{plaf}, \code{spatial_incidence} and
#'   \code{spatial_mosquitoFOI} matrices, the last two \code{NULL} when not
#'   given.
#' @keywords internal
annual_inputs <- function(years, plaf, spatial_incidence = NULL,
                          spatial_mosquitoFOI = NULL) {
  if (length(years) > 1) {
    years <- length(years)
  }
  n_years <- ceiling(years)
  annual <- list(
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
  annual
}

#' Row of a yearly matrix
#'
#' @param m A matrix with one row per year, or \code{NULL}.
#' @param year The year to read.
#' @return The row as a vector, or \code{NULL} when \code{m} is.
#' @keywords internal
annual_row <- function(m, year) {
  if (is.null(m)) {
    return(NULL)
  }
  m[year, ]
}

#' Pad a yearly vector to the length of the simulation
#'
#' Extends a yearly input to \code{n_years} entries by repeating its last
#' value, the padding rule of \code{\link{build_schedule}}, or stops if the
#' input is empty or longer than the simulation.
#'
#' @param x The yearly values.
#' @param n_years Number of years the simulation covers.
#' @param what Name of the input, for the error message.
#' @return \code{x} with \code{n_years} entries.
#' @keywords internal
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

#' Pad a yearly matrix to the length of the simulation
#'
#' The padding rule of \code{\link{pad_annual}} for inputs with one row per
#' year: the last row is repeated to \code{n_years} rows, and more rows than
#' years is an error. A plain vector is one row that applies to every year,
#' except for the spatial inputs, where a plain vector is one value per year
#' as \code{\link{spl_grab}} and \code{\link{pipeline}} have always accepted.
#'
#' @param x A matrix with one row per year, a vector, or \code{NULL}.
#' @param n_years Number of years the simulation covers.
#' @param what Name of the input, for the error message.
#' @return An unnamed matrix with \code{n_years} rows, or \code{NULL} when
#'   \code{x} is.
#' @keywords internal
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

#' @param x,... A schedule, and arguments passed on to the data frame
#'   \code{print()} method.
#' @rdname build_schedule
#' @export
print.magenta_schedule <- function(x, ...) {
  days <- if (is.null(x$day_end)) "" else paste0(" over ", max(x$day_end), " days")
  cat("<magenta_schedule> ", nrow(x), " period(s)", days, "\n", sep = "")
  df <- x
  class(df) <- "data.frame"
  print(df, ...)
  invisible(x)
}
