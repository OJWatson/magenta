# Helpers for the fixed-seed characterisation tests.
#
# The characterisation tests pin the current output of pipeline() for a handful
# of small configurations so that refactors of the R layer can be checked as
# behaviour-preserving. They are not a statement that the pinned values are
# correct. A change that is meant to alter results updates the snapshots
# deliberately (testthat::snapshot_accept()) and gets a NEWS entry.

# Reduce an arbitrary pipeline() output to something small and readable enough
# to live in a snapshot file. Short vectors are kept verbatim; long numeric
# vectors are reduced to a length, a sum and a few order statistics; external
# pointers and functions are replaced by a marker.
snapshot_summary <- function(x, max_keep = 12L) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "externalptr")) {
    return("<externalptr>")
  }
  if (is.function(x) || is.language(x)) {
    return(paste("<", class(x)[1], ">"))
  }
  if (is.data.frame(x)) {
    return(lapply(as.list(x), snapshot_summary, max_keep = max_keep))
  }
  if (is.table(x)) {
    return(list(names = names(x), counts = as.integer(x)))
  }
  if (is.list(x)) {
    out <- lapply(x, snapshot_summary, max_keep = max_keep)
    if (is.null(names(x)) && length(out) > max_keep) {
      return(list(
        n_elements = length(out),
        first = out[seq_len(5)],
        last = out[length(out)]
      ))
    }
    return(out)
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.raw(x)) {
    x <- as.integer(x)
  }
  if (is.matrix(x)) {
    return(list(
      dim = dim(x),
      colsums = round(colSums(x), 6),
      rowsums = snapshot_summary(round(rowSums(x), 6), max_keep)
    ))
  }
  if (is.character(x) || is.logical(x)) {
    if (length(x) <= max_keep) {
      return(x)
    }
    return(list(length = length(x), table = snapshot_summary(table(x))))
  }
  if (is.numeric(x)) {
    if (length(x) <= max_keep) {
      return(round(x, 6))
    }
    return(list(
      length = length(x),
      n_na = sum(is.na(x)),
      sum = round(sum(x, na.rm = TRUE), 6),
      quantiles = round(
        unname(stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)),
        6
      )
    ))
  }
  paste("<", class(x)[1], ">")
}

# Quiet housekeeping so the tests do not print progress bars or seeds.
quiet_housekeeping <- function() {
  housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)
}

# Snapshots record floating point sums from the C++ model. They are only
# compared on Linux, where the CI runners share a compiler and libm.
skip_characterisation_if_unsupported <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_os(c("windows", "mac"))
}
