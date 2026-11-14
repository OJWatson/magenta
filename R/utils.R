# The magrittr pipe is used internally and no longer re-exported.
#' @importFrom magrittr %>%
#' @noRd
NULL

#' Progress messages for the main loop
#'
#' Ticks the progress bar, or on a cluster (where a bar cannot redraw) prints
#' the percentage done every two percent, as \code{pipeline()} has always
#' done.
#'
#' @param housekeeping_list The housekeeping list, for \code{quiet_print} and
#'   \code{cluster}.
#' @param res The results list, for the number of periods.
#' @param progress_bar A \code{progress::progress_bar}.
#' @param i The period just started.
#' @param initial Whether this is the call before the loop.
#' @param p_print The next percentage to print.
#' @return The next percentage to print, to pass back in.
#' @keywords internal
progress_logging <- function(housekeeping_list, res, progress_bar,
                             i = 1, initial = FALSE, p_print = 2) {
  if (initial) {
    if (housekeeping_list$quiet_print) {
      if (housekeeping_list$cluster) {
        p_print <- 2
        message(paste0("||:|"), appendLF = FALSE)
      } else {
        progress_bar$tick()
      }
    }
  } else {
    p_perc <- (i / length(res)) * 100

    if (p_perc > p_print) {
      if (housekeeping_list$quiet_print) {
        if (housekeeping_list$cluster) {
          message(paste0(p_print, "|"), appendLF = FALSE)
          if (p_print == 100) {
            message(paste0(p_print + 2, "|:||"), appendLF = FALSE)
          }
        }
      }
      p_print <- p_print + 2
    }

    if (!housekeeping_list$cluster) {
      progress_bar$tick()
    }
  }

  return(p_print)
}

#' Bind a list of data frames by row
#'
#' @param x List of data frames with the same columns.
#' @return One data frame with the rows renumbered.
#' @keywords internal
rbind_list_base <- function(x) {
  x2 <- do.call(
    rbind.data.frame,
    c(x, stringsAsFactors = FALSE, make.row.names = FALSE)
  )
  rownames(x2) <- seq_len(dim(x2)[1])
  x2
}

#' Consecutive index ranges
#'
#' @param diff Length of each range.
#' @param end Total length.
#' @return List of integer vectors \code{1:diff}, \code{(diff + 1):(2 * diff)}
#'   and so on.
#' @keywords internal
ranges <- function(diff, end) {
  r <- list()
  for (i in 1:(end / diff)) {
    r[[i]] <- (1 + ((i - 1) * diff)):(diff * i)
  }
  return(r)
}

#' Function to convert vector of bools representing up to a 32 bit
#' number into the corresponding integer
#'
#' @param x Vector of raw/logicals
#' @param big_endian Is big_endian. Default = FALSE
#' @keywords internal
bitsToInt <- function(x, big_endian = FALSE) {
  if (big_endian) {
    sum(2^(which(rev(as.logical(x))) - 1))
  } else {
    sum(2^(which((as.logical(x))) - 1))
  }
}

#' Apply a check unless the value is NULL
#'
#' @param x The value.
#' @param func A check such as \code{\link{assert_scalar_character}}.
#' @param name Name used in the error message.
#' @return What \code{func} returns, or \code{NULL}.
#' @keywords internal
assert_null_and_func <- function(x, func, name = deparse(substitute(x))) {
  if (!is.null(x)) {
    func(x, name)
  }
}

#' Check for a single string
#'
#' @param x The value.
#' @param name Name used in the error message.
#' @return \code{x}, invisibly.
#' @keywords internal
assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar character", name))
  }
  invisible(x)
}

#' Check for a single TRUE or FALSE
#'
#' @param x The value.
#' @param name Name used in the error message.
#' @return \code{x}, invisibly.
#' @keywords internal
assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar logical", name))
  }
  invisible(x)
}

#' Check for a single number
#'
#' @param x The value.
#' @param name Name used in the error message.
#' @return \code{x}, invisibly.
#' @keywords internal
assert_scalar_numeric <- function(x, name = deparse(substitute(x))) {
  if (!(is.numeric(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar numeric", name))
  }
  invisible(x)
}
