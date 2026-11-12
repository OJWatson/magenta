#' Pipe operator
#'
#' See \code{\link[magrittr:pipe]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

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

#' @noRd
rbind_list_base <- function(x) {
  x2 <- do.call(
    rbind.data.frame,
    c(x, stringsAsFactors = FALSE, make.row.names = FALSE)
  )
  rownames(x2) <- seq_len(dim(x2)[1])
  x2
}

#' @noRd
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

# Convert integer to binary for a given n
binary <- function(x, n = 24) {
  i <- 0
  string <- numeric(n)
  while (x > 0) {
    string[n - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1
  }
  string
}

# function to put all of one infividuals' parms into one place
person_make <- function(o, n) {
  if (length(names(o)) == 6) {
    df <- lapply(o[1:3], function(x) x[n]) %>% as.data.frame()
    l <- lapply(o[4:6], function(x) x[n])
  } else {
    df <- lapply(o$population_List, function(x) x[n]) %>% as.data.frame()
    oths <- which(lapply(o$populations_event_and_strains_List, class) %>% unlist() != "list")
    df <- cbind(df, lapply(o$populations_event_and_strains_List[oths], function(x) x[n])) %>%
      as.data.frame()
    l <- lapply(o$populations_event_and_strains_List[-unique(c(oths))], function(x) x[n])
  }
  return(list("vars" = df, "l" = l))
}

#' @noRd
func_call <- function(args, call) {
  paste0(
    as.character(call)[1],
    "(",
    paste0(mapply(paste, names(args), args, sep = "="), collapse = ", "),
    ")"
  )
}

## Assertions

#' @noRd
assert_null_and_func <- function(x, func, name = deparse(substitute(x))) {
  if (!is.null(x)) {
    func(x, name)
  }
}

#' @noRd
assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar character", name))
  }
  invisible(x)
}

#' @noRd
assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar logical", name))
  }
  invisible(x)
}

#' @noRd
assert_scalar_numeric <- function(x, name = deparse(substitute(x))) {
  if (!(is.numeric(x) && length(x) == 1L && !is.na(x))) {
    stop(sprintf("'%s' must be a scalar numeric", name))
  }
  invisible(x)
}
