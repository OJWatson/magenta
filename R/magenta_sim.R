# The object pipeline() returns, and saving and restoring a run.

new_magenta_sim <- function(output, params, schedule, meta) {
  structure(output, params = params, schedule = schedule, meta = meta, class = "magenta_sim")
}

new_magenta_params <- function(x) {
  structure(x, class = c("magenta_params", "list"))
}

#' @export
print.magenta_params <- function(x, ...) {
  cat("<magenta_params>\n")
  scalars <- vapply(x, function(v) is.atomic(v) && length(v) == 1, logical(1))
  for (nm in names(x)[scalars]) {
    cat("  ", nm, ": ", format(x[[nm]]), "\n", sep = "")
  }
  for (nm in names(x)[!scalars]) {
    cat("  ", nm, ": <", class(x[[nm]])[1], ">\n", sep = "")
  }
  invisible(x)
}

#' Simulation results
#'
#' \code{\link{pipeline}} returns what the run kept, in the shape it has
#' always had: with \code{update_save = FALSE} the final slot alone (by
#' default \code{list(Loggers = ...)}, or the full model state or the human
#' population if the saver's \code{final} asks for them); with
#' \code{update_save = TRUE} a list with one element per period, the last
#' being that final slot. So \code{out$Loggers}, \code{out[[i]]} and
#' \code{length(out)} work as they did before 1.4.0.
#'
#' The object has class \code{magenta_sim} and carries the run's inputs as
#' attributes: \code{params} (the \code{magenta_params} the run was created
#' with, every parameter list included), \code{schedule} (the
#' \code{\link{build_schedule}} data frame the run followed, with the drug
#' used in each period filled in) and \code{meta} (the call, seed, package
#' version, saver, timings and the state of the random number generator at
#' the end). \code{sim_params()}, \code{sim_schedule()} and
#' \code{sim_meta()} read them; \code{attr(out, "meta")} works as it did
#' before 1.4.0.
#'
#' Printing a \code{magenta_sim} shows a short header rather than the whole
#' list. \code{unclass(out)} prints the list and its attributes, and
#' \code{sim_output(out)} is the bare list with no attributes at all, which
#' is what to compare when checking that two runs gave the same result.
#'
#' \code{as.data.frame()} gives one row per period holding every scalar in
#' that period's output: with \code{save_loggers()} these are the prevalence
#' by state, incidence, time and mosquito and treatment counts; with
#' \code{save_summary()} the prevalence and treatment failure measures.
#' \code{summary()} prints the run's settings and the distribution of those
#' columns.
#'
#' @param x,object A \code{magenta_sim}.
#' @param ... Unused.
#' @return \code{sim_output()} the result as a plain list; \code{sim_params()}
#'   the \code{magenta_params}; \code{sim_schedule()} the schedule;
#'   \code{sim_meta()} the metadata list.
#' @name magenta_sim
#' @examples
#' \donttest{
#' out <- pipeline(
#'   EIR = 2, years = 1, N = 100, seed = 1L,
#'   housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
#' )
#' out
#' names(out)
#' sim_meta(out)$seed
#' sim_schedule(out)
#' }
NULL

#' @rdname magenta_sim
#' @export
sim_output <- function(x) {
  nms <- names(x)
  attributes(x) <- NULL
  names(x) <- nms
  x
}

#' @rdname magenta_sim
#' @export
sim_params <- function(x) {
  attr(x, "params", exact = TRUE)
}

#' @rdname magenta_sim
#' @export
sim_schedule <- function(x) {
  attr(x, "schedule", exact = TRUE)
}

#' @rdname magenta_sim
#' @export
sim_meta <- function(x) {
  attr(x, "meta", exact = TRUE)
}

#' @rdname magenta_sim
#' @export
print.magenta_sim <- function(x, ...) {
  p <- sim_params(x)
  m <- sim_meta(x)
  n_out <- if (isTRUE(m$update_save)) length(x) else 1L
  cat("<magenta_sim>\n")
  cat("  call:     ", abbreviate_call(m$call), "\n", sep = "")
  cat("  N: ", p$N, "   years: ", p$years, "   EIR: ", p$EIR, "\n", sep = "")
  cat("  output:   ", m$saver, ", ", n_out, " period(s), final = ", m$final, "\n", sep = "")
  cat("  seed:     ", m$seed, "   version: ", as.character(m$version), "\n", sep = "")
  if (!is.null(m$elapsed)) {
    cat("  elapsed:  ", format(round(m$elapsed, 1)), "\n", sep = "")
  }
  invisible(x)
}

abbreviate_call <- function(call, width = 70) {
  txt <- paste(deparse(call, width.cutoff = 500), collapse = " ")
  txt <- gsub("\\s+", " ", txt)
  if (nchar(txt) > width) {
    txt <- paste0(substr(txt, 1, width - 3), "...")
  }
  txt
}

#' @rdname magenta_sim
#' @export
summary.magenta_sim <- function(object, ...) {
  structure(
    list(sim = object, periods = as.data.frame(object)),
    class = "summary.magenta_sim"
  )
}

#' @export
print.summary.magenta_sim <- function(x, ...) {
  print(x$sim)
  df <- x$periods
  numeric_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], c("period", "year"))
  if (nrow(df) && length(numeric_cols)) {
    cat("\nPer-period output (", nrow(df), " period(s)):\n", sep = "")
    print(summary(df[numeric_cols]))
  } else {
    cat("\nNo scalar per-period output to summarise.\n")
  }
  invisible(x)
}

#' @rdname magenta_sim
#' @export
as.data.frame.magenta_sim <- function(x, ...) {
  out <- sim_output(x)
  periods <- if (isTRUE(sim_meta(x)$update_save)) out else list(out)
  rows <- lapply(periods, scalar_fields)
  cols <- unique(unlist(lapply(rows, names)))
  df <- data.frame(period = seq_along(periods))
  for (col in cols) {
    df[[col]] <- vapply(rows, function(r) {
      v <- r[[col]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
  }
  sched <- sim_schedule(x)
  if (!is.null(sched) && nrow(sched) == nrow(df)) {
    df <- cbind(
      df[, "period", drop = FALSE],
      day_end = sched$day_end, year = sched$year,
      df[, setdiff(names(df), "period"), drop = FALSE]
    )
  }
  df
}

# Every length-one number in a period's output, looking inside the loggers
# and their Mosquitoes and Treatments sub-lists.
scalar_fields <- function(x) {
  if (!is.list(x)) {
    return(list())
  }
  if (!is.null(x$Loggers)) {
    x <- x$Loggers
  }
  out <- list()
  for (nm in names(x)) {
    v <- x[[nm]]
    if (is.atomic(v) && length(v) == 1 && is.numeric(v)) {
      out[[nm]] <- v
    } else if (is.list(v) && nm %in% c("Mosquitoes", "Treatments")) {
      for (nm2 in names(v)) {
        if (is.atomic(v[[nm2]]) && length(v[[nm2]]) == 1) {
          out[[nm2]] <- v[[nm2]]
        }
      }
    }
  }
  out
}

#' Save a simulation to disk and restore it
#'
#' \code{sim_save} writes the full model state at the end of a run, together
#' with the parameters, schedule and metadata of the run, as one versioned
#' object. \code{sim_restore} reads it back and checks it. The restored
#' object, or the path, can be passed to \code{\link{pipeline}} as
#' \code{saved_state_path} to continue the simulation.
#'
#' The run must have kept its final state: use a saver with
#' \code{final = "full"} (or \code{\link{save_full}}), or
#' \code{set_up_only = TRUE}.
#'
#' @param sim A \code{magenta_sim} from \code{\link{pipeline}}, or the state
#'   returned by \code{pipeline(set_up_only = TRUE)}.
#' @param path File to write to or read from (an \code{.rds} file).
#' @return \code{sim_save} returns \code{path} invisibly. \code{sim_restore}
#'   returns a \code{magenta_saved_state}: the C++ state lists, the
#'   \code{saved_inputs} needed to continue, and the \code{params},
#'   \code{schedule} and \code{meta} of the run that produced it, with a
#'   \code{format_version}.
#' @export
#' @examples
#' \dontrun{
#' out <- pipeline(EIR = 5, N = 1000, years = 2, output = save_loggers(final = "full"))
#' path <- tempfile(fileext = ".rds")
#' sim_save(out, path)
#' saved <- sim_restore(path)
#' out2 <- pipeline(EIR = 5, N = 1000, years = 1, saved_state_path = path)
#' }
sim_save <- function(sim, path) {
  check_scalar_character(path, "path")
  if (inherits(sim, "magenta_sim")) {
    state <- final_full_state(sim)
    extra <- list(params = sim_params(sim), schedule = sim_schedule(sim), meta = sim_meta(sim))
    version <- sim_meta(sim)$version
  } else if (is_full_state(sim)) {
    state <- sim
    extra <- list(params = NULL, schedule = NULL, meta = NULL)
    version <- utils::packageVersion("magenta")
  } else {
    stop("sim must be a magenta_sim from pipeline() or a full model state", call. = FALSE)
  }
  extra$meta$Ptr <- NULL
  obj <- new_saved_state(
    c(unclass(state), extra),
    format_version = 1L, version = as.character(version)
  )
  saveRDS(obj, path)
  invisible(path)
}

#' @rdname sim_save
#' @export
sim_restore <- function(path) {
  check_scalar_character(path, "path")
  if (!file.exists(path)) {
    stop("no file at ", path, call. = FALSE)
  }
  as_saved_state(readRDS(path))
}

check_scalar_character <- function(x, what) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(what, " must be a single string", call. = FALSE)
  }
  invisible(x)
}

is_full_state <- function(x) {
  is.list(x) && all(cpp_contract$saved_state %in% names(x))
}

# The final full state of a run, or an error saying how to get one.
final_full_state <- function(sim) {
  out <- sim_output(sim)
  final <- if (isTRUE(sim_meta(sim)$update_save)) out[[length(out)]] else out
  if (!is_full_state(final)) {
    stop(
      "this run did not keep its final state; run pipeline() with a saver",
      " whose final = \"full\" (for example output = save_loggers(final = \"full\"))",
      call. = FALSE
    )
  }
  final
}

new_saved_state <- function(x, format_version, version) {
  x$format_version <- format_version
  x$version <- version
  structure(x, class = c("magenta_saved_state", "list"))
}

# Check a saved state read from disk or passed directly, accepting both the
# versioned object sim_save() writes and the plain list older full saves
# were.
as_saved_state <- function(x) {
  if (!is_full_state(x)) {
    stop(
      "not a saved magenta state: expected the elements ",
      paste(cpp_contract$saved_state, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(x$saved_inputs)) {
    stop(
      "Saved state missing saved_inputs; please re-run the simulation with a",
      " saver whose final = \"full\"",
      call. = FALSE
    )
  }
  if (!inherits(x, "magenta_saved_state")) {
    x <- new_saved_state(x, format_version = 0L, version = NA_character_)
  }
  if (x$format_version > 1L) {
    stop(
      "saved state format version ", x$format_version,
      " is newer than this version of magenta understands",
      call. = FALSE
    )
  }
  x
}

#' @export
print.magenta_saved_state <- function(x, ...) {
  cat("<magenta_saved_state> format ", x$format_version, ", magenta ",
    x$version, "\n",
    sep = ""
  )
  cat("  N: ", length(x$population_List$Infection_States),
    "   day: ", x$parameters_List$g_current_time, "\n",
    sep = ""
  )
  if (!is.null(x$meta$call)) {
    cat("  from:     ", abbreviate_call(x$meta$call), "\n", sep = "")
  }
  invisible(x)
}
