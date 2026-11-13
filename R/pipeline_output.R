# Savers: what pipeline() extracts from the simulation after each period.
#
# A saver is a function(state, period, params) returning whatever should be
# kept for that period. `state` is the internal simulation state (its `sim`
# element holds the C++ pointer and the loggers), `period` is the schedule
# row for the period with its index `i` and the number of periods `n`, and
# `params` is the magenta_params object recorded with the run. The four
# constructors below wrap such a function with a name and a choice of what
# the final slot of the output holds.

#' Construct a saver
#'
#' A saver is a \code{function(state, period, params)} returning whatever
#' should be kept for a period, wrapped with a name and a choice of what the
#' final slot of the output holds. See \code{\link{savers}} for the four
#' constructors built on this.
#'
#' @param fun The function. \code{state} is the internal simulation state
#'   (its \code{sim} element holds the C++ pointer and the loggers),
#'   \code{period} the schedule row for the period with its index \code{i}
#'   and the number of periods \code{n}, and \code{params} the
#'   \code{magenta_params} recorded with the run.
#' @param name Name shown by \code{print()} and kept in \code{meta$saver}.
#' @param final What the last element of the output holds.
#' @return A \code{magenta_saver}.
#' @keywords internal
new_saver <- function(fun, name, final = c("loggers", "full", "humans")) {
  final <- match.arg(final, c("loggers", "full", "humans"))
  if (!is.function(fun) || length(formals(fun)) < 3) {
    stop("a saver must be a function(state, period, params)", call. = FALSE)
  }
  structure(list(fun = fun, name = name, final = final), class = "magenta_saver")
}

#' Output savers for \code{pipeline()}
#'
#' A saver decides what \code{\link{pipeline}} keeps after each update period
#' and what its final slot holds. Pass one as the \code{output} argument.
#'
#' \describe{
#'   \item{\code{save_loggers()}}{The loggers (prevalence by state, incidence,
#'     immunity, mosquito counts, treatment counts) after each period.}
#'   \item{\code{save_summary()}}{A genetic and epidemiological summary of the
#'     human population after each period: either the complexity of infection
#'     tables from \code{COI_df_create()} (\code{genetics = "coi"}), the
#'     population allele frequencies (\code{"allele_freqs"}, optionally with
#'     the frequency of each lineage) or the raw strain data frame
#'     (\code{"df"}), followed by the treatment and prevalence summaries used
#'     in resistance work.}
#'   \item{\code{save_humans()}}{The human population and its strains after
#'     each period.}
#'   \item{\code{save_full()}}{The full model state after each period, as
#'     needed to restart the simulation.}
#' }
#'
#' Any \code{function(state, period, params)} can be passed instead; it
#' receives the internal state (\code{state$sim$Loggers} and
#' \code{state$sim$Ptr} are the parts of interest), the schedule row for the
#' period (with its index \code{i} and the number of periods \code{n}) and
#' the \code{magenta_params} recorded with the run, and returns what to keep.
#'
#' @param final What the last element of the output holds: the loggers
#'   (\code{"loggers"}, the default), the full model state (\code{"full"})
#'   so that the run can be continued with \code{\link{sim_save}}, or the
#'   human population only (\code{"humans"}).
#' @param sample_size Number of individuals sampled at the end of each period.
#'   \code{Inf} samples everyone. A vector of sample sizes samples at each.
#' @param sample_states Infection states included in the sample. Default
#'   \code{0:5} is every state; \code{1:4} would be infected individuals only.
#' @param sample_reps Number of sample replicates. Default = 1.
#' @param age_breaks Age breaks in years used when summarising the population.
#' @param mean_only Whether to return only the mean when summarising the
#'   population COI, COU and related quantities. Default = TRUE.
#' @param genetics Which genetic summary to keep: \code{"coi"},
#'   \code{"allele_freqs"} or \code{"df"}. See Description.
#' @param lineages With \code{genetics = "allele_freqs"}, also keep the
#'   frequency of each lineage (each distinct barcode). Default = FALSE.
#' @param barcodes_every,barcodes_last With \code{genetics = "coi"}, the
#'   tabled barcodes are kept every \code{barcodes_every} periods during the
#'   last \code{barcodes_last} periods.
#' @return A \code{magenta_saver}.
#' @name savers
#' @examples
#' \dontrun{
#' out <- pipeline(
#'   EIR = 5, N = 1000, years = 3, update_save = TRUE, update_length = 90,
#'   output = save_summary(genetics = "allele_freqs", final = "full")
#' )
#' }
NULL

#' @rdname savers
#' @export
save_loggers <- function(final = c("loggers", "full", "humans")) {
  final <- match.arg(final)
  new_saver(
    function(state, period, params) state$sim$Loggers,
    name = "loggers", final = final
  )
}

#' @rdname savers
#' @export
save_summary <- function(sample_size = Inf, sample_states = 0:5, sample_reps = 1,
                         age_breaks = c(-0.001, 5, 15, 100.1), mean_only = TRUE,
                         genetics = c("coi", "allele_freqs", "df"),
                         lineages = FALSE, barcodes_every = 12,
                         barcodes_last = 180,
                         final = c("loggers", "full", "humans")) {
  genetics <- match.arg(genetics)
  final <- match.arg(final)
  check_number(sample_size, "sample_size", min = 1, n = NULL)
  check_number(sample_states, "sample_states", min = 0, max = 5, n = NULL, integer = TRUE)
  check_number(sample_reps, "sample_reps", min = 1, integer = TRUE)
  check_number(age_breaks, "age_breaks", n = NULL)
  check_flag(mean_only, "mean_only")
  check_flag(lineages, "lineages")
  check_number(barcodes_every, "barcodes_every", min = 1, integer = TRUE)
  check_number(barcodes_last, "barcodes_last", min = 0, integer = TRUE)

  fun <- function(state, period, params) {
    sim_out <- state$sim
    ibd <- params$barcode$barcode_type

    # sample the population's strains
    df <- pop_strains_df(sim_out$Ptr,
      sample_size = if (length(sample_size) > 1) 0 else sample_size * sample_reps,
      sample_states = sample_states, ibd = ibd
    )

    res <- switch(genetics,
      coi = COI_df_create(df,
        barcodes = period$i %% barcodes_every == 0 && period$i >= (period$n - barcodes_last),
        ibd = ibd, breaks = age_breaks, n = sample_size, reps = sample_reps,
        mean_only = mean_only
      ),
      allele_freqs = {
        has_strains <- unlist(lapply(df$barcode_states, length)) > 0
        out <- list(af = pop_alf(df$nums[has_strains], params$num_loci))
        if (lineages) {
          out$lineage <- lineages(df$nums[has_strains], params$num_loci)
        }
        out
      },
      df = list(df = df)
    )
    c(res, treatment_and_prevalence_summary(sim_out, params$core))
  }
  new_saver(fun, name = paste0("summary (", genetics, ")"), final = final)
}

#' @rdname savers
#' @export
save_humans <- function(final = c("loggers", "full", "humans")) {
  final <- match.arg(final)
  new_saver(
    function(state, period, params) {
      collect_output(sim_get_state(state$sim$Ptr), human_only = TRUE)
    },
    name = "humans", final = final
  )
}

#' @rdname savers
#' @export
save_full <- function(final = c("full", "loggers", "humans")) {
  final <- match.arg(final)
  new_saver(
    function(state, period, params) prepare_state_save(state),
    name = "full", final = final
  )
}

#' @param x,... A saver, and arguments ignored by \code{print()}.
#' @rdname savers
#' @export
print.magenta_saver <- function(x, ...) {
  cat("<magenta_saver> ", x$name, ", final = ", x$final, "\n", sep = "")
  invisible(x)
}

#' Treatment outcome and prevalence summaries for a period
#'
#' The summaries appended to every \code{\link{save_summary}} result, in the
#' order the output has always had: treatment counts and failure, the daily
#' population EIR, mutations, and the PCR, patent and microscopy (2 to 10 year
#' olds) prevalence.
#'
#' @param sim_out The C++ handle and loggers after an update.
#' @param mpl Model parameter list, for the detection parameters.
#' @return A named list of the summaries.
#' @keywords internal
treatment_and_prevalence_summary <- function(sim_out, mpl) {
  loggers <- sim_out$Loggers
  res <- list()
  res$succesfull_treatments <- loggers$Treatments$Successful_Treatments
  res$unsuccesful_treatments_lpf <- loggers$Treatments$Unsuccesful_Treatments_LPF
  res$not_treated <- loggers$Treatments$Not_Treated
  res$treatment_failure <- res$unsuccesful_treatments_lpf /
    (res$unsuccesful_treatments_lpf + res$succesfull_treatments)
  res$daily_pop_eir <- loggers$Treatments$daily_infectious_bite_counters /
    loggers$Log_Counter * 365 / length(loggers$Ages)
  res$mutations <- loggers$daily_mutations_per_loci

  res$treatment_failure[is.na(res$treatment_failure)] <- 0
  res$overall_treatment_failure <- period_treatment_failure(sim_out)

  # prevalence measures
  two_ten <- loggers$Ages < 3650 & loggers$Ages > 720
  dt_10 <- sum(loggers$InfectionStates %in% c(1, 4) & two_ten) / sum(two_ten)
  micro_det <- q_fun(
    mpl$d1, loggers$ID, mpl$ID0, mpl$kD,
    sapply(loggers$Ages, fd, mpl$fD0, mpl$aD, mpl$gammaD)
  )
  a_10 <- sum(loggers$InfectionStates %in% 2 & two_ten) / sum(two_ten) *
    mean(micro_det[two_ten])
  res$pcr_prev <- sum(unlist(loggers[c("D", "A", "U", "T")]))
  res$dat <- sum(unlist(loggers[c("D", "A", "T")]))
  res$micro_2_10 <- dt_10 + a_10
  res
}

#' Treatment failure in a period
#'
#' The proportion of treatments in a period that failed (late parasitological
#' failure), from the loggers every update returns; \code{0} when nobody was
#' treated. The sequential cycling rule in \code{\link{drug_list_update}}
#' reads this.
#'
#' @param sim_out The C++ handle and loggers after an update.
#' @return A single number between 0 and 1.
#' @keywords internal
period_treatment_failure <- function(sim_out) {
  t <- sim_out$Loggers$Treatments
  tf <- sum(t$Unsuccesful_Treatments_LPF) /
    sum(t$Unsuccesful_Treatments_LPF + t$Successful_Treatments)
  if (is.na(tf)) 0 else tf
}

#' Wrap an old update_save_func as a saver
#'
#' The deprecated \code{update_save_func} argument took a
#' \code{function(res, i, sim.out, mpl, num_loci)} returning the whole
#' results list. This calls it that way and keeps its slot for the period.
#'
#' @param fun The old-style function.
#' @param final What the last element of the output holds.
#' @return A \code{magenta_saver}.
#' @keywords internal
legacy_saver <- function(fun, final = "loggers") {
  res <- list()
  new_saver(
    function(state, period, params) {
      length(res) <<- period$n
      res <<- fun(res, period$i, state$sim, params$core, params$num_loci)
      res[[period$i]]
    },
    name = "update_save_func", final = final
  )
}

#' The saver for a pipeline() call
#'
#' Turns the \code{output} argument of \code{\link{pipeline}}, or the
#' deprecated flags it replaces, into a saver. Giving both is an error; the
#' flags alone give a warning naming the saver that replaces them.
#'
#' @param output The \code{output} argument: \code{NULL}, a
#'   \code{magenta_saver} or a function.
#' @param flags Named list of the values of the deprecated arguments.
#' @param used Named logical vector saying which of them the caller set.
#' @return A \code{magenta_saver}.
#' @keywords internal
resolve_saver <- function(output, flags, used) {
  if (!is.null(output) && !inherits(output, "magenta_saver")) {
    if (is.function(output)) {
      output <- new_saver(output, name = "user function")
    } else {
      stop(
        "output must be a saver (save_loggers(), save_summary(), save_humans(),",
        " save_full()) or a function(state, period, params)",
        call. = FALSE
      )
    }
  }
  if (!any(used)) {
    return(if (is.null(output)) save_loggers() else output)
  }
  if (!is.null(output)) {
    stop(
      "output was given together with the deprecated argument(s) ",
      paste(names(used)[used], collapse = ", "),
      "; set what they asked for through the saver instead",
      call. = FALSE
    )
  }
  mapped <- saver_from_flags(flags)
  warning(
    "pipeline(): the argument(s) ", paste(names(used)[used], collapse = ", "),
    " are deprecated; use output = ", mapped$replacement, " instead",
    call. = FALSE
  )
  mapped$saver
}

#' The saver that a combination of the old output flags means
#'
#' @param f Named list of the values of the deprecated arguments, as
#'   \code{\link{resolve_saver}} collects them.
#' @return List with the \code{saver} and the \code{replacement} call to
#'   show in the deprecation warning.
#' @keywords internal
saver_from_flags <- function(f) {
  final <- if (isTRUE(f$human_only_full_save)) {
    "humans"
  } else if (isTRUE(f$full_save)) {
    "full"
  } else {
    "loggers"
  }
  final_arg <- if (final == "loggers") "" else sprintf('final = "%s"', final)
  args <- function(...) {
    a <- c(..., final_arg)
    paste(a[nzchar(a)], collapse = ", ")
  }

  if (!is.null(f$update_save_func)) {
    saver <- legacy_saver(f$update_save_func, final = final)
    return(list(saver = saver, replacement = "function(state, period, params)"))
  }
  if (!isTRUE(f$human_update_save)) {
    return(list(
      saver = save_loggers(final = final),
      replacement = sprintf("save_loggers(%s)", args())
    ))
  }
  if (!isTRUE(f$summary_saves_only)) {
    if (isTRUE(f$full_update_save)) {
      return(list(
        saver = save_full(final = final),
        replacement = sprintf("save_full(%s)", args())
      ))
    }
    return(list(
      saver = save_humans(final = final),
      replacement = sprintf("save_humans(%s)", args())
    ))
  }
  genetics <- if (!isTRUE(f$genetics_df_without_summarising)) {
    "coi"
  } else if (isTRUE(f$only_allele_freqs)) {
    "allele_freqs"
  } else {
    "df"
  }
  saver <- save_summary(
    sample_size = f$sample_size, sample_states = f$sample_states,
    sample_reps = f$sample_reps, age_breaks = f$age_breaks,
    mean_only = f$mean_only, genetics = genetics,
    lineages = isTRUE(f$save_lineages), final = final
  )
  shown <- c(
    sprintf("sample_size = %s", deparse1(f$sample_size)),
    sprintf("sample_states = %s", deparse1(f$sample_states)),
    if (f$sample_reps != 1) sprintf("sample_reps = %s", f$sample_reps),
    if (!isTRUE(f$mean_only)) "mean_only = FALSE",
    sprintf('genetics = "%s"', genetics),
    if (isTRUE(f$save_lineages)) "lineages = TRUE"
  )
  list(saver = saver, replacement = sprintf("save_summary(%s)", args(shown)))
}

# The object pipeline() returns, and saving and restoring a run.

#' Construct the object pipeline() returns
#'
#' @param output What the run kept: the per-period results, or the final slot
#'   alone.
#' @param params The \code{magenta_params} of the run.
#' @param schedule The schedule the run followed.
#' @param meta The metadata list.
#' @return \code{output} with the three as attributes and class
#'   \code{magenta_sim}; see \code{\link{magenta_sim}}.
#' @keywords internal
new_magenta_sim <- function(output, params, schedule, meta) {
  structure(output, params = params, schedule = schedule, meta = meta, class = "magenta_sim")
}

#' Construct the record of a run's inputs
#'
#' @param x Named list of the inputs: the scalar settings and every parameter
#'   list.
#' @return \code{x} with class \code{magenta_params}.
#' @keywords internal
new_magenta_params <- function(x) {
  structure(x, class = c("magenta_params", "list"))
}

#' @rdname magenta_sim
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
#' @param x,object A \code{magenta_sim} (for \code{print.magenta_params()},
#'   a \code{magenta_params}; for \code{print.summary.magenta_sim()}, the
#'   summary).
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

#' A call on one line, cut to a width
#'
#' @param call The call.
#' @param width Characters to keep.
#' @return A single string, ending in \code{...} if cut.
#' @keywords internal
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

#' @rdname magenta_sim
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

#' Every length-one number in a period's output
#'
#' Looks inside the loggers (a period kept as \code{list(Loggers = ...)} or
#' as the loggers themselves) and their \code{Mosquitoes} and
#' \code{Treatments} sub-lists, which is what \code{as.data.frame()} turns
#' into columns.
#'
#' @param x One period's output.
#' @return A named list of single numbers, empty if \code{x} is not a list.
#' @keywords internal
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

#' Check for a single string
#'
#' @param x The value.
#' @param what Name used in the error message.
#' @return \code{x}, invisibly.
#' @keywords internal
check_scalar_character <- function(x, what) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(what, " must be a single string", call. = FALSE)
  }
  invisible(x)
}

#' Is this a full model state?
#'
#' @param x Any object.
#' @return \code{TRUE} if \code{x} is a list with the five lists the C++
#'   returns, named in \code{cpp_contract$saved_state}.
#' @keywords internal
is_full_state <- function(x) {
  is.list(x) && all(cpp_contract$saved_state %in% names(x))
}

#' The final full state of a run
#'
#' @param sim A \code{magenta_sim}.
#' @return The full model state in its final slot, or an error saying how to
#'   get one.
#' @keywords internal
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

#' Stamp a saved state with its format and version
#'
#' @param x The saved state as a list.
#' @param format_version Format of the object: \code{0} for the plain full
#'   saves of earlier versions, \code{1} for what \code{\link{sim_save}}
#'   writes.
#' @param version The magenta version that wrote it.
#' @return \code{x} with the two fields and class
#'   \code{magenta_saved_state}.
#' @keywords internal
new_saved_state <- function(x, format_version, version) {
  x$format_version <- format_version
  x$version <- version
  structure(x, class = c("magenta_saved_state", "list"))
}

#' Check a saved state read from disk or passed directly
#'
#' Accepts both the versioned object \code{\link{sim_save}} writes and the
#' plain list older full saves were, upgrading the latter.
#'
#' @param x The object.
#' @return A \code{magenta_saved_state}, or an error saying what is wrong
#'   with it.
#' @keywords internal
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

#' @param x,... A \code{magenta_saved_state}, and arguments ignored by
#'   \code{print()}.
#' @rdname sim_save
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
