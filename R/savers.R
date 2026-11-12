# Savers: what pipeline() extracts from the simulation after each period.
#
# A saver is a function(state, period, params) returning whatever should be
# kept for that period. `state` is the internal simulation state (its `sim`
# element holds the C++ pointer and the loggers), `period` is the schedule
# row for the period with its index `i` and the number of periods `n`, and
# `params` is the magenta_params object recorded with the run. The four
# constructors below wrap such a function with a name and a choice of what
# the final slot of the output holds.

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
      sample_states = sample_states, ibd = ibd, seed = params$seed
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
      collect_output(sim_get_state(state$sim$Ptr, params$seed), human_only = TRUE)
    },
    name = "humans", final = final
  )
}

#' @rdname savers
#' @export
save_full <- function(final = c("full", "loggers", "humans")) {
  final <- match.arg(final)
  new_saver(
    function(state, period, params) prepare_state_save(state, params$seed),
    name = "full", final = final
  )
}

#' @export
print.magenta_saver <- function(x, ...) {
  cat("<magenta_saver> ", x$name, ", final = ", x$final, "\n", sep = "")
  invisible(x)
}

# The treatment outcome and prevalence summaries appended to every
# save_summary() result, in the order the output has always had.
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
  res$overall_treatment_failure <- sum(res$unsuccesful_treatments_lpf) /
    (sum(res$unsuccesful_treatments_lpf + res$succesfull_treatments))
  if (is.na(res$overall_treatment_failure)) {
    res$overall_treatment_failure <- 0
  }

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

# Wrap a function written for the old update_save_func argument,
# function(res, i, sim.out, mpl, num_loci) returning the whole results list,
# as a saver.
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

# Turn the `output` argument of pipeline(), or the deprecated flags it
# replaces, into a saver. `flags` holds the values of the old arguments and
# `used` which of them the caller set.
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

# The old flag combinations and the saver each one means.
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
