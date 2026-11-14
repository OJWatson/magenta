#' Create housekeeping parameter list
#'
#' List for simulation housekeeping vars, e.g. quiet prints,
#'
#' @param quiet Boolean for quiet simulation. Default = TRUE
#' @param quiet_test Boolean for quiet testing statement. Default = TRUE
#' @param cluster Boolean for simulation being on cluster. Default = TRUE
#' @param clear_up Boolean for whether to clear up the memory used by the
#'   simulation. Default = TRUE
#' @return A \code{magenta_housekeeping_list}: a list with the fields the
#'   simulation reads, see \code{\link{cpp_contract}}.
#' @examples
#' hk <- housekeeping_list_create(quiet = TRUE, cluster = FALSE)
#' hk
#' validate_housekeeping_list(hk)
#' @export
housekeeping_list_create <- function(quiet = TRUE,
                                     quiet_test = TRUE,
                                     cluster = FALSE,
                                     clear_up = TRUE) {
  l <- list(
    "quiet_print" = quiet,
    "quiet_test_print" = quiet_test,
    "cluster" = cluster,
    "clear_up" = clear_up
  )

  validate_housekeeping_list(new_param_list(l, "magenta_housekeeping_list"))
}

#' @rdname housekeeping_list_create
#' @param x List to validate.
#' @export
validate_housekeeping_list <- function(x) {
  what <- "housekeeping_list"
  check_fields(x, cpp_contract$housekeeping, what)
  for (field in cpp_contract$housekeeping$fields) {
    check_flag(x[[field]], paste0(what, "$", field))
  }
  invisible(x)
}


#' Create drug list
#'
#' List for simulating drug usage for resistance/mft variables
#'
#' @param resistance_flag Boolean are we simulating resistance
#' @param absolute_fitness_cost_flag Boolean are we simulating fitness costs as
#'   absolute (i.e. impactive onward transmission chance) or relative (resistant
#'   strains have decreases chance of being onwardly transmitted in mixed
#'   infections). Default = FALSE.
#' @param number_of_resistance_loci Numeric for number of res. loci. With
#'   resistance on it must be at least the number of distinct barcode
#'   positions the \code{drugs} name, or the list is an error.
#' @param cost_of_resistance Numeric vector for fitness cost of each resistance loci.
#' @param epistatic_logic Is there compensatory relationships. i.e. what loci
#'   need to be true for resistance cost to exist. Default of NULL means that
#'   this becomes seq_len(number_of_resistance_loci), i.e. only dependent on
#'   their own loci. (TODO: Change this to be a list of length norl)
#' @param artemisinin_loci Numerics for barcode positions that confer artemisinin
#'   resistance
#' @param drugs List of drugs that are being used, with each list element being
#'   created by \code{drug_create}.
#' @param mft_flag Boolean are we doing mft
#' @param temporal_cycling Numeric for when in years a drug switch occurs
#' @param sequential_cycling Numeric for what perc. treatment failure before switch
#' @param sequential_update How long does it take in years for sequential to be implemented
#' @param number_of_drugs Numeric for number of drugs used
#' @param drug_choice What's the default drug choice to begin. Default = 0
#' @param partner_drug_ratios Numeric vector for ratio of first line drugs used
#' @return A \code{magenta_drug_list}: a list with the fields the simulation
#'   reads, see \code{\link{cpp_contract}}.
#' @examples
#' # a single, fully effective drug with no resistance (the default)
#' drug_list_create()
#'
#' # artemether-lumefantrine with resistance at six loci and a fitness cost
#' dl <- drug_list_create(
#'   resistance_flag = TRUE, number_of_resistance_loci = 6,
#'   cost_of_resistance = rep(0.99, 6), drugs = list(drug_create_al()),
#'   number_of_drugs = 1, partner_drug_ratios = 1
#' )
#' dl$resistance_loci
#' validate_drug_list(dl)
#' @export
drug_list_create <- function(resistance_flag = FALSE,
                             number_of_resistance_loci = 2,
                             artemisinin_loci = c(0),
                             cost_of_resistance = c(0.99, 0.99),
                             absolute_fitness_cost_flag = FALSE,
                             epistatic_logic = NULL,
                             number_of_drugs = 1,
                             drugs = list(drug_create_default_no_resistance()),
                             mft_flag = FALSE,
                             temporal_cycling = -1,
                             sequential_cycling = -1,
                             sequential_update = 3,
                             drug_choice = 0,
                             partner_drug_ratios = rep(1 / number_of_drugs, number_of_drugs)) {
  # TODO: More check for formatting here that you have enough prophylactic etc

  if (temporal_cycling > 0 && sequential_cycling > 0) {
    stop("Both sequential and temporal cycling can't be greater than 0")
  }

  resistance_cost_table_create <- function(number_of_resistance_loci,
                                           cost_of_resistance,
                                           epistatic_logic = NULL) {
    if (length(cost_of_resistance) > number_of_resistance_loci) {
      cost_of_resistance <- cost_of_resistance[seq_len(number_of_resistance_loci)]
    }

    if (is.null(epistatic_logic)) {
      epistatic_logic <- seq_len(number_of_resistance_loci)
    }

    max_bits <- (2^(number_of_resistance_loci)) - 1
    barcode_poss <- matrix(intToBits(0:max_bits), ncol = 32, byrow = TRUE)
    barcode_poss <- barcode_poss[, seq_len(number_of_resistance_loci), drop = FALSE]
    results <- rep(1, nrow(barcode_poss))
    count <- 1

    for (j in seq_len(nrow(barcode_poss))) {
      comp <- epistatic_logic[as.logical(barcode_poss[j, ])]
      res <- as.logical(barcode_poss[j, ])
      res[which(res)] <- comp %in% which(res)
      results[count] <- prod(cost_of_resistance[res])
      count <- count + 1
    }

    return(results)
  }

  cost_of_resistance <- resistance_cost_table_create(
    number_of_resistance_loci,
    cost_of_resistance,
    epistatic_logic
  )

  resistance_loci <- sort(unique(unlist(
    lapply(drugs, "[[", "barcode_positions")
  )))


  l <- list(
    "resistance_flag" = resistance_flag,
    "number_of_resistance_loci" = number_of_resistance_loci,
    "resistance_loci" = resistance_loci,
    "artemisinin_loci" = artemisinin_loci,
    "cost_of_resistance" = cost_of_resistance,
    "absolute_fitness_cost_flag" = absolute_fitness_cost_flag,
    "drugs" = drugs,
    "mft_flag" = mft_flag,
    "temporal_cycling" = temporal_cycling,
    "next_temporal_cycle" = temporal_cycling,
    "sequential_cycling" = sequential_cycling,
    "sequential_update" = sequential_update,
    "number_of_drugs" = number_of_drugs,
    "drug_choice" = drug_choice,
    "partner_drug_ratios" = partner_drug_ratios
  )

  validate_drug_list(new_param_list(l, "magenta_drug_list"))
}

#' @rdname drug_list_create
#' @param x List to validate.
#' @export
validate_drug_list <- function(x) {
  what <- "drug_list"
  check_fields(x, cpp_contract$drug, what)
  # resistance_flag may be one value or one per simulation year; the pipeline
  # slices it to the current year before the list reaches C++
  check_flag(x$resistance_flag, paste0(what, "$resistance_flag"), n = NULL)
  check_number(x$number_of_resistance_loci, paste0(what, "$number_of_resistance_loci"),
    min = 1, integer = TRUE
  )
  check_number(x$resistance_loci, paste0(what, "$resistance_loci"),
    min = 0, n = NULL, integer = TRUE
  )
  check_number(x$artemisinin_loci, paste0(what, "$artemisinin_loci"),
    min = 0, n = NULL, integer = TRUE
  )
  check_number(x$cost_of_resistance, paste0(what, "$cost_of_resistance"),
    min = 0, n = 2^x$number_of_resistance_loci
  )
  # the cost table has one entry per genotype at the drugs' resistance
  # positions, so a run with resistance on needs a locus for every position
  # the drugs name; any() because the flag may be one value per year and the
  # pipeline switches it on when that year starts
  if (any(x$resistance_flag) && x$number_of_resistance_loci < length(x$resistance_loci)) {
    stop(sprintf(
      paste0(
        "%s$number_of_resistance_loci is %d but the drugs have resistance at ",
        "%d barcode positions (%s): set number_of_resistance_loci = %d, with a ",
        "cost_of_resistance for each"
      ),
      what, x$number_of_resistance_loci, length(x$resistance_loci),
      paste(x$resistance_loci, collapse = ", "), length(x$resistance_loci)
    ), call. = FALSE)
  }
  check_flag(x$absolute_fitness_cost_flag, paste0(what, "$absolute_fitness_cost_flag"))
  check_flag(x$mft_flag, paste0(what, "$mft_flag"))
  check_number(x$temporal_cycling, paste0(what, "$temporal_cycling"))
  check_number(x$next_temporal_cycle, paste0(what, "$next_temporal_cycle"))
  check_number(x$sequential_cycling, paste0(what, "$sequential_cycling"))
  check_number(x$sequential_update, paste0(what, "$sequential_update"))
  check_number(x$number_of_drugs, paste0(what, "$number_of_drugs"), min = 1, integer = TRUE)
  check_number(x$drug_choice, paste0(what, "$drug_choice"),
    min = 0, max = x$number_of_drugs - 1, integer = TRUE
  )
  check_number(x$partner_drug_ratios, paste0(what, "$partner_drug_ratios"),
    min = 0, max = 1, n = x$number_of_drugs
  )
  if (!is.list(x$drugs) || length(x$drugs) < x$number_of_drugs) {
    stop(sprintf(
      "%s$drugs must be a list of at least number_of_drugs (%d) drugs", what,
      x$number_of_drugs
    ), call. = FALSE)
  }
  for (i in seq_along(x$drugs)) {
    validate_drug(x$drugs[[i]], what = sprintf("%s$drugs[[%d]]", what, i))
  }
  invisible(x)
}


#' Create list describing parameters for drug efficacy and prophlyaxis
#'
#' @param prob_of_lpf Vector of probabilities for the chance of late
#'   parasitological failure (lpf). The vector gives the prob of lpf for
#'   each relevent barcode combination for a given drug. E.g. The default is:
#'
#'   \code{c(1.0, 0.97, 0.80, 0.55)}
#'
#'   The vector is length 4, and so 2 barcode positions change the
#'   prob of lpf. If the parasite is 0,0 then the prob of lpf is 0 (1-1),
#'   but if it was 0,1 it would be 0.2.
#'
#' @param barcode_res_pos Vector for which barcode positions correspond to which
#'   drug resistance mechanism. E.g. the default is:
#'
#'   \code{c(0,1)}
#'
#'   Resistance to drug is encoded at barcode position 0 and 1
#'
#' @param prophylactic_pos Vector for which barcode positions determine the
#'   impact of drug resistance in shortening the effective prophylactic period.
#'   E.g. the default is \code{c(1)}, which shows that for the drug,
#'   the prophylactic position is encoded in barcode position 1.
#' @param dur_P Duration of prophylaxis in days. Default = 25.
#' @param dur_SPC Duration of slow parasite clearance. Default = 5.
#' @param drug_clearance_max_time Maximum number of days to which to consider
#'   waning prophylaxis. Default = 60 days.
#' @param prophylactic_probability Vector of changing probability of reinfection
#'   due to waning prophylaxis. The last element reflects the probability after
#'   \code{drug_clearance_max_time} and the first element is the probability at
#'   time = 0 days.
#' @param prophylactic_resistant_probability Vector of changing probability of
#'   reinfection due to waning prophylaxis when challenged by a parasite that is
#'  resistant to the partner drug. The last element reflects the probability after
#'   \code{drug_clearance_max_time} and the first element is the probability at
#'   time = 0 days.
#' @return A \code{magenta_drug}: a list describing one therapy, see
#'   \code{\link{cpp_contract}}.
#' @examples
#' # a two-locus drug with 28 day failure rates for each resistance genotype
#' d <- drug_create(prob_of_lpf = c(1, 0.9, 0.8, 0.5), barcode_res_pos = c(0, 1))
#' names(d)
#' d$hill_n
#' @importFrom stats pgamma dexp optim
#' @export
drug_create <- function(prob_of_lpf = c(1.0, 0.97, 0.80, 0.55),
                        barcode_res_pos = c(0, 1),
                        prophylactic_pos = 1,
                        dur_P = 25,
                        dur_SPC = 5,
                        drug_clearance_max_time = 60,
                        prophylactic_probability = 1 - pgamma(
                          seq(0, drug_clearance_max_time, 0.2),
                          shape = 16.8, rate = 16.8 / 17.9
                        ),
                        prophylactic_resistant_probability = 1 - pgamma(
                          seq(0, drug_clearance_max_time, 0.2),
                          shape = 16.8, rate = 16.8 / 8.7
                        )) {
  # Set up for working out hill function parameters for the drug
  time_vec_proph <- seq(0, drug_clearance_max_time,
    length.out = length(prophylactic_probability)
  )
  time_vec_proph_res <- seq(0, drug_clearance_max_time,
    length.out = length(prophylactic_resistant_probability)
  )

  # median of exponential based on the mean for drug elimination
  median_of_prophylactic_duration <- log(2) / (1 / dur_P)
  drug_elimination <- dexp(time_vec_proph, 1 / median_of_prophylactic_duration)
  drug_elimination <- drug_elimination / max(drug_elimination)
  drug_elimination_res <- dexp(time_vec_proph_res, 1 / median_of_prophylactic_duration)
  drug_elimination_res <- drug_elimination_res / max(drug_elimination_res)

  minimise <- function(param, proph, drug_elimination) {
    ## parameters to fit
    n <- param[1]
    conc_fix <- param[2]
    concentration <- drug_elimination

    ## this is the hill function
    pred <- (concentration^n) / ((concentration^n) + (conc_fix^n))
    data <- proph
    sum((data - pred)^2)
  }

  # starting conditions
  param_guess <- c(10, 0.5 * max(drug_elimination))
  hill_model_fit <- optim(param_guess, minimise,
    method = "Nelder-Mead",
    hessian = FALSE,
    proph = prophylactic_probability,
    drug_elimination = drug_elimination
  )


  hill_model_fit_res <- optim(param_guess, minimise,
    method = "Nelder-Mead",
    hessian = FALSE,
    proph = prophylactic_resistant_probability,
    drug_elimination = drug_elimination_res
  )

  # results list that describes properties of our drug
  res_list <- list(
    "lpf" = prob_of_lpf,
    "barcode_positions" = barcode_res_pos,
    "prophylactic_positions" = prophylactic_pos,
    "dur_P" = dur_P,
    "dur_SPC" = dur_SPC,
    "hill_n" = hill_model_fit$par[1],
    "hill_kA" = hill_model_fit$par[2],
    "hill_res_n" = hill_model_fit_res$par[1],
    "hill_res_kA" = hill_model_fit_res$par[2]
  )

  validate_drug(new_param_list(res_list, "magenta_drug"))
}

#' @rdname drug_create
#' @param x List to validate.
#' @param what Name used for the list in error messages.
#' @export
validate_drug <- function(x, what = "drug") {
  check_fields(x, cpp_contract$drug_entry, what)
  check_number(x$lpf, paste0(what, "$lpf"), min = 0, max = 1, n = NULL)
  check_number(x$barcode_positions, paste0(what, "$barcode_positions"),
    min = 0, n = NULL, integer = TRUE
  )
  check_number(x$prophylactic_positions, paste0(what, "$prophylactic_positions"),
    min = 0, n = NULL, integer = TRUE
  )
  for (field in c("dur_P", "dur_SPC", "hill_n", "hill_kA", "hill_res_n", "hill_res_kA")) {
    check_number(x[[field]], paste0(what, "$", field), min = 0)
  }
  invisible(x)
}


#' DHA-PPQ Drug Create
#' @return A \code{magenta_drug}, see \code{\link{drug_create}}.
#' @examples
#' drug_create_dhappq()$lpf
#' @export
drug_create_dhappq <- function() {
  drug_table <- magenta::drug_table

  # https://www.ncbi.nlm.nih.gov/pubmed/25425081
  w_scale <- 28.1
  w_slope <- 4.4
  pip <- exp(-((seq(0, 60, 0.2) / w_scale)^w_slope))

  # based on extrapolating the impact of partner drug resitsance in absence of Kelch
  # i.e. for ASAQ, AL, DHAPPQ difference is c(0.085858, 0.134567, 0.204002) from drug tables.
  # The curves for ASAQ and AL give difference in T0.5 for resistance as c(6, 9.2).
  # Extrapolating then predicts a reduction in time till 0.5 of 13.7612 days
  pip_res <- exp(-((seq(0, 60, 0.2) / w_scale * 2.1)^w_slope))

  drug <- drug_create(
    prob_of_lpf = drug_table$DHAPPQ,
    barcode_res_pos = c(0:5),
    prophylactic_pos = c(5),
    dur_P = 25,
    dur_SPC = 9, # based off average from slater 2016
    drug_clearance_max_time = 60,
    prophylactic_probability = pip,
    prophylactic_resistant_probability = pip_res
  )
  return(drug)
}

#' AL Drug Create
#'
#' @note
#' We have curves of the longest and shortest duration of AL prophylaxis from
#' Bretscher et al. However, we have multiple types of partner drug resistance.
#' At the moment we will say that any loci associated with lumefantrine resistance
#' yield the resistant curve but it may be better to have a distinct curve for each
#' lumefantrine resistant genotype with the T0.5 aligned to the prob of lpf.
#' @return A \code{magenta_drug}, see \code{\link{drug_create}}.
#' @examples
#' drug_create_al()$lpf
#' @export
drug_create_al <- function() {
  drug_table <- magenta::drug_table

  dur_lum_long <- 17.9
  dur_lum_short <- 8.7
  shape_lum <- 93.5


  lum <- 1 - pgamma(seq(0, 60, 0.2), shape = shape_lum, rate = shape_lum / dur_lum_long)
  lum_res <- 1 - pgamma(seq(0, 60, 0.2), shape = shape_lum, rate = shape_lum / dur_lum_short)


  drug <- drug_create(
    prob_of_lpf = drug_table$AL,
    barcode_res_pos = 0:5,
    prophylactic_pos = 0:3,
    dur_P = seq(0, 60, 0.2)[which.min(abs(lum - 0.5))],
    dur_SPC = 6,
    drug_clearance_max_time = 60,
    prophylactic_probability = lum,
    prophylactic_resistant_probability = lum_res
  )
  return(drug)
}

#' ASAQ Drug Create
#'
#' @note
#' We have curves of the longest and shortest duration of ASAQ prophylaxis from
#' Bretscher et al. However, we have multiple types of partner drug resistance.
#' At the moment we will say that any loci associated with AQ resistance
#' yield the resistant curve but it may be better to have a distinct curve for each
#' AQ resistant genotype with the T0.5 aligned to the prob of lpf.
#' @return A \code{magenta_drug}, see \code{\link{drug_create}}.
#' @examples
#' drug_create_asaq()$lpf
#' @export
drug_create_asaq <- function() {
  drug_table <- magenta::drug_table

  dur_aq_long <- 17.8
  dur_aq_short <- 11.6
  shape_aq <- 16.8
  aq <- 1 - pgamma(seq(0, 60, 0.2), shape = shape_aq, rate = shape_aq / dur_aq_long)
  aq_res <- 1 - pgamma(seq(0, 60, 0.2), shape = shape_aq, rate = shape_aq / dur_aq_short)


  drug <- drug_create(
    prob_of_lpf = drug_table$ASAQ,
    barcode_res_pos = 0:5,
    prophylactic_pos = c(0:3),
    dur_P = seq(0, 60, 0.2)[which.min(abs(aq - 0.5))],
    dur_SPC = 6,
    drug_clearance_max_time = 60,
    prophylactic_probability = aq,
    prophylactic_resistant_probability = aq_res
  )
  return(drug)
}

#' Perfect Drug Create
#'
#' @note
#' Perfect Efficacy Drug. Used as default to match deterministic model easily
#' @return A \code{magenta_drug}, see \code{\link{drug_create}}.
#' @examples
#' drug_create_default_no_resistance()$lpf
#' @export
drug_create_default_no_resistance <- function() {
  drug <- list(
    "lpf" = 1,
    "barcode_positions" = 0,
    "prophylactic_positions" = 0,
    "dur_P" = 25,
    "dur_SPC" = 9,
    "hill_n" = 4.415815,
    "hill_kA" = 0.01307829,
    "hill_res_n" = 9.274472,
    "hill_res_kA" = 0.02846058
  )
  validate_drug(new_param_list(drug, "magenta_drug"))
}


#' Create vector adaptation list
#'
#' List for vector adaptations relating to oocyst success
#'
#' @param vector_adaptation_loci Vector of integers detailing which loci
#'   in the barcode correspond to the vector adaptation phenotype
#' @param vector_adaptation_flag Boolean are we doing vector adaptation.
#' @param local_oocyst_advantage Numeric for probability that non adapated
#'   parasites will get through. Default = 0.5
#' @param gametocyte_sterilisation_flag Boolean for whether we are doing
#'   gametocye sterilisation as a result from artemisinin. Default = FALSE
#' @param gametocyte_sterilisation Numeric for theimpact of artemisinin on
#'   male gametocytes. Default = 0.5, which causes the probability that a
#'   wild type male gametocyte will be chosen to be in oocysts is halved.
#' @param oocyst_reduction_by_artemisinin Numeric for reduction in oocyst
#'   under artemisinin drug pressure. default = 0.2
#' @return A \code{magenta_vector_adaptation_list}: a list with the fields the
#'   simulation reads, see \code{\link{cpp_contract}}.
#' @examples
#' va <- vector_adaptation_list_create(
#'   vector_adaptation_loci = 5, vector_adaptation_flag = TRUE,
#'   local_oocyst_advantage = 0.6
#' )
#' va$vector_adaptation_flag
#' @export
vector_adaptation_list_create <- function(vector_adaptation_loci,
                                          vector_adaptation_flag = FALSE,
                                          local_oocyst_advantage = 0.5,
                                          gametocyte_sterilisation_flag = FALSE,
                                          gametocyte_sterilisation = 0.5,
                                          oocyst_reduction_by_artemisinin = 0.2) {
  l <- list(
    "vector_adaptation_flag" = vector_adaptation_flag,
    "vector_adaptation_loci" = vector_adaptation_loci,
    "local_oocyst_advantage" = local_oocyst_advantage,
    "gametocyte_sterilisation_flag" = gametocyte_sterilisation_flag,
    "gametocyte_sterilisation" = gametocyte_sterilisation,
    "oocyst_reduction_by_artemisinin" = oocyst_reduction_by_artemisinin
  )

  validate_vector_adaptation_list(
    new_param_list(l, "magenta_vector_adaptation_list")
  )
}

#' @rdname vector_adaptation_list_create
#' @param x List to validate.
#' @export
validate_vector_adaptation_list <- function(x) {
  what <- "vector_adaptation_list"
  check_fields(x, cpp_contract$vector_adaptation, what)
  check_flag(x$vector_adaptation_flag, paste0(what, "$vector_adaptation_flag"))
  check_number(x$vector_adaptation_loci, paste0(what, "$vector_adaptation_loci"),
    min = 0, n = NULL, integer = TRUE
  )
  check_number(x$local_oocyst_advantage, paste0(what, "$local_oocyst_advantage"),
    min = 0, max = 1
  )
  check_flag(x$gametocyte_sterilisation_flag, paste0(what, "$gametocyte_sterilisation_flag"))
  check_number(x$gametocyte_sterilisation, paste0(what, "$gametocyte_sterilisation"),
    min = 0, max = 1
  )
  check_number(x$oocyst_reduction_by_artemisinin,
    paste0(what, "$oocyst_reduction_by_artemisinin"),
    min = 0, max = 1
  )
  invisible(x)
}


#' Create nmf list
#'
#' List for simulating non malarial fever
#'
#' @param nmf_flag Boolean are we doing non malarial fevers
#' @param mean_nmf_frequency Vector for mean number of days between fevers for
#'   the age bracket considered
#' @param nmf_age_brackets Vector for age brackets
#' @param prob_of_testing_nmf Numeric for probability that a NMF is tested by
#'   RDT before being treated with antimalarials.
#' @return A \code{magenta_nmf_list}: a list with the fields the simulation
#'   reads, see \code{\link{cpp_contract}}.
#' @examples
#' nmf <- nmf_list_create(nmf_flag = TRUE, prob_of_testing_nmf = 0.5)
#' nmf$nmf_flag
#' @export
nmf_list_create <- function(nmf_flag = FALSE,
                            mean_nmf_frequency = c(
                              148.578, 139.578, 141.564, 155.874, 179.364,
                              216.192, 233.478, 268.056, 312.858, 315.564,
                              285.156, 255.246, 238.302, 216.618
                            ),
                            nmf_age_brackets = c(
                              -0.1, 365.0, 730.0, 1095.0, 1460.0, 1825.0,
                              2555.0, 3285.0, 4015.0, 4745.0, 5475.0,
                              7300.0, 9125.0, 10950.0, 36850.0
                            ),
                            prob_of_testing_nmf = 0.5) {
  l <- list(
    "nmf_flag" = nmf_flag,
    "mean_nmf_frequency" = mean_nmf_frequency,
    "nmf_age_brackets" = nmf_age_brackets,
    "prob_of_testing_nmf" = prob_of_testing_nmf
  )

  validate_nmf_list(new_param_list(l, "magenta_nmf_list"))
}

#' @rdname nmf_list_create
#' @param x List to validate.
#' @export
validate_nmf_list <- function(x) {
  what <- "nmf_list"
  check_fields(x, cpp_contract$nmf, what)
  check_flag(x$nmf_flag, paste0(what, "$nmf_flag"))
  check_number(x$mean_nmf_frequency, paste0(what, "$mean_nmf_frequency"),
    min = 0, n = NULL
  )
  check_number(x$nmf_age_brackets, paste0(what, "$nmf_age_brackets"),
    n = length(x$mean_nmf_frequency) + 1
  )
  check_number(x$prob_of_testing_nmf, paste0(what, "$prob_of_testing_nmf"),
    min = 0, max = 1
  )
  invisible(x)
}


#' Barcode parameter list creation
#'
#' \code{barcode_list_create} creates list detialing the barcode/genetic
#'   parameters for the model
#'
#' @param num_loci Number of loci. Default = 24
#' @param ibd_length If we are simulating IBD dynamics, each loci is now
#'   represented by a bitset of ibd_length. Thus ibd_length needs to be long
#'   enough to ensure that as new identity relationships occur, i.e.
#'   an importation barcode will be a new identity. e.g. If your population is
#'   1000, we may expect at 80% prevalence, with a mean COI of 3 we will need
#'   2400 different identities, i.e. 2^ibd_length > 2400. However, keep in mind
#'   importations as these need to be continually new, i.e. if we are simulating
#'   for 30 years, with 3 importations a day, then we will need at least length
#'   to ensure that 2^ibd_length > 2400 + (30*365*3). This will probably be
#'   automatically calculated in the future. If we are not ding IBD, then this
#'   should be 1, which is the defalt.
#' @param plaf Vector of population level allele frequencies for the barcode.
#'   Default = rep(0.5, num_loci)
#' @param prob_crossover Vector of probabilities for crossover events for the
#'   barcode. Default = rep(0.5, num_loci)
#' @param starting_ibd Starting ibd. Default = 0
#' @param mutation_flag Boolean for simulating mutations
#' @param mutation_rate Probability of mutation occuring and fixing
#' @param mutation_treated_modifier Multiplier for how much more likely mutations
#'   are to occur in treated individuals with respect to resistance. Default = 1, i.e
#'   no difference
#' @return A \code{magenta_barcode_list}: a list with the fields the
#'   simulation reads, see \code{\link{cpp_contract}}.
#' @examples
#' bl <- barcode_list_create(num_loci = 6, plaf = rep(0.2, 6), mutation_flag = TRUE)
#' bl$num_loci
#' validate_barcode_list(bl)
#' @export
barcode_list_create <- function(num_loci = 24,
                                ibd_length = 1,
                                plaf = rep(0.5, num_loci),
                                prob_crossover = rep(0.5, num_loci),
                                starting_ibd = 0,
                                mutation_flag = FALSE,
                                mutation_rate = rep(1e-7, num_loci),
                                mutation_treated_modifier = 1.0) {
  # are we doing ibd or not
  if (ibd_length > 1) {
    barcode_type <- TRUE
  } else {
    barcode_type <- FALSE
  }

  res <- list(
    "num_loci" = num_loci,
    "ibd_length" = ibd_length,
    "starting_ibd" = starting_ibd,
    "plaf" = plaf,
    "prob_crossover" = prob_crossover,
    "barcode_type" = barcode_type,
    "mutation_flag" = mutation_flag,
    "mutation_rate" = mutation_rate,
    "mutation_treated_modifier" = mutation_treated_modifier
  )

  # catch for changing mutation rates etc over time
  if (length(res$mutation_flag) > 1) {
    res$mutation_flag <- res$mutation_flag[1]
  }

  validate_barcode_list(new_param_list(res, "magenta_barcode_list"))
}

#' @rdname barcode_list_create
#' @param x List to validate.
#' @export
validate_barcode_list <- function(x) {
  what <- "barcode_list"
  check_fields(x, cpp_contract$barcode, what)
  check_number(x$num_loci, paste0(what, "$num_loci"), min = 1, integer = TRUE)
  check_number(x$ibd_length, paste0(what, "$ibd_length"), min = 1, integer = TRUE)
  check_number(x$starting_ibd, paste0(what, "$starting_ibd"), min = 0, max = 1)
  check_number(x$plaf, paste0(what, "$plaf"), min = 0, max = 1, n = x$num_loci)
  check_number(x$prob_crossover, paste0(what, "$prob_crossover"),
    min = 0, max = 1, n = x$num_loci
  )
  check_flag(x$barcode_type, paste0(what, "$barcode_type"))
  check_flag(x$mutation_flag, paste0(what, "$mutation_flag"))
  check_number(x$mutation_rate, paste0(what, "$mutation_rate"),
    min = 0, max = 1, n = x$num_loci
  )
  check_number(x$mutation_treated_modifier, paste0(what, "$mutation_treated_modifier"),
    min = 0
  )
  invisible(x)
}

# The contract between the R parameter lists and the C++ entry points.
#
# `cpp_contract` records, for every list R hands to C++, the field names in
# the order the R constructor produces them (`fields`) and the subset that the
# C++ side reads (`cpp`; see src/main_init.cpp, src/main_update.cpp and
# src/main_saved_init.cpp). The validators below check a list against it, and
# tests/testthat/test-params.R checks that each constructor produces exactly
# the `fields` set, so a change on either side fails a test rather than
# failing inside C++.

#' Field names read by the C++ simulation
#'
#' \code{cpp_contract} lists, for every parameter list passed from R to the
#' C++ simulation, the fields the R constructor produces (\code{fields}, in
#' order) and the subset the C++ side reads (\code{cpp}). The
#' \code{validate_*()} functions check a list against it and the package
#' tests check that every constructor produces exactly the recorded set, so a
#' change on either side is caught by a test rather than inside C++.
#'
#' @format A named list. The entries \code{init}, \code{update},
#'   \code{saved_init}, \code{get}, \code{finalizer} and \code{saved_state}
#'   are character vectors naming the top-level fields of each entry point;
#'   \code{housekeeping}, \code{barcode}, \code{drug}, \code{drug_entry},
#'   \code{vector_adaptation}, \code{nmf}, \code{spatial}, \code{eqss} and
#'   \code{core} are lists with \code{fields} and \code{cpp}.
#' @export
#' @name cpp_contract
#' @examples
#' cpp_contract$drug$cpp
#' setequal(names(drug_list_create()), cpp_contract$drug$fields)

NULL

cpp_contract <- list(
  # top-level lists for each entry point
  init = c(
    "N", "eqSS", "barcode_list", "spatial_list", "housekeeping_list",
    "drug_list", "nmf_list", "vector_adaptation_list", "core_parameter_list"
  ),
  update = c(
    "years", "ft", "mu_vec", "fv_vec", "statePtr",
    "spatial_list", "drug_list", "barcode_list"
  ),
  saved_init = "savedState",
  get = "statePtr",
  finalizer = c("statePtr", "finalizer"),
  saved_state = c(
    "population_List", "populations_event_and_strains_List", "scourge_List",
    "parameters_List", "parameter_housekeeping_List"
  ),
  # the parameter lists themselves
  housekeeping = list(
    fields = c("quiet_print", "quiet_test_print", "cluster", "clear_up"),
    cpp = c("quiet_print", "quiet_test_print")
  ),
  barcode = list(
    fields = c(
      "num_loci", "ibd_length", "starting_ibd", "plaf", "prob_crossover",
      "barcode_type", "mutation_flag", "mutation_rate",
      "mutation_treated_modifier"
    ),
    cpp = c(
      "num_loci", "ibd_length", "starting_ibd", "plaf", "prob_crossover",
      "barcode_type", "mutation_flag", "mutation_rate",
      "mutation_treated_modifier"
    )
  ),
  drug = list(
    fields = c(
      "resistance_flag", "number_of_resistance_loci", "resistance_loci",
      "artemisinin_loci", "cost_of_resistance", "absolute_fitness_cost_flag",
      "drugs", "mft_flag", "temporal_cycling", "next_temporal_cycle",
      "sequential_cycling", "sequential_update", "number_of_drugs",
      "drug_choice", "partner_drug_ratios"
    ),
    cpp = c(
      "resistance_flag", "number_of_resistance_loci", "resistance_loci",
      "artemisinin_loci", "cost_of_resistance", "absolute_fitness_cost_flag",
      "drugs", "mft_flag", "number_of_drugs", "drug_choice",
      "partner_drug_ratios"
    )
  ),
  drug_entry = list(
    fields = c(
      "lpf", "barcode_positions", "prophylactic_positions", "dur_P", "dur_SPC",
      "hill_n", "hill_kA", "hill_res_n", "hill_res_kA"
    ),
    cpp = c(
      "lpf", "barcode_positions", "prophylactic_positions", "dur_P", "dur_SPC",
      "hill_n", "hill_kA", "hill_res_n", "hill_res_kA"
    )
  ),
  vector_adaptation = list(
    fields = c(
      "vector_adaptation_flag", "vector_adaptation_loci",
      "local_oocyst_advantage", "gametocyte_sterilisation_flag",
      "gametocyte_sterilisation", "oocyst_reduction_by_artemisinin"
    ),
    cpp = c(
      "vector_adaptation_flag", "vector_adaptation_loci",
      "local_oocyst_advantage", "gametocyte_sterilisation_flag",
      "gametocyte_sterilisation", "oocyst_reduction_by_artemisinin"
    )
  ),
  nmf = list(
    fields = c(
      "nmf_flag", "mean_nmf_frequency", "nmf_age_brackets",
      "prob_of_testing_nmf"
    ),
    cpp = c(
      "nmf_flag", "mean_nmf_frequency", "nmf_age_brackets",
      "prob_of_testing_nmf"
    )
  ),
  spatial = list(
    fields = c(
      "spatial_type", "cotransmission_freq_vector", "oocyst_freq_vector",
      "imported_cotransmissions_events", "imported_oocyst_events", "plaf",
      "island_imports_plaf_linked_flag"
    ),
    cpp = c(
      "spatial_type", "cotransmission_freq_vector", "oocyst_freq_vector",
      "imported_cotransmissions_events", "imported_oocyst_events", "plaf",
      "island_imports_plaf_linked_flag"
    )
  ),
  eqss = list(
    fields = c(
      "age_brackets", "het_brackets", "Smat", "Dmat", "Amat", "Umat", "Tmat",
      "Pmat", "IBmat", "ICAmat", "ICMmat", "IDmat", "FOI", "phi", "Sv", "Ev",
      "Iv", "ICM_Init", "theta"
    ),
    cpp = c(
      "age_brackets", "het_brackets", "Smat", "Dmat", "Amat", "Umat", "Tmat",
      "Pmat", "IBmat", "ICAmat", "ICMmat", "IDmat", "FOI", "phi", "Sv", "Ev",
      "Iv", "ICM_Init", "theta"
    )
  ),
  # the core (epidemiological) parameter list. C++ reads phi0 only; the rest
  # is consumed by equilibrium_init_create() and the odin model.
  core = list(
    fields = c(
      "DY", "eta", "rho", "a0", "sigma2", "max_age",
      "rA", "rT", "rD", "rU", "rP", "dE", "delayGam",
      "cD", "cT", "cU", "gamma1",
      "d1", "dID", "ID0", "kD", "uD", "aD", "fD0", "gammaD",
      "alphaA", "alphaU",
      "b0", "b1", "dB", "IB0", "kB", "uB",
      "phi0", "phi1", "dCA", "IC0", "kC", "uCA", "PM", "dCM",
      "delayMos", "tau1", "tau2", "mu0", "Q0", "chi", "bites_Bed",
      "bites_Indoors", "fv0", "av0", "Surv0", "p10", "p2",
      "muEL", "muLL", "muPL", "dEL", "dLL", "dPL", "gammaL", "km", "cm",
      "betaL", "eov", "b_lambda", "lambda",
      "num_int", "itn_cov", "irs_cov", "ITN_IRS_on", "pop_split", "cov",
      "d_ITN0", "r_ITN0", "r_ITN1", "r_IRS0", "d_IRS0",
      "irs_half_life", "itn_half_life", "IRS_interval", "ITN_interval",
      "irs_loss", "itn_loss"
    ),
    cpp = "phi0"
  )
)

# ---- generic checks ---------------------------------------------------------
#
# Each check reports the list and field it is looking at, e.g.
# "drug_list$drug_choice must be a single number", and returns the value
# invisibly so the validators read as a list of requirements.

#' Check a parameter list against the contract
#'
#' The generic checks used by every \code{validate_*()}. Each reports the
#' list and field it is looking at, e.g.
#' \code{"drug_list$drug_choice must be a single number"}, and returns the
#' value invisibly so the validators read as a list of requirements.
#'
#' @param x The list.
#' @param contract The entry of \code{\link{cpp_contract}} for it.
#' @param what Name used for the list in error messages.
#' @return \code{x}, invisibly; an error names the missing or unexpected
#'   fields.
#' @keywords internal
check_fields <- function(x, contract, what) {
  if (!is.list(x)) {
    stop(sprintf("%s must be a list", what), call. = FALSE)
  }
  missing <- setdiff(contract$fields, names(x))
  extra <- setdiff(names(x), contract$fields)
  if (length(missing)) {
    stop(sprintf(
      "%s is missing the field(s): %s", what, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (length(extra)) {
    stop(sprintf(
      "%s has unexpected field(s): %s", what, paste(extra, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(x)
}

#' Check for a logical without NAs
#'
#' @param x The value.
#' @param what Name used in the error message.
#' @param n Required length, or \code{NULL} for any.
#' @return \code{x}, invisibly.
#' @keywords internal
check_flag <- function(x, what, n = 1L) {
  if (!is.logical(x) || anyNA(x) || (!is.null(n) && length(x) != n)) {
    stop(sprintf(
      "%s must be %s", what,
      if (identical(n, 1L)) "TRUE or FALSE" else "a logical vector without NAs"
    ), call. = FALSE)
  }
  invisible(x)
}

#' Check for a number, or a vector of them, in a range
#'
#' @param x The value.
#' @param what Name used in the error message.
#' @param min,max Bounds, inclusive.
#' @param n Required length, or \code{NULL} for any.
#' @param integer Whether whole numbers are required (logicals pass too).
#' @return \code{x}, invisibly.
#' @keywords internal
check_number <- function(x, what, min = -Inf, max = Inf, n = 1L,
                         integer = FALSE) {
  kind <- if (integer) "whole number" else "number"
  if (!(is.numeric(x) || (integer && is.logical(x))) || anyNA(x)) {
    stop(sprintf(
      "%s must be %s", what,
      if (identical(n, 1L)) paste0("a single ", kind) else paste0("a ", kind, " vector")
    ), call. = FALSE)
  }
  if (!is.null(n)) {
    check_length(x, what, n)
  }
  in_range <- all(x >= min) && all(x <= max) && (!integer || all(x == round(x)))
  if (!in_range) {
    range <- if (is.finite(min) && is.finite(max)) {
      sprintf(" between %s and %s", format(min), format(max))
    } else if (is.finite(min)) {
      sprintf(" of at least %s", format(min))
    } else if (is.finite(max)) {
      sprintf(" of at most %s", format(max))
    } else {
      ""
    }
    stop(sprintf(
      "%s must be %s%s", what,
      if (identical(n, 1L)) paste0("a single ", kind) else paste0("a ", kind, " vector"),
      range
    ), call. = FALSE)
  }
  invisible(x)
}

#' Check a length
#'
#' @param x The value.
#' @param what Name used in the error message.
#' @param n Required length.
#' @return \code{x}, invisibly.
#' @keywords internal
check_length <- function(x, what, n) {
  if (length(x) != n) {
    stop(sprintf("%s must have length %d, not %d", what, n, length(x)),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Check that a list came from its constructor
#'
#' @param x The value.
#' @param class The class the constructor gives, e.g.
#'   \code{"magenta_drug_list"}.
#' @param what Name used in the error message.
#' @return \code{x}, invisibly; the error names the constructor to use.
#' @keywords internal
check_class <- function(x, class, what) {
  if (!inherits(x, class)) {
    stop(sprintf(
      "%s must be created by %s()", what,
      sub("^magenta_", "", class) |> paste0("_create")
    ), call. = FALSE)
  }
  invisible(x)
}

#' Stamp a parameter list with its class
#'
#' The list stays a plain list underneath, so that the C++ side and existing
#' code that treats it as a list are unaffected.
#'
#' @param x The list.
#' @param class The class, e.g. \code{"magenta_drug_list"}.
#' @return \code{x} with class \code{c(class, "magenta_param_list", "list")}.
#' @keywords internal
new_param_list <- function(x, class) {
  structure(x, class = c(class, "magenta_param_list", "list"))
}

#' Print a parameter list
#'
#' @param x A parameter list from one of the \code{*_list_create()}
#'   constructors.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @keywords internal
#' @export
print.magenta_param_list <- function(x, ...) {
  cat("<", class(x)[1], "> with ", length(x), " fields\n", sep = "")
  utils::str(unclass(x), max.level = 1, give.attr = FALSE, no.list = TRUE)
  invisible(x)
}
