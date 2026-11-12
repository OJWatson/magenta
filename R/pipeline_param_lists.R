#' Create housekeeping parameter list
#'
#' List for simulation housekeeping vars, e.g. quiet prints,
#'
#' @param quiet Boolean for quiet simulation. Default = TRUE
#' @param quiet_test Boolean for quiet testing statement. Default = TRUE
#' @param cluster Boolean for simulation being on cluster. Default = TRUE
#' @param clear_up Boolean for whether to clear up the memory used by the
#'   simulation. Default = TRUE

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

  return(l)
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
#' @param number_of_resistance_loci Numeric for number of res. loci
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
#'
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
    barcode_poss <- matrix(intToBits(0:max_bits), ncol = 32, byrow = T)[, 1:(number_of_resistance_loci)]
    results <- rep(1, nrow(barcode_poss))
    count <- 1

    for (j in 1:nrow(barcode_poss)) {
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

  return(l)
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
#' @importFrom stats pgamma dexp optim
drug_create <- function(prob_of_lpf = c(1.0, 0.97, 0.80, 0.55),
                        barcode_res_pos = c(0, 1),
                        prophylactic_pos = 1,
                        dur_P = 25,
                        dur_SPC = 5,
                        drug_clearance_max_time = 60,
                        prophylactic_probability = 1 - pgamma(seq(0, drug_clearance_max_time, 0.2), shape = 16.8, rate = 16.8 / 17.9),
                        prophylactic_resistant_probability = 1 - pgamma(seq(0, drug_clearance_max_time, 0.2), shape = 16.8, rate = 16.8 / 8.7)) {
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
    output <- sum((data - pred)^2)
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

  return(res_list)
}


#' DHA-PPQ Drug Create
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
drug_create_default_no_resistance <- function() {
  data.frame(
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
#'

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

  return(l)
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
#'
#'

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

  return(l)
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
#'   Default = rep(0.5, barcode_length)
#' @param prob_crossover Vector of probabilities for crossover events for the
#'   barcode. Default = rep(0.5, barcode_length)
#' @param starting_ibd Starting ibd. Default = 0
#' @param mutation_flag Boolean for simulating mutations
#' @param mutation_rate Probability of mutation occuring and fixing
#' @param mutation_treated_modifier Multiplier for how much more likely mutations
#'   are to occur in treated individuals with respect to resistance. Default = 1, i.e
#'   no difference
#'
#' @keywords internal
barcode_list_create <- function(num_loci = 24,
                                ibd_length = 1,
                                plaf = rep(0.5, 24),
                                prob_crossover = rep(0.5, 24),
                                starting_ibd = 0,
                                mutation_flag = FALSE,
                                mutation_rate = rep(1e-7, num_loci),
                                mutation_treated_modifier = 1.0) {
  if (length(mutation_rate) != num_loci) {
    stop("Not enough mutation rates. Should be num_loci length")
  }

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

  return(res)
}
