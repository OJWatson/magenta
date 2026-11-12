#------------------------------------------------
#' Brute force conversion of PfPR and ft to an EIR
#'
#' @param ft Treatment seeking. Default = 0
#' @param PfPR Parasite prevaence by PCR
#' @param PfPR_micro Parasite prevaence by microscopy. Default = NULL
#' @param mv Mosquito size if it's known. Default = NULL
#' @param age_range Numeric vecotr length 2 for the age range to be included.
#'   Default = NULL
#' @param all_a Boolean for whether all asymptomatics in state A are included
#'   in microscopy based prevalence estimates. Default = NULL
#' @param ... Any other params to be fed to the model parameter list
#'
#' @keywords internal
pfpr_to_eir_heuristic <- function(ft = 0,
                                  PfPR = NULL, PfPR_micro = NULL, mv = NULL,
                                  age_range = NULL, all_a = FALSE,
                                  ...) {
  mpl <- model_param_list_create(eta = 1 / (21 * 365), ...)

  num_het_brackets <- 5
  num_age_brackets <- 20
  max_age <- 100
  ## Create the geometric age brackets
  ratio <- (max_age / 0.1)^(1 / num_age_brackets)
  age.vector <- 0.1 * ratio**(1:num_age_brackets)
  age.vector[1] <- 0

  bm <- read.csv(system.file("extdata/bm.txt", package = "magenta"),
    sep = ",",
    header = T
  )

  if (!is.null(PfPR)) {
    PR <- PfPR
    EIR <- bm$EIRY_eq[which.min(abs(bm$pcr_pos_all.final..1 - PR))]
  } else if ((!is.null(PfPR_micro))) {
    PR <- PfPR_micro
    EIR <- bm$EIRY_eq[which.min(abs(bm$slide_pos_2_10 - PR))]
  } else {
    PR <- mv
    EIR <- bm$EIRY_eq[which.min(abs(bm$mv - PR))]
  }


  not_found <- TRUE
  if (is.null(age_range)) {
    age_low <- 1
    age_high <- length(age.vector)
  } else {
    age_low <- which(eqInit$age > (age_range[1] * 365))[1]
    age_high <- which(eqInit$age > (age_range[2] * 365))[1]
  }

  # while(not_found){
  for (i in 1:100) {
    ## Create a near equilibirum initial condition
    eqInit <- equilibrium_init_create(
      age_vector = age.vector,
      het_brackets = num_het_brackets,
      ft = ft,
      EIR = EIR,
      model_param_list = mpl
    )

    if (!is.null(PfPR)) {
      diff <- PR - (1 - sum(eqInit$init_S[age_low:age_high, , 1] + eqInit$init_P[age_low:age_high, , 1]) / sum(eqInit$den[age_low:age_high]))
    } else if ((!is.null(PfPR_micro))) {
      if (all_a) {
        diff <- PR - (sum(eqInit$init_D[age_low:age_high, , 1] +
          eqInit$init_T[age_low:age_high, , 1] +
          eqInit$init_A[age_low:age_high, , 1]) /
          sum(eqInit$den[age_low:age_high]))
      } else {
        diff <- PR - (sum(eqInit$init_D[age_low:age_high, , 1] +
          eqInit$init_T[age_low:age_high, , 1] +
          eqInit$init_A[age_low:age_high, , 1] * eqInit$p_det_eq[age_low:age_high, ]) /
          sum(eqInit$den[age_low:age_high]))
      }
    } else {
      diff <- PR - eqInit$mv0
    }

    if (abs(diff) < 0.001) {
      not_found <- FALSE
    } else {
      if (diff < 0) {
        if (abs(diff) > 0.01) {
          EIR <- EIR - 0.2
        } else if (abs(diff) > 0.005) {
          EIR <- EIR - 0.1
        } else if (abs(diff) > 0.0025) {
          EIR <- EIR - 0.05
        } else {
          EIR <- EIR - 0.001
        }
      } else {
        if (abs(diff) > 0.01) {
          EIR <- EIR + 0.2
        } else if (abs(diff) > 0.005) {
          EIR <- EIR + 0.1
        } else if (abs(diff) > 0.0025) {
          EIR <- EIR + 0.05
        } else {
          EIR <- EIR + 0.001
        }
      }
    }
  }

  message(abs(diff))
  return(EIR)
}
