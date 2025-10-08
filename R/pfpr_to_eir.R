#------------------------------------------------
#' Brute force conversion of PfPR and ft to an EIR
#'
#' @param ft Treatment seeking. Default = 0
#' @param PfPR Parasite prevaence by PCR
#' @param PfPR_micro Parasite prevaence by microscopy. Default = NULL
#' @param mv Mosquito size if it's known. Default = NULL
#' @param age_range Numeric vector of length 2 giving the lower and upper ages
#'   (in years) to include when matching the target prevalence. Use `NULL` to
#'   include all ages. Default = NULL.
#' @param all_a Boolean for whether all asymptomatics in state A are included  
#'   in microscopy based prevalence estimates. Default = NULL
#' @param ... Any other params to be fed to the model parameter list
#' 
#' @keywords internal
pfpr_to_eir_heuristic <- function(ft = 0,
                                  PfPR = NULL, PfPR_micro = NULL, mv = NULL,
                                  age_range = NULL, all_a = FALSE,
                                  ...){

  mpl <- model_param_list_create(eta = 1/(21*365), ...)

  num_het_brackets <- 5
  num_age_brackets <- 20
  max_age <- 100
  ## Create the geometric age brackets
  ratio <- (max_age/0.1)^(1/num_age_brackets)
  age.vector <- 0.1 * ratio ** (1:num_age_brackets)
  age.vector[1] <- 0
  
  bm <- read.csv(system.file("extdata/bm.txt", package = "magenta"),
                 sep=",",
                 header = T)
    
  if (!is.null(age_range)) {
    if (!is.numeric(age_range) || length(age_range) != 2 || anyNA(age_range)) {
      stop("`age_range` must be a numeric vector of length 2 without missing values")
    }
    if (any(age_range < 0)) {
      stop("`age_range` must contain non-negative ages")
    }
    age_range <- sort(age_range)
    if (age_range[1] == age_range[2]) {
      warning("`age_range` lower and upper bounds are identical; a single age band will be matched")
    }
  }

  if (!is.null(PfPR)) {
    PR <- PfPR
    EIR <- bm$EIRY_eq[which.min(abs(bm$pcr_pos_all.final..1 - PR))]
  } else if (!is.null(PfPR_micro)) {
    PR <- PfPR_micro
    EIR <- bm$EIRY_eq[which.min(abs(bm$slide_pos_2_10 - PR))]
  } else {
    PR <- mv
    EIR <- bm$EIRY_eq[which.min(abs(bm$mv - PR))]
  }

  age_idx <- NULL
  diff <- NA_real_

  for (i in seq_len(100)) {

    ## Create a near equilibirum initial condition
    eqInit <- equilibrium_init_create(age_vector = age.vector,
                                      het_brackets = num_het_brackets,
                                      ft = ft,
                                      EIR = EIR,
                                      model_param_list = mpl)

    if (is.null(age_idx)) {
      if (is.null(age_range)) {
        age_idx <- seq_len(length(eqInit$age))
      } else {
        age_days <- eqInit$age
        lower_idx <- which(age_days >= (age_range[1] * mpl$DY))[1]
        if (is.na(lower_idx)) {
          lower_idx <- length(age_days)
        }
        upper_candidates <- which(age_days <= (age_range[2] * mpl$DY))
        if (length(upper_candidates) == 0) {
          upper_idx <- lower_idx
        } else {
          upper_idx <- max(upper_candidates)
          if (upper_idx < lower_idx) {
            upper_idx <- lower_idx
          }
        }
        age_idx <- seq.int(lower_idx, upper_idx)
      }
    }

    if (!is.null(PfPR)) {
      diff <- PR - (1 - sum(eqInit$init_S[age_idx,,1] + eqInit$init_P[age_idx,,1])/
                      sum(eqInit$den[age_idx]))
    } else if (!is.null(PfPR_micro)) {
      if (all_a) {
        diff <- PR - (sum(eqInit$init_D[age_idx,,1] +
                            eqInit$init_T[age_idx,,1] +
                            eqInit$init_A[age_idx,,1]) /
                        sum(eqInit$den[age_idx]))
      } else {
        diff <- PR - (sum(eqInit$init_D[age_idx,,1] +
                            eqInit$init_T[age_idx,,1] +
                            eqInit$init_A[age_idx,,1] * eqInit$p_det_eq[age_idx,]) /
                        sum(eqInit$den[age_idx]))
      }
    } else {
      diff <- PR - eqInit$mv0
    }

    if (abs(diff) < 0.001) {
      break
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

  message(sprintf("Final prevalence difference: %.4f", abs(diff)))
  return(EIR)

}