###
#' Get ppopulation strain information
#'
#' @details R side wrapper for c++ functions that  loop through the population
#'   retrieving information about their strains using the universe state ptr.
#'   The resulting data frame is then passed to \code{\link{COI_df_create}}.
#'
#' @param statePtr Universe state pointer
#' @param seed The random seed. Probably not needed but was here for debugging
#'   radnomness.
#' @param sample_size Number of indiviudals to sample. Default = 0, which
#'   samples everyone.
#' @param sample_states Numeric vector for which sample states are to be sampled
#'   from Default = 0:5, which is all states.
#' @param ibd Boolean for whether the siulation is an IBD one.
#'
#' @keywords internal
pop_strains_df <- function(statePtr, seed, sample_size = 0, sample_states = 0:5, ibd = FALSE) {
  param_list <- list(statePtr = statePtr, sample_size = sample_size, sample_states = sample_states)
  # message('pop_strains\n')
  if (!ibd) {
    # Get the info
    set.seed(seed)
    list <- population_get_genetics_df_n(param_list)
    df <- as.data.frame.list(list[1:4])
  } else {
    # Get the info
    set.seed(seed)
    list <- population_get_genetics_ibd_df_n(param_list)
    df <- as.data.frame.list(list[1:6])
  }

  # bring barcode states in
  df$barcode_states <- list$barcode_states

  if (dim(df)[1] == 0) {
    df$nums <- df$barcode_states
  } else {
    # split the barcodes into a list and then move
    seqs <- unlist(lapply(list$barcode_states, length))
    breaks <- cumsum(seqs)
    nums <- list()
    length(nums) <- length(df$barcode_states)

    if (seqs[1] > 0) {
      nums[[1]] <- list$barcodes[1:breaks[1], , drop = FALSE]
    }

    if (length(seqs) > 1) {
      for (i in 2:(length(breaks))) {
        if (seqs[i] > 1) {
          nums[[i]] <- list$barcodes[(breaks[i - 1] + 1):breaks[i], ]
        }
        if (seqs[i] == 1) {
          nums[[i]] <- list$barcodes[(breaks[i - 1] + 1):breaks[i], , drop = FALSE]
        }
      }
    }
    df$nums <- nums
  }

  return(df)
}


###
#---
#' Create COI dataframe
#'
#' @param df Output of \code{\link{pop_strains_df}}
#' @param groupvars Grouping vars for summarySE.
#'   Default = \code{c('age_bin','clinical')}
#' @param breaks Numberic vector of age breaks.
#'   Default = \code{c(-0.001,5,15,100.1)}
#' @param barcodes Boolean whether to return tabled barcodes. Default = FALSE
#' @param ibd Boolean for ibd simulations or not. Default = FALSE
#' @param n Number of individuals to sample when creating COI data frame.
#'   Default = Inf that means all individuals are included
#' @param nl Number of loci. Default = 24
#' @param reps Numeric for number of population subsamples to make. Default = 1
#' @param mean_only Boolean for whether to just return the mean when summarising
#'   vs all summary statistics. Default = TRUE
#' @keywords internal
COI_df_create <- function(df,
                          groupvars = c("age_bin", "clinical"),
                          breaks = c(-0.001, 5, 15, 100.1),
                          barcodes = FALSE,
                          ibd = FALSE,
                          n = Inf,
                          reps = 1,
                          nl = 24,
                          mean_only = TRUE) {
  # handle df for either genetic type first
  infection_state <- c("S", "D", "A", "U", "T", "P")
  df$age_bin <- cut(df$age / 365, breaks = breaks)
  df$state <- df$clinical <- infection_state[df$state + 1]
  df$clinical[df$clinical %in% c("D", "T")] <- "Clinical"
  df$clinical[df$clinical %in% c("A", "U")] <- "Asymptomatic"
  df$bs <- unlist(lapply(df$barcode_states, length))

  nrow_df <- nrow(df)

  # if doing samples
  if (length(n) == 1) {
    if (n != Inf) {
      # if we have asked for more samples than there are then sample from them
      if ((reps * n) > nrow_df) {
        if (n > nrow_df) {
          ranges <- lapply(seq_len(reps), function(x) sample(nrow_df, n, TRUE))
        } else {
          ranges <- lapply(seq_len(reps), function(x) sample(nrow_df, n, FALSE))
        }
      } else {
        sample <- sample(reps * n)
        ranges <- ranges(diff = n, end = reps * n)
        ranges <- lapply(ranges, function(x) sample[x])
      }
    }
  } else {
    reps <- length(n) + 1
    ranges <- lapply(n, function(x) {
      if (x > nrow_df) {
        sample(nrow_df, x, TRUE)
      } else {
        sample(nrow_df, x, FALSE)
      }
    })
    ranges[[length(ranges) + 1]] <- seq_len(nrow_df)
  }


  # non ibd simulation summary
  if (!ibd) {
    # polygenomic return
    df$polygenom <- df$coi > 1

    # set up results
    results <- list()
    length(results) <- reps

    # create overall clonality by first converting our nums to integers for quick tabulation
    barcode_freq <- apply(rbind_list_base(df$nums), 2, function(x) {
      mean(as.integer(x))
    })

    df$nums[df$bs > 0] <- lapply(df$nums[df$bs > 0], function(x) {
      apply(x, 1, bitsToInt)
    })

    # make samples and reps of them
    for (i in seq_len(length(results))) {
      # make a sample
      if (n[1] == Inf) {
        chosen <- seq_len(nrow_df)
      } else {
        chosen <- ranges[[i]]
      }


      clonality <- table(table(unlist(df$nums[chosen])))
      barcodes_tab <- sort(table(unlist(df$nums[chosen])),
        decreasing = TRUE
      )

      # create summary dfs
      # ------------------------------------------------------------------

      summary_ages <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "age",
        groupvars = groupvars,
        mean_only = mean_only
      )

      summary_coi <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "coi",
        groupvars = groupvars,
        mean_only = mean_only
      )

      summary_coi_detected <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "coi_detected_pcr",
        groupvars = groupvars,
        mean_only = mean_only
      )

      summary_moi_detected <-
        summarySE_mean_only_max_mean(
          df[chosen, ],
          measurevar = "coi_detected_pcr",
          groupvars = groupvars,
          mean_only = mean_only,
          max = 6
        )

      summary_polygenom <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "polygenom",
        groupvars = groupvars,
        mean_only = mean_only
      )

      # unique needs to be calculated per groups
      summary_unique1 <- dplyr::group_by(df[chosen, ], .data$age_bin) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = clonality_from_barcode_list(.data$nums)
        )

      summary_unique2 <- dplyr::group_by(df[chosen, ], .data$clinical) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = clonality_from_barcode_list(.data$nums)
        )

      summary_unique3 <- dplyr::group_by(df[chosen, ]) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = clonality_from_barcode_list(.data$nums)
        )

      unq_list <- list(summary_unique1, summary_unique2, summary_unique3)
      summary_unique <- dplyr::bind_rows(unq_list)


      # cou needs to be calculated per groups
      summary_cou1 <- dplyr::group_by(df[chosen, ], .data$age_bin) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = cou_from_barcode_list(.data$nums)
        )
      summary_cou2 <- dplyr::group_by(df[chosen, ], .data$clinical) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = cou_from_barcode_list(.data$nums)
        )
      summary_cou3 <- dplyr::group_by(df[chosen, ]) %>%
        dplyr::summarise(
          N = sum(!is.na(.data$nums)),
          age = mean(.data$age[!is.na(.data$nums)]),
          mean = cou_from_barcode_list(.data$nums)
        )

      cou_list <- list(summary_cou1, summary_cou2, summary_cou3)
      summary_cou <- dplyr::bind_rows(cou_list)


      # bundle it all up
      res <- list(
        summary_coi = summary_coi,
        summary_coi_detected = summary_coi_detected,
        summary_moi_detected = summary_moi_detected,
        summary_polygenom = summary_polygenom,
        summary_ages = summary_ages,
        summary_unique = summary_unique,
        summary_cou = summary_cou,
        clonality = clonality,
        coi_table = table(df$coi[chosen]),
        barcode_freq = barcode_freq
      )

      if (barcodes) {
        res$barcodes <- barcodes_tab[barcodes_tab > 1]
      }

      results[[i]] <- res
    }
  } else {
    # set up results
    results <- list()
    length(results) <- reps

    # make samples and reps of them
    for (i in seq_len(length(results))) {
      # make a sample
      if (n[1] == Inf) {
        chosen <- seq_len(nrow_df)
      } else {
        chosen <- ranges[[i]]
      }

      # create summary data frame
      summary_within_ibd <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "pibd_within",
        groupvars = groupvars,
        mean_only = mean_only
      )

      summary_ages <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "age",
        groupvars = groupvars,
        mean_only = mean_only
      )

      summary_ibd <- summarySE_mean_only(
        df[chosen, ],
        measurevar = "pibd",
        groupvars = groupvars,
        mean_only = mean_only
      )

      res <- list(
        summary_ibd = summary_ibd,
        summary_within_ibd = summary_within_ibd,
        summary_ages = summary_ages,
        mean_ibd = mean(df$pibd[chosen], na.rm = TRUE),
        mean_ibd_d = mean(df$pibd_d[chosen], na.rm = TRUE),
        mean_within = mean(df$pibd_within[chosen], na.rm = TRUE),
        mean_within_d = mean(df$pibd_within_d[chosen], na.rm = TRUE)
      )

      results[[i]] <- res
    }
  }

  # if no reps then unlist
  if (reps == 1) {
    results <- res
  }
  return(results)
}
