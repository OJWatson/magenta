#---
#' Plot ordered heatmap
#'
#' \code{heatmap_ordered_binary_plot} creates a heatmap of all human barcodes,
#' with the most common barcode grouped at the bottom
#'
#' @param sim_save Saved output of simulation
#' @param years Lenth of simulation
#' @param EIR Numeric for simulation EIR
#' @param ordered Boolean for oredring heatmap. Default = TRUE
#' @param save_path File path where tiff is saved to. Default = NULL
#'
#' \code{heatmap_ordered_binary_plot}
#'
#' @keywords internal
heatmap_ordered_binary_plot <- function(sim_save, years, EIR, ordered = TRUE, save_path = NULL) {
  if ("populations_event_and_strains_List" %in% names(sim_save)) {
    mat <-
      matrix(as.numeric(
        unlist(
          sim_save$populations_event_and_strains_List$Strain_barcode_vectors
        )
      ), ncol = 24, byrow = T)
  } else {
    mat <-
      matrix(as.numeric(unlist(sim_save$Strain_barcode_vectors)),
        ncol = 24, byrow =
          T
      )
  }

  if (ordered) {
    ID <- apply(mat, MARGIN = 1, FUN = bitsToInt)
    tabled <- sort(table(ID), decreasing = TRUE)
    sorted_row_pos <-
      unlist(sapply(names(tabled), function(x) {
        return(which(ID == as.numeric(x)))
      }))

    if (is.null(save_path)) {
      heatmap(
        mat[sorted_row_pos, ],
        Rowv = NA,
        Colv = NA,
        main = paste0(years, " years | EIR = ", EIR),
        scale = "none",
        labRow = ""
      )
    } else {
      tiff(
        save_path,
        width = 620,
        height = 620,
        units = "mm",
        res = 300,
        pointsize = 36,
        compression = "lzw",
        family = "Times New Roman"
      )
      heatmap(
        mat[sorted_row_pos, ],
        Rowv = NA,
        Colv = NA,
        main = paste0(years, " years | EIR = ", EIR),
        scale = "none",
        labRow = ""
      )
      dev.off()
    }
  } else {
    if (is.null(save_path)) {
      heatmap(
        mat,
        Rowv = NA,
        Colv = NA,
        main = paste0(years, " years | EIR = ", EIR),
        scale = "none",
        labRow = ""
      )
    } else {
      tiff(
        save_path,
        width = 620,
        height = 620,
        units = "mm",
        res = 300,
        pointsize = 36,
        compression = "lzw",
        family = "Times New Roman"
      )
      heatmap(
        mat,
        Rowv = NA,
        Colv = NA,
        main = paste0(years, " years | EIR = ", EIR),
        scale = "none",
        labRow = ""
      )
      dev.off()
    }
  }
}


#---
#' Convert human barcodes to numerics
#'
#' \code{convert_barcode_vectors} converts human barcode vectors to nums and
#'   calculated COIs.
#'
#' @param sim_save Saved output of simulation
#' @param ID Vector of detection immunities for the population
#' @param sub_patents_included Boolean as to whether subpatents are included.
#'   Default = TRUE
#' @param ibd Boolean for IBD simulations. Default = FALSE
#' @param nl Numeric for number of loci. Default = 24
#' @param COI_type String for type of COI calculation. Default = "pcr_imperial"
#'
#' \code{convert_barcode_vectors}
#'
#' @export


convert_barcode_vectors <- function(sim_save, ID,
                                    sub_patents_included = TRUE,
                                    ibd = FALSE, nl = 24,
                                    COI_type = "pcr_imperial") {
  if (!(COI_type %in% c("pcr_imperial", "pcr_alternative", "patent"))) {
    stop("COI_type is not correct. Must be one of pcr_imperial, pcr_alternative, patent")
  }

  # grab their ages so we can work out their detection immunity
  ages <- sim_save$Ages
  ages[ages == 0] <- 0.001
  mpl <- model_param_list_create()
  micro_det <- q_fun(mpl$d1, ID, mpl$ID0, mpl$kD, sapply(ages, fd, mpl$fD0, mpl$aD, mpl$gammaD))

  # parasites strains
  out <- sim_save$Strain_barcode_vectors
  n.strains <- lapply(out, length) %>% unlist()

  # covnvert their barcodes into integers
  if (ibd) {
    intout <- population_ibd_barcodes(out, nl = nl)
  } else {
    intout <- lapply(out, function(x) {
      if (length(x) > 0) {
        lapply(x, function(y) {
          return(bitsToInt(y))
        })
      } else {
        return(NA)
      }
    })
  }

  int.out <- lapply(lapply(intout, function(x) {
    return(unlist(x))
  }), unlist)

  # Calcuate COI based on PCR detectability defined in Imperial College Model
  COI_pcr_imperial <- function(i) {
    a <- length(unique(int.out[[i]][which(sim_save$Strain_infection_state_vectors[[i]] == 1)]))
    a <- a + round(length(unique(int.out[[i]][which(sim_save$Strain_infection_state_vectors[[i]] == 2)])) * (micro_det[i]^mpl$alphaA))
    a <- a + round(length(unique(int.out[[i]][which(sim_save$Strain_infection_state_vectors[[i]] == 3)])) * (micro_det[i]^mpl$alphaU))
    a <- a + length(unique(int.out[[i]][which(sim_save$Strain_infection_state_vectors[[i]] == 4)]))
    return(a)
  }

  # Calculate COI based on alternative PCR detectability method
  COI_pcr_alternative <- function(i) {
    st <- sim_save$Strain_infection_state_vectors[[i]]
    stch <- sim_save$Strain_day_of_infection_state_change_vectors[[i]]
    staying <- st == 1 | st == 4 | st == 2
    staying[st == 3] <- as.logical(rbinom(sum(st == 3), size = 1, micro_det[i]^mpl$alphaU))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }

  # Calculate COI based on if parasites are microscopically detected
  COI_patent <- function(i) {
    st <- sim_save$Strain_infection_state_vectors[[i]]
    stch <- sim_save$Strain_day_of_infection_state_change_vectors[[i]]
    staying <- st == 1 | st == 4
    staying[st == 2] <- as.logical(rbinom(sum(st == 2), size = 1, micro_det[i]))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }

  # calculate the COIs
  if (!ibd) {
    if (!sub_patents_included) {
      COI <- rep(0, length(int.out))
      for (i in 1:length(int.out)) {
        if (n.strains[i] != 0) {
          if (COI_type == "pcr_imperial") {
            COI[i] <- COI_pcr_imperial(i)
          } else if (COI_type == "pcr_alternative") {
            COI[i] <- COI_pcr_alternative(i)
          } else {
            COI[i] <- COI_patent(i)
          }
        }
      }
    } else {
      COI <- unlist(lapply(int.out, function(x) {
        return(length(unique(x)))
      }))
    }
  } else {
    COI <- rep(1, length(int.out))
  }

  # Change NAs to 0
  COI[which(n.strains == 0)] <- 0

  return(list("nums" = int.out, "COI" = COI))
}

#---
#' Sample x times from saved simulation and produce age COI data
#'
#' \code{sample_coi} samples from a human save according to some age distirbution
#'
#' @param sim_save Saved output of simulation
#' @param sample_size Numeric for sample size
#' @param age_densities Densitvy vector for age distribution required
#' @param age_breaks Corresponding vector of age breaks for the density
#' @param reps How many samples are made
#' @param COI_type type of COI calculation
#' @inheritParams convert_barcode_vectors
#'
#' @keywords internal
sample_coi <- function(sim_save, ID, sample_size, age_densities, age_breaks = seq(0, 90 * 365, 2 * 365),
                       sub_patents_included = FALSE, reps = 50, COI_type = "pcr_imperial") {
  COI_out <- convert_barcode_vectors(sim_save, ID,
    sub_patents_included = sub_patents_included,
    COI_type = COI_type
  )

  grouped_ages_by_2 <- cut(sim_save$Ages, breaks = age_breaks, labels = 1:45)
  sample_groups <- round(sample_size * (age_densities * (1 / sum(age_densities))))
  picks <- matrix(rep(rep(1:45, sample_groups), reps), ncol = reps)

  ids <- picks

  for (i in 1:reps) {
    tabled <- table(picks[, i])
    id <- c()

    for (j in names(tabled)) {
      x <- which(grouped_ages_by_2 == j & sim_save$Infection_States %in% c(1, 2, 3, 4))
      id <- c(id, sample(x = x, size = tabled[j], replace = T))
    }

    ids[, i] <- id
  }

  return(list("ids" = ids, "COI" = COI_out$COI, "sim_save" = sim_save))
}

#---
#' Plot COI against age from random sample
#'
#' \code{sample_coi} samples from a human save according to some age distirbution
#'
#' @param sample_coi_out Output of sample_coi
#' @param x Which sample to plot
#' @param span Smoothing parameter for loess
#' @param ylimmax ylim max
#' @param xlimmax xlim max
#' @param max_coi max_coi
#' @importFrom ggplot2 ggplot aes geom_smooth theme_bw geom_point xlim ylim
#'
#' @keywords internal
coi_age_plot_sample_x <- function(sample_coi_out,
                                  x,
                                  span = 0.6,
                                  ylimmax = NULL,
                                  xlimmax = NULL,
                                  max_coi = 25) {
  mround <- function(x, base) {
    base * round(x / base)
  }

  ids <- sample_coi_out$ids[, x]
  Ages <- sample_coi_out$sim_save$Ages[ids] / 365
  COI <- sample_coi_out$COI[ids]
  df <- data.frame("Ages" = Ages, "COI" = COI)

  df$COI[df$COI > max_coi] <- sample(df$COI[df$COI < max_coi], sum(df$COI > max_coi))
  df$Ages <- mround(df$Ages, 1)

  if (is.null(ylimmax)) ylimmax <- max(df$COI) + 1
  if (is.null(xlimmax)) xlimmax <- max(df$Ages) + 1


  gg <- ggplot(df, aes(Ages, COI)) +
    geom_point(alpha = 0.3, size = 2, colour = "blue", fill = "blue", shape = 16) +
    geom_smooth(span = span, color = "red", se = F) +
    theme_bw(base_size = 22) +
    ylim(c(0, ylimmax)) +
    xlim(c(0, xlimmax))

  invisible(gg)
}


#---
#' Create summary prev and genetic dataframes
#'
#' @param res_list list of simulation outputs
#' @param update_length Days between sim_saves in simulation. Default = 30
#' @param years Years sims were run for. Default = 35
#' @param age_breaks Character vecotr for the age break levels used when
#'   sumarising the data by age. Default = c("(-0.001,5]","(5,15]","(15,100]")
#' @param coi_detect Boolean whether to report the COI detected, i.e. what we
#'   believe would be detected my sequencing as sensitive as pcr
#'   Default = FALSE, i.e the true COI
#' @param ibd Were we collecting info on IBD
#'
#' @keywords internal
summary_data_frames_from_sims <- function(res_list,
                                          update_length = 30,
                                          years = 35,
                                          age_breaks = c("(-0.001,5]", "(5,15]", "(15,100]"),
                                          coi_detect = FALSE,
                                          ibd = FALSE) {
  if (coi_detect) {
    coi_det_var <- "summary_coi_detected"
  } else {
    coi_det_var <- "summary_coi"
  }

  age_breaks_all <- c("(-0.001,5]", "(5,15]", "(15,100]")
  states_all <- c("Clinical", "Asymptomatic")
  # set up list for final saves as these need to be replaced by COIs
  reps <- length(res_list)

  # grab the final time point and convert them to the rest of the time points
  time_steps <- length(res_list[[1]]) - 1

  # create vector of correspondng times
  full_time <- seq(30, by = 30, to = time_steps * update_length)

  # positions of these times
  full_time_length <- length(full_time)
  positions <- 1:full_time_length
  time_vec <- sort(rep(positions, length(res_list)))

  # summary functions
  # ----------------------------------------------------------------------------
  mean_a_at_states_and_ages_in_c <- function(res_list, a,
                                             states = c("Asymptomatic", "Clinical"),
                                             ages = NULL,
                                             c, i) {
    if (is.null(ages)) {
      ages <- age_breaks_all
    } else if (is.na(ages)) {} else {
      ages <- ages
    }
    if (is.null(states)) {
      states <- states_all
    }
    lapply(res_list, function(x) {
      infs <- c(x[[positions[i]]][[c]]$clinical %in% states)
      age_pos <- c(x[[positions[i]]][[c]]$age_bin %in% ages)
      res <- weighted.mean(x[[positions[i]]][[c]][[a]][infs & age_pos],
        x[[positions[i]]][[c]]$N[infs & age_pos],
        na.rm = TRUE
      )
      return(res)
    }) %>% unlist()
  }

  mean_prev <- function(res_list, full_time_length, prev = "pcr_prev") {
    i <- 1:(full_time_length)
    lapply(i, function(y) {
      lapply(res_list, function(x) {
        ((x[[positions[y]]][[prev]]))
      })
    }) %>% unlist()
  }

  mean_over_time <- function(res_list, a,
                             states = c("Asymptomatic", "Clinical"),
                             age_breaks = NULL,
                             c,
                             full_time_length) {
    i <- 1:full_time_length
    if (is.null(age_breaks) || is.na(age_breaks)) {
      return(unlist(lapply(i, function(y) {
        mean_a_at_states_and_ages_in_c(res_list,
          a = a,
          c = c, states = states,
          ages = age_breaks,
          i = y
        )
      })))
    } else {
      return(lapply(age_breaks, function(x) {
        unlist(lapply(i, function(y) {
          mean_a_at_states_and_ages_in_c(res_list,
            a = a,
            c = c, states = states,
            ages = x,
            i = y
          )
        }))
      }))
    }
  }

  # ----------------------------------------------------------------------------


  list_res <- list()

  if (!ibd) {
    # build our results for the time steps
    rept <- rep(1:length(res_list), length(full_time))

    list_res <- list()

    list_res$mean_polygenomic <- mean_over_time(res_list, a = "mean", c = "summary_polygenom", full_time_length = full_time_length)
    list_res$mean_polygenomic_ages <- mean_over_time(res_list, a = "mean", c = "summary_polygenom", full_time_length = full_time_length, age_breaks = age_breaks)
    list_res$mean_polygenomic_Clinical <- mean_over_time(res_list, a = "mean", c = "summary_polygenom", full_time_length = full_time_length, states = c("Clinical"))
    list_res$mean_polygenomic_Asymptomatic <- mean_over_time(res_list, a = "mean", c = "summary_polygenom", full_time_length = full_time_length, states = c("Asymptomatic"))

    list_res$meancoi <- mean_over_time(res_list, a = "mean", c = coi_det_var, full_time_length = full_time_length)
    list_res$meancoi_ages <- mean_over_time(res_list, a = "mean", c = coi_det_var, full_time_length = full_time_length, age_breaks = age_breaks)
    list_res$meancoi_Clinical <- mean_over_time(res_list, a = "mean", c = coi_det_var, full_time_length = full_time_length, states = c("Clinical"))
    list_res$meancoi_Asymptomatic <- mean_over_time(res_list, a = "mean", c = coi_det_var, full_time_length = full_time_length, states = c("Asymptomatic"))

    list_res$mean_unique <- mean_over_time(res_list, a = "mean", c = "summary_unique", full_time_length = full_time_length, states = NA, age_breaks = NA) %>% unlist()
    list_res$mean_unique_ages <- mean_over_time(res_list, a = "mean", c = "summary_unique", full_time_length = full_time_length, age_breaks = age_breaks, states = NA)
    list_res$mean_unique_Clinical <- mean_over_time(res_list, a = "mean", c = "summary_unique", full_time_length = full_time_length, states = c("Clinical"), age_breaks = NA) %>% unlist()
    list_res$mean_unique_Asymptomatic <- mean_over_time(res_list, a = "mean", c = "summary_unique", full_time_length = full_time_length, states = c("Asymptomatic"), age_breaks = NA) %>% unlist()

    list_res$mean_cou <- mean_over_time(res_list, a = "mean", c = "summary_cou", full_time_length = full_time_length, states = NA, age_breaks = NA) %>% unlist()
    list_res$mean_cou_ages <- mean_over_time(res_list, a = "mean", c = "summary_cou", full_time_length = full_time_length, age_breaks = age_breaks, states = NA)
    list_res$mean_cou_Clinical <- mean_over_time(res_list, a = "mean", c = "summary_cou", full_time_length = full_time_length, states = c("Clinical"), age_breaks = NA) %>% unlist()
    list_res$mean_cou_Asymptomatic <- mean_over_time(res_list, a = "mean", c = "summary_cou", full_time_length = full_time_length, states = c("Asymptomatic"), age_breaks = NA) %>% unlist()

    list_res$mean_prev <- mean_prev(res_list, full_time_length = full_time_length)

    # turn this into a data frame and melt it for plotting
    df <- as.data.frame.list(list_res)
    df$rep <- rept
    df$time <- time_vec
    df$time <- df$time / 12
    interventions <- factor(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
    df$intervention <- interventions[rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), length(df$time) / 30)]

    names(df)[grep("ages", names(df))] <- c(
      paste0("mean_polygenomic_ages_", c("0-5", "5-15", "15+")),
      paste0("mean_coi_ages_", c("0-5", "5-15", "15+")),
      paste0("mean_unique_ages_", c("0-5", "5-15", "15+")),
      paste0("mean_cou_ages_", c("0-5", "5-15", "15+"))
    )

    # renaming
    melted <- reshape2::melt(df, id.vars = c("time", "rep", "intervention"))
    melted$metric <- as.character(melted$variable)
    melted$metric[grepl("coi", melted$metric)] <- "COI"
    melted$metric[grepl("cou", melted$metric)] <- "COU"
    melted$metric[grepl("unique", melted$metric)] <- "% Unique"
    melted$metric[grepl("polygenom", melted$metric)] <- "% Polygenomic"
    melted$sampled <- "All"

    # refactoring and naming
    terms <- c("0-5", "5-15", "15+", "Clinical", "Asymptomatic")
    for (t in terms) {
      melted$sampled[grepl(t, melted$variable, fixed = TRUE)] <- t
    }
    levels(melted$sampled) <- c("All", terms)

    deeped <- dplyr::group_by(melted, .data$intervention, .data$metric, .data$sampled, .data$time) %>%
      dplyr::summarise(value = mean(.data$value, na.rm = TRUE))
    deeped$rep <- deeped$intervention

    return(list("melted" = melted, "deeped" = deeped))
  } else {
    # build our results for the time steps
    rept <- rep(1:length(res_list), length(full_time))

    # summaarise what we need
    list_res$mean_prev <- mean_prev(res_list, full_time_length = full_time_length)

    list_res$mean_ibd <- mean_over_time(res_list, a = "mean", c = "summary_ibd", full_time_length = full_time_length, states = NULL, age_breaks = NULL) %>% unlist()
    list_res$mean_ibd_ages <- mean_over_time(res_list, a = "mean", c = "summary_ibd", full_time_length = full_time_length, age_breaks = age_breaks, states = NULL)
    list_res$mean_ibd_Clinical <- mean_over_time(res_list, a = "mean", c = "summary_ibd", full_time_length = full_time_length, states = c("Clinical"), age_breaks = NULL) %>% unlist()
    list_res$mean_ibd_Asymptomatic <- mean_over_time(res_list, a = "mean", c = "summary_ibd", full_time_length = full_time_length, states = c("Asymptomatic"), age_breaks = NULL) %>% unlist()

    list_res$mean_within <- mean_over_time(res_list, a = "mean", c = "summary_within_ibd", full_time_length = full_time_length, states = NULL, age_breaks = NULL) %>% unlist()
    list_res$mean_within_ages <- mean_over_time(res_list, a = "mean", c = "summary_within_ibd", full_time_length = full_time_length, age_breaks = age_breaks, states = NULL)
    list_res$mean_within_Clinical <- mean_over_time(res_list, a = "mean", c = "summary_within_ibd", full_time_length = full_time_length, states = c("Clinical"), age_breaks = NULL) %>% unlist()
    list_res$mean_within_Asymptomatic <- mean_over_time(res_list, a = "mean", c = "summary_within_ibd", full_time_length = full_time_length, states = c("Asymptomatic"), age_breaks = NULL) %>% unlist()

    # turn this into a data frame and melt it for plotting
    df <- as.data.frame.list(list_res)
    df$rep <- rept
    df$time <- time_vec
    df$time <- df$time / 12
    interventions <- factor(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
    df$intervention <- interventions[rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), length(df$time) / 30)]

    # renaming
    names(df)[grep("ages", names(df))] <- c(
      paste0("mean_ibd_ages_", c("0-5", "5-15", "15+")),
      paste0("mean_within_ages_", c("0-5", "5-15", "15+"))
    )

    melted <- reshape2::melt(df, id.vars = c("time", "rep", "intervention"))
    melted$metric <- as.character(melted$variable)
    melted$metric[grepl("ibd", melted$metric)] <- "pIBD"
    melted$metric[grepl("within", melted$metric)] <- "iIBD"
    melted$sampled <- "All"
    terms <- c("0-5", "5-15", "15+", "Clinical", "Asymptomatic")
    for (t in terms) {
      melted$sampled[grepl(t, melted$variable, fixed = TRUE)] <- t
    }
    levels(melted$sampled) <- c("All", terms)

    melted$value[which(melted$value == 0)] <- 1 # convert iIBD measures to 1 if 0 (i.e. all monogenomic)

    # summarise across simulation reps
    deeped <- dplyr::group_by(melted, .data$intervention, .data$metric, .data$sampled, .data$time) %>%
      dplyr::summarise(value = mean(.data$value, na.rm = TRUE))
    deeped$rep <- deeped$intervention

    return(list("melted" = melted, "deeped" = deeped))
  }
}


#' @noRd
plot_prevalence <- function(res, age_bin = NULL, states = c("D", "A", "U", "T")) {
  infs <- c("S", "D", "A", "U", "T", "P")
  inf_states <- which(infs %in% states)

  l <- length(res)
  prev <- rep(0, l - 1)

  if (is.null(age_bin)) {
    age_bin <- c(0, (365 * 100))
  } else {
    age_bin <- age_bin * 365
  }

  if (sum(names(res[[1]]) %in% c("Ages", "InfectionStates")) == 2) {
    n <- length(res[[1]]$InfectionStates)
    for (i in 1:(l - 1)) {
      prev[i] <- sum(res[[i]]$InfectionStates %in% (inf_states - 1) &
        res[[i]]$Ages > age_bin[1] &
        res[[i]]$Ages < age_bin[2]) /
        sum(res[[i]]$Ages > age_bin[1] &
          res[[i]]$Ages < age_bin[2])
    }
  } else {
    for (i in 1:(l - 1)) {
      prev[i] <- res[[i]][states] %>%
        unlist() %>%
        sum()
    }
  }
  plot(prev, type = "l")
  invisible(prev)
}

#' @noRd
plot_x_loggers <- function(res, x, alpha = 1, scale = FALSE, extra = NULL) {
  l <- list()
  length(l) <- length(x)
  names(l) <- x
  want <- 1:length(res)

  for (i in 1:length(x)) {
    l[[i]] <- lapply(res[want], function(y) y[x[i]]) %>% unlist()
  }

  df <- as.data.frame.list(l)
  if (scale) {
    df <- t(apply(df, 1, function(z) z / sum(z))) %>% as.data.frame()
  }
  df$time <- 1:length(l[[1]])
  if ("summary_ibd" %in% names(res[[1]])) {
    df$prev <- lapply(res[1:(length(res) - 1)], function(x) {
      sum(x$summary_ibd$N[x$summary_ibd$state %in% c("D", "A", "T", "U")])
    }) %>% unlist()
    df$prev <- df$prev / sum(res[[1]]$summary_ibd$N)
  } else {
    df$prev <- lapply(res[1:(length(res) - 1)], function(x) {
      sum(x$summary_coi$N[x$summary_coi$state %in% c("D", "A", "T", "U")])
    }) %>% unlist()
    df$prev <- df$prev / sum(res[[1]]$summary_coi$N)
  }
  melt <- reshape2::melt(df, id.var = "time")

  if (is.null(extra)) {
    gg <- ggplot2::ggplot(melt, ggplot2::aes(x = time, y = .data$value, color = .data$variable)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::geom_line(alpha = alpha)
  } else {
    gg <- ggplot2::ggplot(melt, ggplot2::aes(x = time, y = .data$value, color = .data$variable)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::geom_line(alpha = alpha) +
      eval(parse(text = extra))
  }
  print(gg)
  invisible(list("df" = melt, "gg" = gg))
}

#' @noRd
plot_mean_summary <- function(res, x,
                              group = "clinical",
                              groups = c("Asymptomatic", "Clinical"),
                              alpha = 1,
                              scale = FALSE,
                              extra = NULL) {
  l <- list()
  length(l) <- length(x)
  names(l) <- x
  want <- 1:(length(res) - 1)

  for (i in 1:length(x)) {
    l[[i]] <- lapply(res[want], function(y) {
      N <- y[[x]]$N[y[[x]][[group]] %in% groups]
      mean <- y[[x]]$mean[y[[x]][[group]] %in% groups]
      return(sum(N * mean) / sum(N))
    }) %>% unlist()
  }

  df <- as.data.frame.list(l)
  if (scale) {
    df <- t(apply(df, 1, function(z) z / sum(z))) %>% as.data.frame()
  }
  df$time <- 1:length(l[[1]])
  if ("summary_ibd" %in% names(res[[1]])) {
    df$prev <- lapply(res[1:(length(res) - 1)], function(x) {
      sum(x$summary_ibd$N[x$summary_ibd[[group]] %in% groups])
    }) %>% unlist()
    df$prev <- df$prev / sum(res[[1]]$summary_ibd$N)
  } else {
    df$prev <- lapply(res[1:(length(res) - 1)], function(x) {
      sum(x$summary_coi$N[x$summary_coi[[group]] %in% groups])
    }) %>% unlist()
    df$prev <- df$prev / sum(res[[1]]$summary_coi$N)
  }
  melt <- reshape2::melt(df, id.var = "time")

  if (is.null(extra)) {
    gg <- ggplot2::ggplot(melt, ggplot2::aes(x = time, y = .data$value, color = .data$variable)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::geom_line(alpha = alpha)
  } else {
    gg <- ggplot2::ggplot(melt, ggplot2::aes(x = time, y = .data$value, color = .data$variable)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::geom_line(alpha = alpha) +
      eval(parse(text = extra))
  }
  print(gg)
  invisible(list("df" = melt, "gg" = gg))
}

# Plot resistance frequenceies and lineages over time:
#' @noRd
plot_resistance_lineages <- function(t, legend = TRUE, va = FALSE) {
  l <- length(t) - 1
  bs <- length(t[[1]]$af)
  bst <- 2^bs
  colours <- hcl(seq(360 / bst, 360, length.out = bst),
    c = 100, l = 75, alpha = 1, fixup = TRUE
  )

  if ("lineage" %in% names(t[[1]])) {
    lins <- TRUE
  } else {
    lins <- FALSE
  }

  af <- matrix(unlist(lapply(t[1:l], "[[", "af")), ncol = bs, byrow = TRUE)
  muts <- matrix(unlist(lapply(t[1:l], "[[", "mutations")), ncol = bs, byrow = TRUE)

  if (lins) {
    lin <- matrix(unlist(lapply(t[1:l], "[[", "lineage")), ncol = bst, byrow = TRUE)
  }

  af_df <- data.frame("time" = 1:l, "af" = as.numeric(af), "loci" = as.character(sort(rep(1:ncol(af), l))))
  af_gg <- ggplot(af_df, aes(x = .data$time, y = .data$af, color = .data$loci)) +
    geom_point() +
    geom_line() +
    scale_color_discrete(name = "Loci") +
    ylab("AF") +
    theme(axis.title.x = element_blank())

  mut_df <- data.frame("time" = 1:l, "muts" = as.numeric(muts), "loci" = as.character(sort(rep(1:ncol(af), l))))
  mut_gg <- ggplot(mut_df, aes(x = .data$time, y = .data$muts, color = .data$loci)) +
    geom_point() +
    geom_line() +
    scale_color_discrete(name = "Loci") +
    ylab("Mutations") +
    theme(axis.title.x = element_blank())

  if (lins) {
    lin_names <- apply(
      sapply(0:((bst) - 1), function(x) {
        intToBits(x)
      }),
      2,
      function(x) {
        paste0(as.numeric(x[1:bs]), collapse = "")
      }
    )
    vad <- apply(
      sapply(0:((bst) - 1), function(x) {
        intToBits(x)
      }),
      2,
      function(x) {
        as.logical(x[bs])
      }
    )

    lin_df <- data.frame(
      "time" = 1:l, "af" = as.numeric(lin),
      "lin" = factor(as.character(mapply(rep, lin_names, l)), levels = lin_names),
      "vad" = as.logical(mapply(rep, vad, l))
    )

    if (va) {
      va <- c("", "VA")
    } else {
      va <- c("", "")
    }
    va_group <- paste(lin_names, va[as.numeric(vad) + 1])
    lin_df$group <- factor(as.character(mapply(rep, va_group, l)), levels = va_group)

    if (length(colours) <= 16) {
      colours <- rep(RColorBrewer::brewer.pal(8, "Paired"), 2)
    }

    lin_gg <- ggplot(lin_df, aes(x = .data$time, y = .data$af, color = .data$group, shape = .data$group)) +
      geom_point() +
      geom_line() +
      scale_color_manual(name = "Barcode", labels = lin_names, values = colours) +
      ylab("N") +
      theme(axis.title.x = element_blank())
    heights <- c(1, 1, 1, 1, 1)

    if (length(colours) <= 16) {
      lin_gg <- lin_gg +
        scale_shape_manual(name = "Barcode", labels = lin_names, values = c(rep(19, 8), rep(17, 8)))
    } else {
      lin_gg <- lin_gg + theme(legend.position = "none")
    }
  } else {
    lin_gg <- NULL
    heights <- c(1, 1, 0, 1, 1)
  }


  prev <- lapply(t[1:l], function(x) {
    x[c("pcr_prev", "dat", "micro_2_10")]
  }) %>%
    rbind_list_base() %>%
    reshape2::melt() %>%
    dplyr::mutate(time = rep(1:l, 3)) %>%
    ggplot(aes(x = .data$time, y = .data$value, color = .data$variable)) +
    geom_line() +
    scale_color_discrete(name = "Type") +
    ylab("Prevalence") +
    theme(axis.title.x = element_blank())

  ntf <- lapply(tail(t, l / 2), function(x) {
    sum(x$unsuccesful_treatments_lpf, x$not_treated)
  }) %>%
    unlist() %>%
    sum()

  ntf <- ntf / (length(t[[l + 1]]$Loggers$Ages) / 100) / round(t[[l + 1]]$Loggers$Time / (2 * 365))

  tf <- lapply(t[1:l], function(x) {
    x[c("overall_treatment_failure")]
  }) %>%
    rbind_list_base() %>%
    reshape2::melt() %>%
    dplyr::mutate(time = rep(1:l, 1)) %>%
    ggplot(aes(x = .data$time, y = .data$value, color = .data$variable)) +
    geom_line() +
    ylab("Treatment Failure") +
    theme(axis.title.x = element_blank(), legend.position = "none")

  if (legend) {
    cowplot::plot_grid(af_gg, mut_gg,
      lin_gg, tf,
      prev,
      rel_heights = heights, ncol = 1
    )
  } else {
    if (length(colours) <= 16) {
      ling_gg <- lin_gg + theme(legend.position = "top") +
        guides(col = guide_legend(ncol = bst))
    }
    cowplot::plot_grid(af_gg + theme(legend.position = "top"),
      mut_gg + theme(legend.position = "top"),
      lin_gg,
      tf + theme(legend.position = "top"),
      prev + theme(legend.position = "top"),
      rel_heights = heights, ncol = 1
    )
  }
}


# function to calculate probability of strain being detected by microscopy
q_fun <- function(d1, ID, ID0, kD, fd) {
  return(d1 + ((1 - d1) / (1 + ((ID / ID0)^kD) * fd)))
}

fd <- function(age, fD0, aD, gammaD) {
  return(1 - ((1 - fD0) / (1 + (age / aD)^gammaD)))
}


#' @noRd
microscopy_detected <- function(ages, IDs, states, age_brackets = c(2 * 365, 10 * 365)) {
  age_in <- ages < age_brackets[2] & ages > age_brackets(1)
  dt <- sum(states %in% c(1, 4) & age_in) / sum(age_in)
  mpl <- model_param_list_create()

  micro_det <- q_fun(
    mpl$d1, IDs, mpl$ID0, mpl$kD,
    sapply(ages, fd, mpl$fD0, mpl$aD, mpl$gammaD)
  )

  a <- sum(states %in% 2 & age_in) / sum(age_in) * mean(micro_det[age_in])

  return(dt + a)
}

#' Function to generate ggplot colours
#'
#' @param n Number of colors
#' @keywords internal
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
