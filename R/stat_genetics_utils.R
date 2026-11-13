#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars Characters of grouping variables
#' @param conf.interval numeric for CI
#'
#' @keywords internal
summarySE <- function(data = NULL, measurevar = "COI",
                      groupvars = c("Age_Bin", "State"),
                      conf.interval = .95) {
  res <- dplyr::group_by_at(data, dplyr::vars(dplyr::one_of(groupvars))) %>%
    dplyr::summarise(
      N = length(!!dplyr::sym(measurevar)),
      mean = mean(!!dplyr::sym(measurevar), na.rm = TRUE),
      sd = sd(!!dplyr::sym(measurevar), na.rm = TRUE),
      se = sd / sqrt(.data$N),
      ci = .data$se * suppressWarnings(qt(conf.interval / 2 + .5, .data$N - 1))
    )

  return(as.data.frame(res))
}

#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars What are we summarise by
#' @param mean_only Are we just calculating the mean. Default = TRUE
#' @keywords internal
summarySE_mean_only <- function(data = NULL, measurevar = "COI",
                                groupvars = c("Age_Bin", "State"), mean_only = TRUE) {
  if (!mean_only) {
    return(summarySE(data, measurevar, groupvars))
  } else {
    var <- dplyr::sym(measurevar)

    res <- dplyr::group_by_at(data, groupvars) %>%
      dplyr::summarise(
        N = sum(!is.na(!!var)),
        mean = mean(!!var, na.rm = TRUE)
      )

    return(as.data.frame(res))
  }
}

#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars What are we summarise by
#' @param mean_only Are we just calculating the mean. Default = TRUE
#' @param max Maximum value that can be obsered in the data beign summarised.
#'   If greater than zero, all values in data will be set to the max before
#'   summarising. Default = 0.
#' @keywords internal
summarySE_mean_only_max_mean <- function(data = NULL,
                                         measurevar = "COI",
                                         groupvars = c("Age_Bin", "State"),
                                         mean_only = TRUE,
                                         max = 0) {
  if (nrow(data) != 0) {
    if (max) {
      data[measurevar][data[measurevar] > max] <- max
    }
  }
  if (!mean_only) {
    return(summarySE(data, measurevar, groupvars))
  } else {
    var <- dplyr::sym(measurevar)

    res <- dplyr::group_by_at(data, groupvars) %>%
      dplyr::summarise(
        N = length(!!var),
        mean = mean(!!var, na.rm = TRUE)
      )

    return(as.data.frame(res))
  }
}

#' Zero truncated geometric integers
#'
#' Draws geometric numbers with the given mean, scales them by \code{p} and
#' keeps those above zero, redrawing until \code{n} are in hand. Used for the
#' number of co-transmitted strains per infectious bite.
#'
#' @param n Number of values wanted.
#' @param mean Mean of the geometric draw.
#' @param p Proportion of the draw that survives.
#' @return An integer vector of length \code{n}.
#' @keywords internal
ztrgeomintp <- function(n, mean, p) {
  ngs <- rgeom(n * 1.5, prob = 1 / mean)
  ngs <- round(ngs * p)
  ngs <- ngs[ngs > 0]
  while (length(ngs) < n) {
    ng2 <- rgeom(n * 1.5, prob = 1 / mean)
    ngs <- c(ngs, round(ng2 * p))
    ngs <- ngs[ngs > 0]
  }

  ngs <- sample(ngs, size = n, replace = FALSE)
  return(ngs)
}

#' Zero truncated negative binomial integers
#'
#' Used for the number of oocysts formed from a bite.
#'
#' @param n Number of values wanted.
#' @param mean Mean of the negative binomial.
#' @param size Shape parameter.
#' @return An integer vector of length \code{n}.
#' @keywords internal
ztrnbinom <- function(n, mean, size) {
  nbs <- rnbinom(n * 1.2, size = size, mu = mean)
  nbs <- nbs[nbs > 0]
  while (length(nbs) < n) {
    nbs <- c(nbs, rnbinom(n * .5, size = size, mu = mean))
    nbs <- nbs[nbs > 0]
  }

  nbs <- sample(nbs, size = n, replace = FALSE)
  return(nbs)
}

#' Clonality from a list of barcode numbers
#'
#' The proportion of distinct barcodes seen in only one individual.
#'
#' @param barcode_list List of integer barcodes, one element per individual.
#' @return A single number, \code{0} when every barcode is shared.
#' @keywords internal
clonality_from_barcode_list <- function(barcode_list) {
  tbl <- table(table(unlist(lapply(barcode_list, unique))))

  if ("1" %in% names(tbl)) {
    return(tbl[1] / sum(tbl * as.numeric(names(tbl))))
  } else {
    return(0)
  }
}

#' Complexity of uniqueness from a list of barcode numbers
#'
#' @param barcode_list List of integer barcodes, one element per individual.
#' @return A single number between 0 and 1.
#' @keywords internal
cou_from_barcode_list <- function(barcode_list) {
  tbl <- table(unlist(lapply(barcode_list, unique)))
  samp_size <- sum(tbl)
  z <- sum((tbl / samp_size)^2)
  return((z - (1 / samp_size)) / (1 - (1 / samp_size)))
}


#' Convert an IBD barcode to one integer per locus
#'
#' In an IBD simulation each locus is a bitset of \code{ibd_length} bits, so
#' the barcode is \code{nl} such bitsets end to end.
#'
#' @param b The barcode as a logical vector.
#' @param nl Number of loci.
#' @return A numeric vector of \code{nl} identities.
#' @keywords internal
convert_ibd_barcode <- function(b, nl) {
  ib <- length(b) / nl

  br <- seq_len(ib)
  ibd <- rep(0, nl)
  for (i in seq_len(nl)) {
    ibd[i] <- bitsToInt(b[br + (ib * (i - 1))])
  }

  return(ibd)
}

#' Update saved state barcode vectors to match a target PLAF
#'
#' This function adjusts the barcode vectors stored in a saved simulation
#' state so that the population level allele frequency (PLAF) matches a user
#' supplied target. The adjustment is applied to each of the supported barcode
#' containers individually (human strain barcodes, infection realisations, and
#' male and female mosquito oocyst barcodes) as well as to the combined set of
#' barcodes across all containers.
#'
#' @param saved_state A list representing the saved simulation state returned
#'   by [saveRDS()] on the final element of a `run_sim()` output.
#' @param target_plaf Numeric vector giving the desired PLAF for each locus.
#' @param seed Optional integer to set the random number generator state before
#'   modifying barcodes. This is useful to obtain reproducible adjustments.
#'
#' @return The input `saved_state` with barcode vectors modified in place so
#'   that their PLAF equals `target_plaf` (within numerical tolerance).
#'
#' @export
#' @examples
#' \donttest{
#' state <- pipeline(
#'   EIR = 1, years = 1, N = 100, num_loci = 2, set_up_only = TRUE, seed = 1L,
#'   housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
#' )
#' state <- update_saved_state_barcode_plaf(state, target_plaf = c(0.9, 0.1), seed = 1L)
#' }

update_saved_state_barcode_plaf <- function(saved_state, target_plaf, seed = NULL) {
  if (!is.list(saved_state)) {
    stop("`saved_state` must be a list produced by readRDS on a saved simulation state.")
  }

  validate_plaf_vector(target_plaf)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  list_paths <- list(
    c("populations_event_and_strains_List", "Strain_barcode_vectors"),
    c("populations_event_and_strains_List", "Infection_barcode_realisation_vectors"),
    c("scourge_List", "Mosquito_Oocyst_barcode_male_vectors"),
    c("scourge_List", "Mosquito_Oocyst_barcode_female_vectors")
  )

  adjusted_matrices <- list()

  for (path in list_paths) {
    if (!has_nested_element(saved_state, path)) {
      stop(sprintf("Saved state does not contain `%s`.", paste(path, collapse = "$")))
    }

    barcode_list <- get_nested_element(saved_state, path)
    flattened <- flatten_barcode_vectors(barcode_list)

    if (length(flattened$flat) == 0) {
      stop(sprintf(
        "Cannot adjust PLAF for `%s` because it does not contain any barcode vectors.",
        paste(path, collapse = "$")
      ))
    }

    barcode_length <- length(flattened$flat[[1]])
    if (length(target_plaf) != barcode_length) {
      stop(sprintf(
        "`target_plaf` has length %d but barcodes in `%s` have length %d.",
        length(target_plaf), paste(path, collapse = "$"), barcode_length
      ))
    }

    barcode_matrix <- flat_to_matrix(flattened$flat, barcode_length)
    adjusted_matrix <- adjust_matrix_to_plaf(barcode_matrix, target_plaf)

    saved_state <- set_nested_element(
      saved_state,
      path,
      rebuild_barcode_vectors(adjusted_matrix, flattened$counts)
    )

    adjusted_matrices[[paste(path, collapse = "$")]] <- adjusted_matrix
  }

  invisible(saved_state)
}

#' Check a target allele frequency vector
#'
#' @param target_plaf The vector.
#' @return \code{target_plaf}, invisibly.
#' @keywords internal
validate_plaf_vector <- function(target_plaf) {
  if (!is.numeric(target_plaf)) {
    stop("`target_plaf` must be a numeric vector.")
  }
  if (length(target_plaf) == 0) {
    stop("`target_plaf` must have positive length.")
  }
  if (any(target_plaf < 0 | target_plaf > 1)) {
    stop("`target_plaf` values must lie between 0 and 1.")
  }
  invisible(target_plaf)
}

#' Does a nested list hold an element at a path?
#'
#' @param x The list.
#' @param path Character vector of names, outermost first.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
has_nested_element <- function(x, path) {
  current <- x
  for (name in path) {
    if (!is.list(current) || is.null(current[[name]])) {
      return(FALSE)
    }
    current <- current[[name]]
  }
  TRUE
}

#' Read an element of a nested list by path
#'
#' @param x The list.
#' @param path Character vector of names, outermost first.
#' @return The element.
#' @keywords internal
get_nested_element <- function(x, path) {
  current <- x
  for (name in path) {
    current <- current[[name]]
  }
  current
}

#' Replace an element of a nested list by path
#'
#' @param x The list.
#' @param path Character vector of names, outermost first.
#' @param value The new element.
#' @return The list with the element replaced.
#' @keywords internal
set_nested_element <- function(x, path, value) {
  if (length(path) == 1) {
    x[[path[[1]]]] <- value
    return(x)
  }

  head_path <- path[[1]]
  x[[head_path]] <- set_nested_element(x[[head_path]], path[-1], value)
  x
}

#' Flatten a list of barcode vectors
#'
#' Turns the per-individual list of barcodes into one flat list, remembering
#' how many each individual had so that
#' \code{\link{rebuild_barcode_vectors}} can undo it.
#'
#' @param barcode_list List with one element per individual, each a list of
#'   logical barcodes (or \code{NULL}).
#' @return List with \code{flat}, the barcodes in order, and \code{counts},
#'   the number per individual.
#' @keywords internal
flatten_barcode_vectors <- function(barcode_list) {
  counts <- integer(length(barcode_list))
  flat <- list()
  flat_idx <- 0L

  for (i in seq_along(barcode_list)) {
    entry <- barcode_list[[i]]
    if (is.null(entry)) {
      entry <- list()
    }

    if (!is.list(entry)) {
      entry <- list(entry)
    }

    counts[i] <- length(entry)

    if (counts[i] > 0) {
      for (j in seq_len(counts[i])) {
        flat_idx <- flat_idx + 1L
        flat[[flat_idx]] <- as.logical(entry[[j]])
      }
    }
  }

  list(flat = flat, counts = counts)
}

#' Stack barcodes into an integer matrix
#'
#' @param flat List of logical barcodes of equal length.
#' @param barcode_length Their length.
#' @return An integer matrix with one row per barcode.
#' @keywords internal
flat_to_matrix <- function(flat, barcode_length) {
  if (length(flat) == 0) {
    return(matrix(NA_real_, nrow = 0, ncol = barcode_length))
  }

  mat <- t(vapply(flat, function(x) {
    if (length(x) != barcode_length) {
      stop("All barcode vectors must have the same length.")
    }
    as.integer(x)
  }, integer(barcode_length)))

  storage.mode(mat) <- "integer"
  mat
}

#' Rebuild the per-individual barcode list from a matrix
#'
#' The inverse of \code{\link{flatten_barcode_vectors}} and
#' \code{\link{flat_to_matrix}}.
#'
#' @param barcode_matrix Integer matrix with one row per barcode.
#' @param counts Number of barcodes per individual.
#' @return List with one element per individual, each a list of logical
#'   barcodes.
#' @keywords internal
rebuild_barcode_vectors <- function(barcode_matrix, counts) {
  barcode_matrix <- as.matrix(barcode_matrix)
  storage.mode(barcode_matrix) <- "integer"
  idx <- 0L
  res <- vector("list", length(counts))

  for (i in seq_along(counts)) {
    n <- counts[i]
    if (n == 0L) {
      res[[i]] <- list()
    } else {
      res_i <- vector("list", n)
      for (j in seq_len(n)) {
        idx <- idx + 1L
        res_i[[j]] <- as.logical(barcode_matrix[idx, ])
      }
      res[[i]] <- res_i
    }
  }

  res
}

#' Flip alleles so that a barcode matrix matches a target frequency
#'
#' For each locus, flips as many randomly chosen entries as needed for the
#' column's allele count to equal the target frequency times the number of
#' barcodes, then checks the result.
#'
#' @param barcode_matrix Integer matrix with one row per barcode.
#' @param target_plaf Target frequency per locus.
#' @return The adjusted matrix.
#' @keywords internal
adjust_matrix_to_plaf <- function(barcode_matrix, target_plaf) {
  n <- nrow(barcode_matrix)
  if (n == 0) {
    return(barcode_matrix)
  }

  barcode_matrix <- as.matrix(barcode_matrix)
  storage.mode(barcode_matrix) <- "integer"

  for (col_idx in seq_along(target_plaf)) {
    target_count <- round(target_plaf[col_idx] * n)

    if (target_count < 0 || target_count > n) {
      stop("`target_plaf` results in impossible allele counts.")
    }

    current_count <- sum(barcode_matrix[, col_idx])
    diff <- target_count - current_count

    if (diff > 0) {
      zero_idx <- which(barcode_matrix[, col_idx] == 0L)
      if (length(zero_idx) < diff) {
        stop("Not enough 0 entries available to reach the target PLAF.")
      }
      if (length(zero_idx) == 1L && diff == 1L) {
        barcode_matrix[zero_idx, col_idx] <- 1L
      } else {
        to_flip <- sample(zero_idx, diff)
        barcode_matrix[to_flip, col_idx] <- 1L
      }
    } else if (diff < 0) {
      one_idx <- which(barcode_matrix[, col_idx] == 1L)
      if (length(one_idx) < abs(diff)) {
        stop("Not enough 1 entries available to reach the target PLAF.")
      }
      if (length(one_idx) == 1L && diff == -1) {
        barcode_matrix[one_idx, col_idx] <- 0L
      } else {
        to_flip <- sample(one_idx, abs(diff))
        barcode_matrix[to_flip, col_idx] <- 0L
      }
    }
  }

  validate_plaf(barcode_matrix, target_plaf, "barcode vectors")
  barcode_matrix
}

#' Check that a barcode matrix matches a target frequency
#'
#' @param barcode_matrix Integer matrix with one row per barcode.
#' @param target_plaf Target frequency per locus.
#' @param context What the matrix holds, for the error message.
#' @return \code{NULL}, invisibly; an error reports the counts expected and
#'   found.
#' @keywords internal
validate_plaf <- function(barcode_matrix, target_plaf, context) {
  if (nrow(barcode_matrix) == 0) {
    return(invisible(NULL))
  }

  if (!identical(round(target_plaf * nrow(barcode_matrix)), colSums(barcode_matrix))) {
    stop(sprintf(
      "Failed to match target PLAF for %s. Expected %s but obtained %s.",
      context,
      paste(
        paste(round(target_plaf * nrow(barcode_matrix)), collapse = ", "),
        "out of", nrow(barcode_matrix)
      ),
      paste(colSums(barcode_matrix), collapse = ", ")
    ))
  }

  invisible(NULL)
}

#' Probability of a strain being detected by microscopy
#'
#' @param d1 Minimum detection probability.
#' @param ID Detection immunity of each individual.
#' @param ID0 Immunity scale.
#' @param kD Immunity shape.
#' @param fd Age effect from \code{\link{fd}}.
#' @return A probability per individual.
#' @keywords internal
q_fun <- function(d1, ID, ID0, kD, fd) {
  return(d1 + ((1 - d1) / (1 + ((ID / ID0)^kD) * fd)))
}

#' Age effect on microscopy detection
#'
#' @param age Age in days.
#' @param fD0 Effect at birth.
#' @param aD Age scale.
#' @param gammaD Age shape.
#' @return A value per age.
#' @keywords internal
fd <- function(age, fD0, aD, gammaD) {
  return(1 - ((1 - fD0) / (1 + (age / aD)^gammaD)))
}


#' Microscopy prevalence in an age bracket
#'
#' The proportion in the clinical and treated states, plus the proportion of
#' asymptomatic infections that microscopy would detect given each
#' individual's detection immunity, within an age bracket.
#'
#' @param ages Ages in days.
#' @param IDs Detection immunities.
#' @param states Infection states.
#' @param age_brackets Lower and upper age in days. Default is 2 to 10 years.
#' @return A single prevalence.
#' @keywords internal
microscopy_detected <- function(ages, IDs, states, age_brackets = c(2 * 365, 10 * 365)) {
  age_in <- ages < age_brackets[2] & ages > age_brackets[1]
  dt <- sum(states %in% c(1, 4) & age_in) / sum(age_in)
  mpl <- model_param_list_create()

  micro_det <- q_fun(
    mpl$d1, IDs, mpl$ID0, mpl$kD,
    sapply(ages, fd, mpl$fD0, mpl$aD, mpl$gammaD)
  )

  a <- sum(states %in% 2 & age_in) / sum(age_in) * mean(micro_det[age_in])

  return(dt + a)
}

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
#' @examples
#' \donttest{
#' out <- pipeline(
#'   EIR = 2, years = 1, N = 200, num_loci = 4, seed = 1L,
#'   update_save = TRUE, update_length = 365,
#'   output = save_loggers(final = "humans"),
#'   housekeeping_list = housekeeping_list_create(quiet = TRUE, cluster = TRUE)
#' )
#' humans <- out[[length(out)]]
#' coi <- convert_barcode_vectors(humans, humans$ID, nl = 4)
#' table(coi$COI)
#' }

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
    states_i <- sim_save$Strain_infection_state_vectors[[i]]
    n_unique <- function(state) length(unique(int.out[[i]][which(states_i == state)]))
    a <- a + round(n_unique(2) * (micro_det[i]^mpl$alphaA))
    a <- a + round(n_unique(3) * (micro_det[i]^mpl$alphaU))
    a <- a + length(unique(int.out[[i]][which(sim_save$Strain_infection_state_vectors[[i]] == 4)]))
    return(a)
  }

  # Calculate COI based on alternative PCR detectability method
  COI_pcr_alternative <- function(i) {
    st <- sim_save$Strain_infection_state_vectors[[i]]
    staying <- st == 1 | st == 4 | st == 2
    staying[st == 3] <- as.logical(rbinom(sum(st == 3), size = 1, micro_det[i]^mpl$alphaU))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }

  # Calculate COI based on if parasites are microscopically detected
  COI_patent <- function(i) {
    st <- sim_save$Strain_infection_state_vectors[[i]]
    staying <- st == 1 | st == 4
    staying[st == 2] <- as.logical(rbinom(sum(st == 2), size = 1, micro_det[i]))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }

  # calculate the COIs
  if (!ibd) {
    if (!sub_patents_included) {
      COI <- rep(0, length(int.out))
      for (i in seq_along(int.out)) {
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

#' Integer identities of a population's IBD barcodes
#'
#' Converts the last strain of each individual with at least one strain, by
#' \code{\link{convert_ibd_barcode}}, dropping individuals with none.
#'
#' @param barcode_vec List with one element per individual, each a list of
#'   logical barcodes.
#' @param nl Number of loci.
#' @return List of numeric identity vectors.
#' @keywords internal
population_ibd_barcodes <- function(barcode_vec, nl) {
  n.strains <- lapply(barcode_vec, length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains > 0)],
    function(x) {
      lapply(tail(x, 1), convert_ibd_barcode, nl) %>% unlist()
    }
  )
}
