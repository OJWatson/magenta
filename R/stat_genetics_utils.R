
#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars Characters of grouping variables
#' @param conf.interval numeric for CI
#' 
#' @keywords internal
summarySE <- function(data=NULL, measurevar="COI", 
                      groupvars=c("Age_Bin","State"),
                      conf.interval=.95) {
  
  
  res <- dplyr::group_by_at(data,dplyr::vars(dplyr::one_of(groupvars))) %>% 
    dplyr::summarise(N=length(!!dplyr::sym(measurevar)),
                     mean=mean(!!dplyr::sym(measurevar),na.rm=TRUE),
                     sd = sd(!!dplyr::sym(measurevar),na.rm=TRUE),
                     se = sd/sqrt(.data$N),
                     ci = .data$se * suppressWarnings(qt(conf.interval/2 + .5, .data$N-1))
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
summarySE_mean_only <- function(data=NULL, measurevar="COI", groupvars=c("Age_Bin","State"),mean_only=TRUE) {
  
  if(!mean_only) {
    return(summarySE(data, measurevar,groupvars))
  } else {
    var <-  dplyr::sym(measurevar)
    
    res <- dplyr::group_by_at(data,groupvars) %>% 
      dplyr::summarise(N=sum(!is.na(!!var)),
                       mean=mean(!!var,na.rm=TRUE))
    
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
summarySE_mean_only_max_mean <- function(data=NULL, 
                                         measurevar="COI", 
                                         groupvars=c("Age_Bin","State"),
                                         mean_only=TRUE, 
                                         max = 0) {
  
  if (nrow(data) != 0) {
    if (max) {
      data[measurevar][data[measurevar]>max] <- max
    }
  }
  if(!mean_only) {
    return(summarySE(data, measurevar,groupvars))
  } else {
    var <-  dplyr::sym(measurevar)
    
    res <- dplyr::group_by_at(data,groupvars) %>% 
      dplyr::summarise(N=length(!!var),
                       mean=mean(!!var,na.rm=TRUE))
    
    return(as.data.frame(res))
  }
}

# zero truncated geometric integers
ztrgeomintp <- function(n, mean, p){
  
  ngs <- rgeom(n*1.5,prob = 1/mean)
  ngs <- round(ngs*p)
  ngs <- ngs[ngs>0]
  while(length(ngs) < n) {
    ng2 <- rgeom(n*1.5,prob = 1/mean)
    ngs <- c(ngs,round(ng2*p))
    ngs <- ngs[ngs>0]
  }
  
  ngs <- sample(ngs,size = n,replace=FALSE)
  return(ngs)
  
  
}

# zero truncated negative binomial
ztrnbinom <- function(n,mean,size) {
  
  
  nbs <- rnbinom(n*1.2,size = size, mu = mean)
  nbs <- nbs[nbs>0]
  while(length(nbs) < n) {
    nbs <- c(nbs,rnbinom(n*.5,size = size, mu = mean))
    nbs <- nbs[nbs>0]
  }
  
  nbs <- sample(nbs,size = n,replace=FALSE)
  return(nbs)
}

# get the clonality from a list of barcode numbers
clonality_from_barcode_list <- function(barcode_list){
  
  tbl <- table(table(unlist(lapply(barcode_list,unique))))
  
  if ("1" %in% names(tbl)){
    return(tbl[1]/sum(tbl*as.numeric(names(tbl))))
  } else {
    return(0)
  }
  
}

# get the cou from a list of barcode numbers
cou_from_barcode_list <- function(barcode_list){
  
  tbl <- table(unlist(lapply(barcode_list,unique)))
  samp_size <- sum(tbl)
  z <- sum((tbl/samp_size)^2)
  return((z - (1/samp_size))/(1-(1/samp_size)))
  
}


# convert ibd barcode of logicals to a vector of integers
convert_ibd_barcode <- function(b, nl){
  
  ib <- length(b)/nl
  
  br <- seq_len(ib)
  ibd <- rep(0,nl)
  for(i in seq_len(nl)){
    ibd[i] <- bitsToInt(b[br+(ib*(i-1))])
  }  
  
  return(ibd)
  
}

# convert population's ibd barcodes to list of int vectors
population_ibd_barcodes <- function(barcode_vec,nl){
  
  n.strains <- lapply(barcode_vec,length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains>0)],
    function(x){
      lapply(tail(x,1),convert_ibd_barcode,nl) %>% unlist
    }
  )
  
}

# covert populations' ibd barcodes to list of int vectors using c++
population_ibd_barcodes_c <- function(barcode_vec,bl,nl,ib){
  
  n.strains <- lapply(barcode_vec,length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains>0)],
    function(x){
      lapply(as.raw(tail(x[[1]],1)),test_ibd_conversion,bl,nl,ib) %>% unlist
    }
  )
  
}


# get the pIBD from a list of barcode numbers
pibd_from_barcode_list <- function(barcode_list, l_factor_i){
  
  t <- rbind_list_base(barcode_list)
  if(nrow(t)>1){
    z <- mean(apply(t,2,function(x) sum((tabulate(x,l_factor_i)/length(x))^2)))
    return((z - (1/nrow(t)))/(1-(1/nrow(t))))
  } else {
    return(NA)
  }
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

update_saved_state_barcode_haplotype <- function(saved_state, target_haplotype, seed = NULL) {

  if (!is.list(saved_state)) {
    stop("`saved_state` must be a list produced by readRDS on a saved simulation state.")
  }

  validate_haplotype_vector(target_haplotype)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  list_paths <- list(
    c("populations_event_and_strains_List", "Strain_barcode_vectors"),
    c("populations_event_and_strains_List", "Infection_barcode_realisation_vectors"),
    c("scourge_List", "Mosquito_Oocyst_barcode_male_vectors"),
    c("scourge_List", "Mosquito_Oocyst_barcode_female_vectors")
  )

  for (path in list_paths) {
    if (!has_nested_element(saved_state, path)) {
      stop(sprintf("Saved state does not contain `%s`.", paste(path, collapse = "$")))
    }

    barcode_list <- get_nested_element(saved_state, path)
    flattened <- flatten_barcode_vectors(barcode_list)

    if (length(flattened$flat) == 0) {
      stop(sprintf(
        "Cannot adjust haplotype frequencies for `%s` because it does not contain any barcode vectors.",
        paste(path, collapse = "$")
      ))
    }

    barcode_length <- length(flattened$flat[[1]])
    haplotype_count <- length(target_haplotype)
    if (haplotype_count != 2^barcode_length) {
      stop(sprintf(
        "`target_haplotype` has length %d but barcodes in `%s` have length %d (requiring %d haplotypes).",
        haplotype_count, paste(path, collapse = "$"), barcode_length, 2^barcode_length
      ))
    }

    barcode_matrix <- flat_to_matrix(flattened$flat, barcode_length)
    adjusted_matrix <- adjust_matrix_to_haplotype(barcode_matrix, target_haplotype)

    saved_state <- set_nested_element(
      saved_state,
      path,
      rebuild_barcode_vectors(adjusted_matrix, flattened$counts)
    )
  }

  invisible(saved_state)
}

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

validate_haplotype_vector <- function(target_haplotype) {
  if (!is.numeric(target_haplotype)) {
    stop("`target_haplotype` must be a numeric vector.")
  }
  if (length(target_haplotype) == 0) {
    stop("`target_haplotype` must have positive length.")
  }
  if (any(target_haplotype < 0)) {
    stop("`target_haplotype` values must be non-negative.")
  }
  total <- sum(target_haplotype)
  if (!isTRUE(all.equal(total, 1, tolerance = .Machine$double.eps^0.5))) {
    stop("`target_haplotype` values must sum to 1.")
  }
  invisible(target_haplotype)
}

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

get_nested_element <- function(x, path) {
  current <- x
  for (name in path) {
    current <- current[[name]]
  }
  current
}

set_nested_element <- function(x, path, value) {
  if (length(path) == 1) {
    x[[path[[1]]]] <- value
    return(x)
  }
  
  head_path <- path[[1]]
  x[[head_path]] <- set_nested_element(x[[head_path]], path[-1], value)
  x
}

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

adjust_matrix_to_haplotype <- function(barcode_matrix, target_haplotype) {
  n <- nrow(barcode_matrix)
  if (n == 0) {
    return(barcode_matrix)
  }

  barcode_matrix <- as.matrix(barcode_matrix)
  storage.mode(barcode_matrix) <- "integer"

  barcode_length <- ncol(barcode_matrix)
  haplotypes <- length(target_haplotype)
  if (haplotypes != 2^barcode_length) {
    stop("Number of haplotypes does not match barcode length.")
  }

  target_counts <- compute_target_haplotype_counts(target_haplotype, n)

  assignments <- rep.int(seq_len(haplotypes), target_counts)
  if (length(assignments) != n) {
    stop("Failed to construct haplotype assignments of the correct size.")
  }

  if (n > 1) {
    assignments <- sample(assignments)
  }

  for (i in seq_len(n)) {
    barcode_matrix[i, ] <- haplotype_index_to_barcode(assignments[i] - 1L, barcode_length)
  }

  validate_haplotype(barcode_matrix, target_haplotype, "barcode vectors")

  barcode_matrix
}

compute_target_haplotype_counts <- function(target_haplotype, n) {
  raw <- target_haplotype * n
  target_counts <- floor(raw)
  remainder <- n - sum(target_counts)

  if (remainder < 0) {
    stop("`target_haplotype` results in impossible haplotype counts.")
  }

  if (remainder > 0) {
    order_remainder <- order(raw - target_counts, decreasing = TRUE)
    for (idx in seq_len(remainder)) {
      target_counts[order_remainder[[idx]]] <- target_counts[order_remainder[[idx]]] + 1L
    }
  }

  storage.mode(target_counts) <- "integer"
  target_counts
}

haplotype_index_to_barcode <- function(index, barcode_length) {
  bits <- intToBits(as.integer(index))
  as.integer(bits[seq_len(barcode_length)])
}

validate_haplotype <- function(barcode_matrix, target_haplotype, context) {
  if (nrow(barcode_matrix) == 0) {
    return(invisible(NULL))
  }

  barcode_length <- ncol(barcode_matrix)
  haplotypes <- 2^barcode_length
  observed_counts <- integer(haplotypes)

  powers <- 2^(seq_len(barcode_length) - 1L)
  for (i in seq_len(nrow(barcode_matrix))) {
    idx <- as.integer(sum(barcode_matrix[i, ] * powers) + 1L)
    observed_counts[[idx]] <- observed_counts[[idx]] + 1L
  }

  expected_counts <- compute_target_haplotype_counts(target_haplotype, nrow(barcode_matrix))

  if (!identical(observed_counts, expected_counts)) {
    stop(sprintf(
      "Failed to match target haplotype frequencies for %s.",
      context
    ))
  }

  invisible(NULL)
}

validate_plaf <- function(barcode_matrix, target_plaf, context) {
  if (nrow(barcode_matrix) == 0) {
    return(invisible(NULL))
  }
  
  if (!identical(round(target_plaf*nrow(barcode_matrix)), colSums(barcode_matrix))) {
    stop(sprintf(
      "Failed to match target PLAF for %s. Expected %s but obtained %s.",
      context,
      paste(paste(round(target_plaf*nrow(barcode_matrix)), collapse = ", "), "out of", nrow(barcode_matrix)),
      paste(colSums(barcode_matrix), collapse = ", ")
    ))
  }
  
  invisible(NULL)
}
