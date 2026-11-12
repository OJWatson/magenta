# function for pooling resistance results from didehpc cluster
#' @noRd
pool_res <- function(grp, nl = 3, N = 100000, reps = 10,
                     methods = c(
                       rep("Sequential Cycling", reps),
                       rep("Multiple First Line Therapies", reps),
                       rep("Sequential Cycling \n with 100% RDT testing", reps)
                     ),
                     meta = NULL,
                     non_cluster = FALSE,
                     list_returns = 1:5,
                     remove_zeros = TRUE) {
  pop_alf <- function(nums, nl = 3) {
    if (class(nums %>% unlist()) == "raw") {
      colMeans(matrix(as.numeric(do.call(rbind, nums)), ncol = nl))
    } else {
      rep(NA, nl)
    }
  }

  make_allele_df <- function(r, i, method, allele_only = TRUE, nl = 3) {
    if (allele_only) {
      af <- lapply(r[1:(length(r) - 1)], function(x) x$af) %>% unlist()
    } else {
      af <- lapply(r[1:(length(r) - 1)], function(x) pop_alf(x$df$nums[unlist(lapply(x$df$barcode_states, length)) > 0], nl)) %>% unlist()
    }
    times <- sort(rep(seq(30, by = 30, length.out = (length(r) - 1)), nl))
    loci <- rep(c(paste0("L", 1:(nl))), (length(r) - 1))
    df <- data.frame("AF" = af, "Time" = times / 365, "Drug_Locus" = loci)
    df$rep <- i
    df$Method <- method

    if (remove_zeros) {
      df <- na.omit(df)
      df <- df[!df$AF == 0, ]
    }
    return(df)
  }

  make_treat_fail_df <- function(r, i, method) {
    treatment_failures <- lapply(r[1:(length(r) - 1)], function(x) x$overall_treatment_failure) %>% unlist()
    unsuccesful_treatments_lpf <- lapply(r[1:(length(r) - 1)], function(x) x$unsuccesful_treatments_lpf) %>% unlist()
    times <- sort(rep(seq(30, by = 30, length.out = (length(r) - 1)), 1))
    df <- data.frame(
      "Treatment_Failure" = treatment_failures,
      "Num_Treatment_Failures" = unsuccesful_treatments_lpf,
      "Time" = times / 365
    )
    df$rep <- i
    df$Method <- method

    if (remove_zeros) {
      df <- na.omit(df)
      df <- df[!df$Treatment_Failure == 0, ]
    }
    return(df)
  }

  make_prevs_df <- function(r, i, method) {
    pcr <- c(lapply(r[1:(length(r) - 1)], function(x) x[["pcr_prev"]]) %>% unlist())
    micro <- c(lapply(r[1:(length(r) - 1)], function(x) x[["micro_2_10"]]) %>% unlist())
    dta <- c(lapply(r[1:(length(r) - 1)], function(x) x[["dat"]]) %>% unlist())
    times <- sort(rep(seq(30, by = 30, length.out = (length(r) - 1)), 1))
    df <- data.frame("pcr" = pcr, "micro" = micro, "dta" = dta, "Time" = times / 365)
    df$rep <- i
    df$Method <- method
    return(df)
  }

  make_lineage_df <- function(r, i, method) {
    lins <- c(lapply(r[1:(length(r) - 1)], function(x) x$lineage / sum(x$lineage)) %>% unlist())
    times <- sort(rep(seq(30, by = 30, length.out = (length(r) - 1)), 2^nl))
    df <- data.frame("Strains" = lins, "Time" = times / 365, "Lineage" = rep(0:((2^nl) - 1)))
    df$rep <- i
    df$Method <- method

    if (remove_zeros) {
      df <- na.omit(df)
      df <- df[!df$Strains == 0, ]
    }
    return(df)
  }

  make_mutations_df <- function(r, i, method, nl = 3) {
    muts <- lapply(r[1:(length(r) - 1)], function(x) x$mutations) %>% unlist()
    times <- sort(rep(seq(30, by = 30, length.out = (length(r) - 1)), nl))
    loci <- rep(c(paste0("L", 1:(nl))), (length(r) - 1))
    df <- data.frame("Mutations" = muts, "Time" = times / 365, "Drug_Locus" = loci)
    df$rep <- i
    df$Method <- method
    df$Year <- as.numeric(as.factor(cut(df$Time, 0:100)))

    if (remove_zeros) {
      df <- na.omit(df)
      df <- df[!df$Mutations == 0, ]
    }
    return(df)
  }

  af_list <- list()
  treat_fail_list <- list()
  prevs_list <- list()
  lineage_list <- list()
  mutation_list <- list()

  for (i in 1:length(methods)) {
    # message(paste0(i," "),appendLF = FALSE)
    found <- FALSE
    if (non_cluster) {
      r <- grp[[i]]
      found <- TRUE
    } else {
      if (grp$tasks[[i]]$status() == "COMPLETE") {
        x <- grp$tasks[[i]]
        path <- x$root$db$driver$name_hash(x$root$db$driver$get_hash(x$id, "task_results"))
        if (file.exists(path)) {
          r <- readRDS(path)
          found <- TRUE
        }
      }
    }

    if (found) {
      if (1 %in% list_returns) {
        af_list[[i]] <- make_allele_df(r, i, methods[i], nl = nl)
      }
      if (2 %in% list_returns) {
        mutation_list[[i]] <- make_mutations_df(r, i, methods[i], nl = nl)
      }
      if (3 %in% list_returns) {
        treat_fail_list[[i]] <- make_treat_fail_df(r, i, methods[i])
      }
      if (4 %in% list_returns) {
        prevs_list[[i]] <- make_prevs_df(r, i, methods[i])
      }
      if (5 %in% list_returns) {
        lineage_list[[i]] <- make_lineage_df(r, i, methods[i])
      }
    }
  }

  if (1 %in% list_returns) {
    af_df <- rbind_list_base(af_list)
    rm(af_list)
  } else {
    af_df <- NULL
  }
  if (2 %in% list_returns) {
    mutation_df <- rbind_list_base(mutation_list)
    rm(mutation_list)
  } else {
    mutation_df <- NULL
  }
  if (3 %in% list_returns) {
    fail_df <- rbind_list_base(treat_fail_list)
    rm(treat_fail_list)
  } else {
    fail_df <- NULL
  }
  if (4 %in% list_returns) {
    prev_df <- rbind_list_base(prevs_list)
    rm(prevs_list)
  } else {
    prev_df <- NULL
  }
  if (5 %in% list_returns) {
    lineage_df <- rbind_list_base(lineage_list)
    rm(lineage_list)
  } else {
    lineage_df <- NULL
  }

  if (!is.null(meta)) {
    for (m in seq_len(length(meta))) {
      if (1 %in% list_returns) {
        if (nrow(af_df) > 0) {
          af_df[names(meta[m])] <- meta[[m]]
        }
      }
      if (2 %in% list_returns) {
        if (nrow(mutation_df) > 0) {
          mutation_df[names(meta[m])] <- meta[[m]]
        }
      }
      if (3 %in% list_returns) {
        if (nrow(fail_df) > 0) {
          fail_df[names(meta[m])] <- meta[[m]]
        }
      }
      if (4 %in% list_returns) {
        if (nrow(prev_df) > 0) {
          prev_df[names(meta[m])] <- meta[[m]]
        }
      }
      if (5 %in% list_returns) {
        if (nrow(lineage_df) > 0) {
          lineage_df[names(meta[m])] <- meta[[m]]
        }
      }
    }
  }

  l <- list("af_df" = af_df, "mutation_df" = mutation_df, "fail_df" = fail_df, "prev_df" = prev_df, "lineage_df" = lineage_df)
  return(l)
}

# function for looping over didehpc cluster results
#' @noRd
create_res_list <- function(gs, nms, nums, nl = 2, reps = 50, methods = rep(reps, NA), list_returns = 1:5, remove_zeros = FALSE) {
  grab_number <- function(x, y) {
    na.omit(sapply(strsplit(x, "_|N_|starting")[[1]], as.numeric))[y]
  }

  pb <- progress::progress_bar$new(total = length(gs))
  res_list <- list()
  for (i in seq_len(length(gs))) {
    pb$tick()
    nm <- names(gs[i])
    if (all(gs[[i]]$status() == "COMPLETE")) {
      meta <- list()
      for (m in seq_along(nms)) {
        meta[[nms[m]]] <- grab_number(nm, nums[m])
      }
      res_list[[nm]] <- pool_res(gs[[i]], nl = nl, reps = reps, methods = methods, meta = meta, remove_zeros = remove_zeros)
      res_list[[nm]] <- res_list[[nm]][list_returns]
    }
  }

  return(res_list)
}
