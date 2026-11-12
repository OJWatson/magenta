# metric functions
## --------------------------------------

#' @noRd
unique_function <- function(r, year0, year, ss, unphased, metric, samp, pos, tab_func) {
  if (!unphased) {
    t <- c(unlist(r[[year0]][["ints"]][samp[pos <= ss]]), unlist(r[[year]][["ints"]][samp[pos > ss]]))
    return(sum(!(duplicated(t) | duplicated(t, fromLast = TRUE))) / length(t))
  } else {
    t <- c(unlist(r[[year0]][["ints_max"]][samp[pos <= ss]]), unlist(r[[year]][["ints_max"]][samp[pos > ss]]))
    return(sum(!(duplicated(t) | duplicated(t, fromLast = TRUE))) / ss)
  }
}

#' @noRd
cou_function <- function(r, year0, year, ss, unphased, metric, samp, pos, tab_func) {
  l_factor_i <- attr(r[[year0]], "l") + attr(r[[year]], "l")
  if (!unphased) {
    t <- c(unlist(r[[year0]][["ints"]][samp[pos <= ss]]), unlist(r[[year]][["ints"]][samp[pos > ss]]))
    z <- sum(tab_func(t, l_factor_i)^2)
    return((z - (1 / length(t))) / (1 - (1 / length(t))))
  } else {
    t <- c(unlist(r[[year0]][["ints_max"]][samp[pos <= ss]]), unlist(r[[year]][["ints_max"]][samp[pos > ss]]))
    z <- sum(tab_func(t, l_factor_i)^2)
    return((z - (1 / length(t))) / (1 - (1 / length(t))))
  }
}

#' @noRd
pibd_function <- function(r, year0, year, ss, unphased, metric, samp, pos, tab_func) {
  l_factor_i <- attr(r[[year0]], "l") + attr(r[[year]], "l")
  if (!unphased) {
    t <- rbind(rbind_list_base(r[[year0]][["ints"]][samp[pos <= ss]]), rbind_list_base(r[[year]][["ints"]][samp[pos > ss]]))
    z <- mean(apply(t, 2, function(x) sum(tab_func(x, l_factor_i)^2)))
    return((z - (1 / nrow(t))) / (1 - (1 / nrow(t))))
  } else {
    t <- rbind(rbind_list_base(r[[year0]][["ints_max"]][samp[pos <= ss]]), rbind_list_base(r[[year]][["ints_max"]][samp[pos > ss]]))
    z <- mean(apply(t, 2, function(x) sum(tab_func(x, l_factor_i)^2)))
    return((z - (1 / nrow(t))) / (1 - (1 / nrow(t))))
  }
}

#' @noRd
mean_function <- function(r, year0, year, ss, unphased, metric, samp, pos, tab_func) {
  mean(c(r[[year0]][[metric]][samp[pos <= ss]], r[[year]][[metric]][samp[pos > ss]]), na.rm = TRUE)
}


bootstrap <- function(sim = "M:/OJ/magenta_Results/scripts/fig3_r1.rds", year0 = 21, year = 22,
                      ss = 10, sub_sample = 100, less = TRUE, conf = .95, permutations = 100,
                      metric = "coi", unphased = FALSE, ages = c(5, 15), states = c(1, 2, 4)) {
  # possible metric list
  metrics <- c(
    "% Polygenomic" = "polygenom", "% Unique" = "unique",
    "COI" = "coi", "COU" = "cou",
    "iIBD" = "pibd_within", "pIBD" = "pibd"
  )
  if (!(metric %in% metrics)) stop("Wrong metric")

  # tabulation function
  tab_func <- function(thing, l_factor_i) {
    (tabulate(thing, l_factor_i)) / length(thing)
  }
  # tab_func <- table

  # read in the sim
  r <- readRDS(sim)
  res <- data.frame(
    "ss" = ss,
    "power" = 0,
    "prev" = r[[year]]$pcr_prev,
    "Year" = year - year0
  )

  # number of positions available to us for sampling
  npos0 <- which(r[[year0]]$state %in% states & r[[year0]]$age < ages[2] * 365 & r[[year0]]$age > ages[1] * 365)
  npos1 <- which(r[[year]]$state %in% states & r[[year]]$age < ages[2] * 365 & r[[year]]$age > ages[1] * 365)

  # replace sampling decision
  if (any(c(length(npos0), length(npos1)) < ss * 1.5)) {
    message("Low sample size possibility warning")
    replace <- TRUE
  } else {
    replace <- FALSE
  }

  # special functions
  metric_func <- match.fun("mean_function")
  metric_func_opts <- c("cou_function", "unique_function", "pibd_function")
  metric_func_choice <- grep(metric, metric_func_opts)
  if (length(metric_func_choice) == 1) {
    metric_func <- match.fun(metric_func_opts[metric_func_choice])
  }

  ## --------------------------------------
  power <- 0
  for (i in seq_len(sub_sample)) {
    samp0 <- sample(npos0, size = ss, replace = replace)
    samp1 <- sample(npos1, size = ss, replace = replace)

    null <- metric_func(r = r, year0 = year0, year = year, ss = ss, unphased = unphased, metric = metric, samp = c(samp0), pos = 1:ss, tab_func = tab_func)
    alt <- metric_func(r = r, year0 = year0, year = year, ss = ss, unphased = unphased, metric = metric, samp = c(samp1), pos = (ss + 1):(2 * ss), tab_func = tab_func)
    diff <- mean(alt - null)

    all <- c(samp0, samp1)
    permute <- t(replicate(permutations, sample(ss * 2, size = ss * 2, FALSE)))
    permuteds <- apply(permute, 1, function(x) {
      null <- metric_func(r = r, year0 = year0, year = year, ss = ss, unphased = unphased, metric = metric, samp = all[x[1:ss]], pos = x[1:ss], tab_func = tab_func)
      alt <- metric_func(r = r, year0 = year0, year = year, ss = ss, unphased = unphased, metric = metric, samp = all[x[(ss + 1):(2 * ss)]], pos = x[(ss + 1):(2 * ss)], tab_func = tab_func)
      return(mean(alt - null))
    })
    if (less) {
      if (sum(permuteds > diff, na.rm = TRUE) > conf * permutations) {
        power <- power + 1
      }
    } else {
      if (sum(permuteds < diff, na.rm = TRUE) > conf * permutations) {
        power <- power + 1
      }
    }
  }

  res$power <- power / sub_sample
  return(res)
}
