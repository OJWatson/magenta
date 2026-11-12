# kullback-leibler calculation
kl <- function(res_path, loc, df_real, densities, sample_sizes, sub_patents_included = TRUE, COI_type = "pcr_imperial") {
  # calcualation of KL entropy as taken from entropy package
  kl_ent <- function(t) {
    kl_p <- function(freqs1, freqs2, unit = c("log", "log2", "log10")) {
      unit <- match.arg(unit)
      freqs1 <- freqs1 / sum(freqs1)
      freqs2 <- freqs2 / sum(freqs2)
      if (any(!(freqs2 > 0))) {
        warning("Vanishing value(s) in argument freqs2!")
      }
      LR <- ifelse(freqs1 > 0, log(freqs1 / freqs2), 0)
      KL <- sum(freqs1 * LR)
      if (unit == "log2") {
        KL <- KL / log(2)
      }
      if (unit == "log10") {
        KL <- KL / log(10)
      }
      return(KL)
    }

    mi_p <- function(freqs2d, unit = c("log", "log2", "log10")) {
      unit <- match.arg(unit)
      freqs2d <- as.matrix(freqs2d / sum(freqs2d))
      freqs.x <- rowSums(freqs2d)
      freqs.y <- colSums(freqs2d)
      freqs.null <- freqs.x %o% freqs.y
      MI <- kl_p(freqs2d, freqs.null, unit = unit)
      return(MI)
    }

    return(mi_p(t / sum(t), unit = c("log", "log2", "log10")))
  }

  # load our dataset and get its length
  x <- readRDS(res_path)
  pos <- length(x)

  # sample from the simulation the COIs
  s <- sample_coi(x[[pos]],
    ID = x[[pos - 1]]$ID,
    sample_size = sample_sizes[loc],
    age_densities = densities[, loc],
    reps = 50,
    sub_patents_included = sub_patents_included,
    COI_type = COI_type
  )

  # generate the comparison and real dataframe
  df_comp <- data.frame(
    "ages" = df_real$ages[df_real$loc == names(sample_sizes[loc])],
    "coi" = df_real$coi[df_real$loc == names(sample_sizes[loc])]
  )

  df <- data.frame(
    "ages" = cut(s$sim_save$Ages[s$ids[, 1]] / 365, seq(0, 90, 2), include.lowest = TRUE),
    "coi" = cut(s$COI[s$ids[, 1]], seq(1, 25, 1), include.lowest = TRUE)
  )


  # now loop through our data frame
  kl_reps <- rep(0, nrow(table(df)))
  kls <- rep(0, 50)

  for (j in 1:50) {
    # create iteration comparison df
    df <- data.frame(
      "ages" = cut(s$sim_save$Ages[s$ids[, j]] / 365, seq(0, 90, 2), include.lowest = TRUE),
      "coi" = cut(s$COI[s$ids[, j]], seq(1, 25, 1), include.lowest = TRUE)
    )

    # loop through each each age bracket
    kl <- rep(0, nrow(table(df)))
    tb <- table(df)
    tb_real <- table(df_comp)

    for (i in 1:nrow(table(df))) {
      if (sum(sum(tb[i, ])) != 0) {
        kl_reps[i] <- kl_ent(t(rbind(tb[i, ], tb_real[i, ])))
      }
    }

    # normalise across the age brackets
    kls[j] <- sum(kl_reps * table(df_comp$ages)) / length(df_comp$ages)
  }

  return(kls)
}
