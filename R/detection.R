# function to calculate probability of strain being detected by microscopy
q_fun <- function(d1, ID, ID0, kD, fd) {
  return(d1 + ((1 - d1) / (1 + ((ID / ID0)^kD) * fd)))
}

fd <- function(age, fD0, aD, gammaD) {
  return(1 - ((1 - fD0) / (1 + (age / aD)^gammaD)))
}


#' @noRd
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
