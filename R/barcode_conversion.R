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

# convert a population's IBD barcodes (last strain per individual) to integer
# vectors, dropping individuals with no strains
population_ibd_barcodes <- function(barcode_vec, nl) {
  n.strains <- lapply(barcode_vec, length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains > 0)],
    function(x) {
      lapply(tail(x, 1), convert_ibd_barcode, nl) %>% unlist()
    }
  )
}
