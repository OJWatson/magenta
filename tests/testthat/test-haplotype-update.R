test_that("update_saved_state_barcode_haplotype adjusts haplotype frequencies", {
  saved <- list(
    populations_event_and_strains_List = list(
      Strain_barcode_vectors = list(
        list(c(FALSE, TRUE)),
        list(c(FALSE, TRUE)),
        list(c(TRUE, TRUE), c(FALSE, FALSE)),
        list(c(FALSE, TRUE)),
        list(c(TRUE, FALSE)),
        list()
      ),
      Infection_barcode_realisation_vectors = list(
        list(c(FALSE, FALSE)),
        list(c(TRUE, TRUE)),
        list(c(FALSE, TRUE)),
        list(c(FALSE, TRUE))
      )
    ),
    scourge_List = list(
      Mosquito_Oocyst_barcode_male_vectors = list(
        list(c(TRUE, FALSE)),
        list(c(FALSE, FALSE)),
        list(c(TRUE, TRUE)),
        list(c(FALSE, TRUE))
      ),
      Mosquito_Oocyst_barcode_female_vectors = list(
        list(c(FALSE, TRUE)),
        list(c(FALSE, TRUE)),
        list(c(TRUE, TRUE)),
        list(c(TRUE, FALSE))
      )
    )
  )

  target <- c(0.1829, 0.1954, 0.1403, 0.4814)

  updated <- update_saved_state_barcode_haplotype(saved, target, seed = 123)

  get_haplotype_freq <- function(barcode_list) {
    flattened <- magenta:::flatten_barcode_vectors(barcode_list)
    mat <- magenta:::flat_to_matrix(flattened$flat, 2)
    if (nrow(mat) == 0) {
      return(rep(0, length(target)))
    }
    idx <- apply(mat, 1, function(row) sum(row * 2^(seq_along(row) - 1))) + 1
    tabulate(idx, nbins = length(target)) / nrow(mat)
  }

  strain_freq <- get_haplotype_freq(updated$populations_event_and_strains_List$Strain_barcode_vectors)
  infection_freq <- get_haplotype_freq(updated$populations_event_and_strains_List$Infection_barcode_realisation_vectors)
  male_freq <- get_haplotype_freq(updated$scourge_List$Mosquito_Oocyst_barcode_male_vectors)
  female_freq <- get_haplotype_freq(updated$scourge_List$Mosquito_Oocyst_barcode_female_vectors)

  expect_equal(strain_freq, target)
  expect_equal(infection_freq, target)
  expect_equal(male_freq, target)
  expect_equal(female_freq, target)

  combined_matrix <- rbind(
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$populations_event_and_strains_List$Strain_barcode_vectors)$flat, 2),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$populations_event_and_strains_List$Infection_barcode_realisation_vectors)$flat, 2),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$scourge_List$Mosquito_Oocyst_barcode_male_vectors)$flat, 2),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$scourge_List$Mosquito_Oocyst_barcode_female_vectors)$flat, 2)
  )

  combined_idx <- apply(combined_matrix, 1, function(row) sum(row * 2^(seq_along(row) - 1))) + 1
  combined_freq <- tabulate(combined_idx, nbins = length(target)) / nrow(combined_matrix)

  expect_equal(combined_freq, target)
})

test_that("update_saved_state_barcode_haplotype validates inputs", {
  saved <- list(
    populations_event_and_strains_List = list(
      Strain_barcode_vectors = list(list(c(TRUE, FALSE))),
      Infection_barcode_realisation_vectors = list(list(c(TRUE, FALSE)))
    ),
    scourge_List = list(
      Mosquito_Oocyst_barcode_male_vectors = list(list(c(TRUE, FALSE))),
      Mosquito_Oocyst_barcode_female_vectors = list(list(c(TRUE, FALSE)))
    )
  )

  expect_error(
    update_saved_state_barcode_haplotype(saved, c(0.5, 0.5, 0.5)),
    "length"
  )

  expect_error(
    update_saved_state_barcode_haplotype(saved, c(0.5, 0.5, 0.1, -0.1)),
    "non-negative"
  )

  expect_error(
    update_saved_state_barcode_haplotype(saved, c(0.25, 0.25, 0.25, 0.2)),
    "sum to 1"
  )

  saved$populations_event_and_strains_List$Strain_barcode_vectors <- list()
  expect_error(
    update_saved_state_barcode_haplotype(saved, c(0.25, 0.25, 0.25, 0.25)),
    "does not contain any barcode vectors"
  )
})
