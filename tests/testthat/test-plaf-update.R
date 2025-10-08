test_that("update_saved_state_barcode_plaf adjusts barcode PLAF", {
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
  
  target <- c(0.5, 1)
  
  updated <- update_saved_state_barcode_plaf(saved, target, seed = 42)
  
  get_plaf <- function(barcode_list) {
    flattened <- magenta:::flatten_barcode_vectors(barcode_list)
    mat <- magenta:::flat_to_matrix(flattened$flat, length(target))
    colMeans(mat)
  }
  
  strain_plaf <- get_plaf(updated$populations_event_and_strains_List$Strain_barcode_vectors)
  infection_plaf <- get_plaf(updated$populations_event_and_strains_List$Infection_barcode_realisation_vectors)
  male_plaf <- get_plaf(updated$scourge_List$Mosquito_Oocyst_barcode_male_vectors)
  female_plaf <- get_plaf(updated$scourge_List$Mosquito_Oocyst_barcode_female_vectors)
  
  expect_equal(strain_plaf, target)
  expect_equal(infection_plaf, target)
  expect_equal(male_plaf, target)
  expect_equal(female_plaf, target)
  
  combined_matrix <- rbind(
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$populations_event_and_strains_List$Strain_barcode_vectors)$flat, length(target)),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$populations_event_and_strains_List$Infection_barcode_realisation_vectors)$flat, length(target)),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$scourge_List$Mosquito_Oocyst_barcode_male_vectors)$flat, length(target)),
    magenta:::flat_to_matrix(magenta:::flatten_barcode_vectors(updated$scourge_List$Mosquito_Oocyst_barcode_female_vectors)$flat, length(target))
  )
  
  expect_equal(colMeans(combined_matrix), target)
})

test_that("update_saved_state_barcode_plaf validates inputs", {
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
    update_saved_state_barcode_plaf(saved, c(0.5)),
    "length"
  )
  
  saved$populations_event_and_strains_List$Strain_barcode_vectors <- list()
  expect_error(
    update_saved_state_barcode_plaf(saved, c(0, 1)),
    "does not contain any barcode vectors"
  )
})