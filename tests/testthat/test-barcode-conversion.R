make_sim_save <- function(strains, states) {
  list(
    Ages = rep(2000, length(strains)),
    Strain_barcode_vectors = strains,
    Strain_infection_state_vectors = states
  )
}

test_that("convert_barcode_vectors counts unique barcodes per individual", {
  strains <- list(
    list(as.raw(c(1, 0, 0, 0)), as.raw(c(1, 0, 0, 0)), as.raw(c(0, 1, 0, 0))),
    list(as.raw(c(1, 1, 0, 0))),
    list()
  )
  states <- list(c(1, 1, 1), 2, integer(0))
  sim_save <- make_sim_save(strains, states)

  res <- convert_barcode_vectors(sim_save, ID = c(0, 0, 0), nl = 4)
  expect_named(res, c("nums", "COI"))
  expect_equal(res$nums[[1]], c(1, 1, 2))
  expect_equal(res$nums[[2]], 3)
  expect_equal(res$COI, c(2, 1, 0))

  # with zero detection immunity every strain is detected, so excluding
  # sub-patents gives the same answer for all three COI definitions
  for (type in c("pcr_imperial", "pcr_alternative", "patent")) {
    res_type <- convert_barcode_vectors(
      sim_save,
      ID = c(0, 0, 0), nl = 4,
      sub_patents_included = FALSE, COI_type = type
    )
    expect_equal(res_type$COI, c(2, 1, 0), info = type)
  }

  expect_error(
    convert_barcode_vectors(sim_save, ID = c(0, 0, 0), COI_type = "nonsense"),
    "COI_type is not correct"
  )
})

test_that("convert_barcode_vectors converts IBD barcodes locus by locus", {
  # two loci of two bits each: locus 1 = bits 1:2, locus 2 = bits 3:4
  strains <- list(
    list(as.raw(c(0, 0, 0, 0)), as.raw(c(1, 0, 0, 1))),
    list(as.raw(c(1, 1, 1, 0))),
    list()
  )
  states <- list(c(1, 1), 1, integer(0))
  sim_save <- make_sim_save(strains, states)

  res <- convert_barcode_vectors(sim_save, ID = c(0, 0, 0), ibd = TRUE, nl = 2)

  # only the last strain of each infected individual is converted
  expect_length(res$nums, 2)
  expect_equal(res$nums[[1]], c(1, 2))
  expect_equal(res$nums[[2]], c(3, 1))
  expect_equal(res$COI, c(1, 1, 0))
})
