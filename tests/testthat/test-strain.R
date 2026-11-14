test_that("random from plaf", {
  # with 1 and 0 should be no variation
  plaf <- c(rep(1, 18), rep(0, 18))
  y <- test_barcode_from_PLAF(plaf, length(plaf))
  x <- as.raw(plaf)
  expect_identical(y, x)

  # from a differing plaf
  plaf <- dnorm(1:47, mean = 12, sd = 5)
  y <- replicate(100000, as.numeric(test_barcode_from_PLAF(plaf, length(plaf))))
  x <- rowMeans(y)
  expect_true(sum(abs(plaf - x) < mean(plaf)) == length(plaf))
})

test_that("depedent plaf", {
  # with 1 should be no variation
  plaf <- c(1, 0.5, 0, 0, 0.5)
  y <- test_dependent_barcode_from_PLAF(plaf, length(plaf))
  x <- as.raw(c(1, 1, 0, 0, 1))
  expect_identical(y, x)
})

test_that("ibd recombinants", {
  # with 1 and 0 should be no variation
  x <- as.raw(c(1, 1, 1, 0, 0, 1, 0, 0))
  y <- as.raw(c(0, 0, 0, 1, 1, 0, 1, 1))
  z <- test_recombinant_with_ibd(x, y, length(x), 4, 2, rep(1, 4))
  expected <- c(y[1:2], x[3:8])
  expect_identical(expected, z)
})


test_that("ibd next generation", {
  z <- test_generate_next_ibd(bl = 48, nl = 24, ib = 2, pc = rep(1, 24), id = 2)
  y <- rep(intToBits(2)[1:2], 24)
  expect_identical(y, z)
})


test_that("ibd conversion", {
  barcode <- as.raw(c(1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0))
  z <- test_ibd_conversion(barcode = barcode, nl = 3, ib = 4, bl = 12)
  expect_identical(z$vec, c(7, 2, 5))
})

test_that("the fitness cost is read at the resistance loci, the last one included", {
  # three resistance loci spread along an eight-locus barcode; the cost table
  # from drug_list_create() has bit j for the j-th locus
  loci <- c(2, 5, 7)
  dl <- drug_list_create(
    number_of_resistance_loci = 3, cost_of_resistance = c(0.9, 0.8, 0.7),
    drugs = list(drug_create(barcode_res_pos = loci))
  )
  expect_identical(dl$resistance_loci, loci)
  cost <- function(mutant_at) {
    barcode <- rep(0, 8)
    barcode[mutant_at + 1] <- 1
    test_relative_contribution(as.raw(barcode), 8L, loci, dl$cost_of_resistance)
  }
  expect_equal(cost(integer(0)), 1)
  expect_equal(cost(2), 0.9)
  expect_equal(cost(5), 0.8)
  # a mutant only at the last resistance locus pays that locus's cost
  expect_equal(cost(7), 0.7)
  expect_equal(cost(c(2, 7)), 0.9 * 0.7)
  expect_equal(cost(c(2, 5, 7)), 0.9 * 0.8 * 0.7)
  # positions that are not resistance loci cost nothing
  expect_equal(cost(c(0, 1, 3, 4, 6)), 1)
})
