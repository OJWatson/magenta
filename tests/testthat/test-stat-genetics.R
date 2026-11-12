test_that("zero truncated draws have the right length and no zeros", {
  set.seed(1)
  g <- magenta:::ztrgeomintp(500, 10, 0.2)
  expect_length(g, 500)
  expect_true(all(g >= 1))
  nb <- magenta:::ztrnbinom(500, mean = 2.5, size = 1)
  expect_length(nb, 500)
  expect_true(all(nb >= 1))
  expect_true(mean(nb) > 1)
})

test_that("the summarySE helpers group and summarise", {
  df <- data.frame(
    Age_Bin = rep(c("a", "b"), each = 4),
    State = rep(c("x", "y"), 4),
    COI = c(1, 2, 3, 4, 5, 6, 7, NA)
  )
  full <- magenta:::summarySE(df, "COI", c("Age_Bin", "State"))
  expect_named(full, c("Age_Bin", "State", "N", "mean", "sd", "se", "ci"))
  expect_equal(full$mean[full$Age_Bin == "a" & full$State == "x"], 2)
  m <- magenta:::summarySE_mean_only(df, "COI", c("Age_Bin", "State"))
  expect_named(m, c("Age_Bin", "State", "N", "mean"))
  expect_equal(m$N[m$Age_Bin == "b" & m$State == "y"], 1)
  expect_equal(
    magenta:::summarySE_mean_only(df, "COI", "Age_Bin", mean_only = FALSE),
    magenta:::summarySE(df, "COI", "Age_Bin")
  )
  capped <- magenta:::summarySE_mean_only_max_mean(df, "COI", "Age_Bin", max = 3)
  expect_equal(capped$mean[capped$Age_Bin == "b"], 3)
  expect_equal(nrow(magenta:::summarySE_mean_only_max_mean(df[0, ], "COI", "Age_Bin", max = 3)), 0)
})

test_that("clonality and cou are computed from barcode lists", {
  bl <- list(c(1, 1, 2), c(2), c(3), c(4))
  # unique barcodes per individual: 1,2 | 2 | 3 | 4 -> counts 1:1, 2:2, 3:1, 4:1
  expect_equal(unname(magenta:::clonality_from_barcode_list(bl)), 3 / 5)
  expect_equal(magenta:::clonality_from_barcode_list(list(c(1), c(1))), 0)
  cou <- magenta:::cou_from_barcode_list(bl)
  expect_true(cou >= 0 && cou <= 1)
  expect_equal(magenta:::cou_from_barcode_list(list(1, 2, 3)), 0)
})

test_that("IBD barcodes convert to integer vectors", {
  # two loci, two bits each: locus 1 = bits 1,2 ; locus 2 = bits 3,4
  b <- c(TRUE, FALSE, FALSE, TRUE)
  expect_equal(magenta:::convert_ibd_barcode(b, nl = 2), c(1, 2))
  expect_equal(magenta:::bitsToInt(c(TRUE, TRUE, FALSE)), 3)
  expect_equal(magenta:::bitsToInt(c(TRUE, TRUE, FALSE), big_endian = TRUE), 6)
})

test_that("barcode flattening and rebuilding round trip", {
  bl <- list(list(c(TRUE, FALSE), c(FALSE, FALSE)), list(), list(c(TRUE, TRUE)))
  flat <- magenta:::flatten_barcode_vectors(bl)
  mat <- magenta:::flat_to_matrix(flat$flat, 2)
  expect_equal(dim(mat), c(3, 2))
  expect_equal(magenta:::rebuild_barcode_vectors(mat, flat$counts), bl)
  expect_error(magenta:::validate_plaf_vector(c(0.5, 1.5)), "between 0 and 1")
  expect_error(magenta:::validate_plaf_vector("a"), "numeric")
  nested <- list(a = list(b = 1))
  expect_true(magenta:::has_nested_element(nested, c("a", "b")))
  expect_false(magenta:::has_nested_element(nested, c("a", "c")))
  expect_equal(magenta:::get_nested_element(nested, c("a", "b")), 1)
  expect_equal(magenta:::set_nested_element(nested, c("a", "b"), 2)$a$b, 2)
})
