test_that("bitset roundtrip", {
  for (i in 1:10) {
    x <- as.raw(sample(c(1, 0), 24, replace = TRUE))
    y <- expect_output(
      test_bitset_serialisation(x, length(x)),
      paste(as.integer(x), collapse = "")
    )
    expect_identical(y, x)
  }

  for (i in 1:10) {
    x <- as.raw(sample(c(1, 0), 48, replace = TRUE))
    y <- expect_output(
      test_bitset_serialisation(x, length(x)),
      paste(as.integer(x), collapse = "")
    )
    expect_identical(y, x)
  }
})
