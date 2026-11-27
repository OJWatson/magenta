# Whether a bite infects is decided when the infection would start, g_dur_E (12) days after the
# bite, as in the Griffin et al. model. A bite on day 2 gives an infection starting on day 14.
# States: 0 S, 1 D, 2 A, 3 U, 4 T, 5 P
follow <- function(state, day_of_change, treated = FALSE) {
  test_pending_infection(state, day_of_change, bite_day = 2L, treated = treated, days = 16L)
}

test_that("a bite in D, T or P infects if the person is S, A or U when it starts", {

  for (state in c(1L, 4L, 5L)) {
    out <- follow(state, day_of_change = 5L)
    # pending from the bite until it starts, and no ICA boost until then
    expect_equal(out$pending[2:13], rep(1, 12))
    expect_equal(out$strains[13], 0)
    expect_equal(out$ICA[13], 0)
    # it starts on day 14
    expect_equal(out$pending[14], 0)
    expect_equal(out$strains[14], 1)
    expect_equal(out$state[14], 2)
    expect_gt(out$ICA[14], 0.99)
  }

  # a treated outcome is a clinical case on the day it starts, in someone treated at the bite
  out <- follow(4L, day_of_change = 5L, treated = TRUE)
  expect_equal(out$state[14], 4)
  expect_equal(out$incidence, replace(rep(0, 16), 14, 1))

})

test_that("pending infections are kept through recovery and treatment", {

  # U -> S on day 5, and T -> P on day 5 then P -> S the day after
  for (state in c(3L, 4L)) {
    out <- follow(state, day_of_change = 5L)
    expect_equal(out$state[6:13], rep(0, 8))
    expect_equal(out$pending[6:13], rep(1, 8))
    expect_equal(out$state[14], 2)
    expect_equal(out$strains[14], 1)
  }

})

test_that("an infection does not start in someone in D, T or P", {

  for (state in c(1L, 4L, 5L)) {
    out <- follow(state, day_of_change = 30L, treated = TRUE)
    expect_equal(out$state, rep(state, 16))
    expect_equal(out$pending[14:16], c(0, 0, 0))
    expect_equal(out$strains, rep(0, 16))
    expect_equal(out$ICA, rep(0, 16))
    expect_equal(out$incidence, rep(0, 16))
  }

  # whereas someone susceptible is infected
  out <- follow(0L, day_of_change = 0L)
  expect_equal(out$state[14], 2)
  expect_equal(out$strains[14], 1)

})
