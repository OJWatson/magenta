# Non-malarial fevers only give drugs to people who seek treatment (ft), so with ft = 0 no one
# should ever be treated or in prophylaxis. States: 0 S, 1 D, 2 A, 3 U, 4 T, 5 P
nmf_states <- function(ft, prob_of_testing_nmf) {
  hk <- housekeeping_list_create(quiet = TRUE, quiet_test = TRUE, cluster = TRUE)
  nmf <- nmf_list_create(nmf_flag = TRUE, mean_nmf_frequency = rep(10, 14),
                         prob_of_testing_nmf = prob_of_testing_nmf)
  out <- pipeline(EIR = 0.1, N = 1000, years = 0.5, num_loci = 2, ft = ft, seed = 11L,
                  full_save = TRUE, nmf_list = nmf, housekeeping_list = hk)
  out$population_List$Infection_States
}

test_that("non-malarial fevers give no prophylaxis without treatment seeking", {
  skip_on_cran()

  expect_false(any(nmf_states(ft = 0, prob_of_testing_nmf = 0.5) %in% c(4, 5)))

  # untested fevers in people who seek treatment are given drugs
  expect_gt(sum(nmf_states(ft = 1, prob_of_testing_nmf = 0) == 5), 100)

  # and tested fevers in uninfected people are negative, so no drugs
  expect_lt(sum(nmf_states(ft = 1, prob_of_testing_nmf = 1) == 5), 50)
})
