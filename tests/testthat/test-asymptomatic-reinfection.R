test_that("reinfection in A keeps the person's and strains' state change days consistent", {

  out <- pipeline(EIR = 50, years = 0.5, N = 1000, ft = 0.4, full_save = TRUE, seed = 7L,
                  housekeeping_list = housekeeping_list_create(quiet = TRUE, quiet_test = TRUE,
                                                               cluster = TRUE))
  pl <- out$populations_event_and_strains_List
  states <- out$population_List$Infection_States
  last_day <- out$parameters_List$g_current_time - 1

  person_day <- pl$Day_of_InfectionStatus_change
  strain_states <- pl$Strain_infection_state_vectors
  strain_days <- pl$Strain_day_of_infection_state_change_vectors

  # everyone in D, A or U carries at least one strain
  expect_true(all(pl$Number_of_Strains[states %in% 1:3] > 0))

  # people in A: their next change is in the future, and no D or A strain outlasts it
  in_a <- which(states == 2)
  expect_gt(length(in_a), 0)
  expect_true(all(person_day[in_a] > last_day))
  patent_ok <- vapply(in_a, function(i) {
    patent <- strain_states[[i]] %in% c(1, 2)
    all(strain_days[[i]][patent] <= person_day[i] & strain_days[[i]][patent] > last_day)
  }, logical(1))
  expect_true(all(patent_ok))

  # people in U carry no patent strains
  in_u <- which(states == 3)
  expect_gt(length(in_u), 0)
  expect_false(any(vapply(in_u, function(i) any(strain_states[[i]] %in% c(1, 2)), logical(1))))

  # strain changes that are still to come are all in the future
  untreated <- unlist(strain_days)[unlist(strain_states) %in% 1:3]
  expect_true(all(untreated > last_day))

  # for people in D, A or U the next strain change is the earliest of the strains' days
  next_ok <- vapply(which(pl$Number_of_Strains > 0 & states %in% 1:3), function(i) {
    d <- strain_days[[i]][strain_states[[i]] != 4 & strain_days[[i]] != 0]
    length(d) == 0 || pl$Day_of_next_strain_state_change[i] == min(d)
  }, logical(1))
  expect_true(all(next_ok))

  # pending infections are all still to start, and are kept through treatment and prophylaxis
  pending <- lapply(seq_along(states), function(i) {
    times <- pl$Infection_time_realisation_vectors[[i]]
    times[seq_along(times) > pl$Number_of_Realised_Infections[i]]
  })
  expect_true(all(unlist(pending) > last_day))
  expect_true(any(lengths(pending)[states %in% 4:5] > 0))

})
