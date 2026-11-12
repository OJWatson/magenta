test_that("equilibrium_init_create checks its inputs and honours seasonality", {
  mpl <- model_param_list_create()
  av <- magenta:::age_brackets(100, 10, TRUE)
  expect_error(equilibrium_init_create("a", 5, 0.4, 10, mpl), "age_vector provided is not numeric")
  expect_error(equilibrium_init_create(av, "5", 0.4, 10, mpl), "het_brackets")
  expect_error(equilibrium_init_create(av, 5, "0.4", 10, mpl), "ft provided")
  expect_error(equilibrium_init_create(av, 5, 0.4, "10", mpl), "EIR provided")
  expect_error(equilibrium_init_create(av, 5, 0.4, 10, mpl, country = 1), "country")
  expect_error(equilibrium_init_create(av, 5, 0.4, 10, mpl, admin = 1), "admin")

  flat <- equilibrium_init_create(av, 5, 0.4, 10, mpl, quiet = TRUE)
  expect_equal(flat$theta_c, 0)
  expect_equal(dim(flat$init_S), c(10, 5, 4))
  seasonal <- equilibrium_init_create(av, 5, 0.4, 10, mpl,
    country = "Democratic Republic of the Congo", admin = "Kinshasa", quiet = TRUE
  )
  expect_gt(seasonal$theta_c, 0)
  expect_equal(seasonal$mv0, flat$mv0)
  expect_message(
    equilibrium_init_create(av, 5, 0.4, 10, mpl,
      country = "Democratic Republic of the Congo", admin = "Kinshasa"
    ),
    "Returned: Kinshasa"
  )
})

test_that("spl_create builds island and closed spatial lists", {
  closed <- magenta:::spl_create(plaf = c(0.5, 0.5))
  expect_s3_class(closed, "magenta_spatial_list")
  expect_equal(closed$imported_cotransmissions_events, 0)
  island <- magenta:::spl_create(
    spatial_type = 1, human_importation_rate_vector = c(0.1, NA, 0.2),
    mosquito_imporation_rate_vector = c(1, 2), plaf = c(0.5, 0.5),
    island_imports_plaf_linked_flag = TRUE
  )
  expect_equal(island$imported_cotransmissions_events, 0.3)
  expect_equal(island$imported_oocyst_events, 3)
  expect_true(island$island_imports_plaf_linked_flag)
  expect_error(magenta:::spl_create(spatial_type = 2), "not implemented")
})
