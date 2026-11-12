test_that("age_brackets gives geometric or linear bounds", {
  geo <- magenta:::age_brackets(100, 10, TRUE)
  expect_length(geo, 10)
  expect_equal(geo[1], 0)
  expect_equal(geo[10], 100)
  expect_true(all(diff(geo) > 0))
  lin <- magenta:::age_brackets(100, 20, FALSE)
  expect_equal(lin, seq(0, 100, 20))
})

test_that("mu_fv_create returns constant rates with no vector interventions", {
  out <- magenta:::mu_fv_create(eqInit = NULL, ft = 0.4, itn_cov = 0, irs_cov = 0, years = 2)
  expect_equal(nrow(out), 730)
  expect_equal(unique(out$mu), 0.132)
  expect_equal(unique(out$fv), 0.333)
  # a padded vector of zeros is still "no interventions"
  out2 <- magenta:::mu_fv_create(NULL, ft = c(0.4, 0.4), itn_cov = c(0, 0), irs_cov = 0, years = 2)
  expect_equal(out, out2)
})

test_that("mu_fv_create runs the deterministic model when coverage is non-zero", {
  skip_on_cran()
  mpl <- model_param_list_create()
  eq <- equilibrium_init_create(
    age_vector = magenta:::age_brackets(100, 40, TRUE), het_brackets = 5,
    ft = 0.4, EIR = 5, model_param_list = mpl, quiet = TRUE
  )
  out <- suppressMessages(magenta:::mu_fv_create(eq,
    ft = 0.4, itn_cov = 0.5, irs_cov = 0, years = 1
  ))
  expect_equal(nrow(out), 365)
  expect_named(out, c("mu", "fv"))
  expect_true(all(is.finite(out$mu)) && all(is.finite(out$fv)))
  # nets raise mosquito mortality and lower the feeding rate on humans
  expect_gt(out$mu[365], 0.132)
  expect_lt(out$fv[365], 1 / 3)
  full <- suppressMessages(magenta:::mu_fv_create(eq,
    ft = 0.4, itn_cov = 0.5, irs_cov = 0.2, years = 1, full = TRUE
  ))
  expect_true(all(c("mu", "fv", "t") %in% names(full)))
})

test_that("drug_list_update switches drugs under sequential and temporal cycling", {
  dl <- drug_list_create(
    resistance_flag = TRUE, number_of_drugs = 2,
    drugs = list(drug_create_al(), drug_create_asaq()),
    sequential_cycling = 0.1, sequential_update = 2, partner_drug_ratios = c(1, 0)
  )
  expect_equal(dl$drug_choice, 0)
  expect_equal(dl$next_temporal_cycle, -1)
  # below the failure threshold: nothing happens
  expect_equal(magenta:::drug_list_update(dl, year = 1, tf = 0.05)$drug_choice, 0)
  # above it: switch, and schedule the next possible switch
  expect_message(
    dl2 <- magenta:::drug_list_update(dl, year = 1, tf = 0.2),
    "switched to drug 1"
  )
  expect_equal(dl2$drug_choice, 1)
  expect_equal(dl2$next_temporal_cycle, 3)
  # too soon for another
  expect_equal(magenta:::drug_list_update(dl2, year = 2, tf = 0.5)$drug_choice, 1)
  # and wraps round to the first drug
  dl3 <- suppressMessages(magenta:::drug_list_update(dl2, year = 4, tf = 0.5))
  expect_equal(dl3$drug_choice, 0)

  tc <- drug_list_create(
    resistance_flag = TRUE, number_of_drugs = 2,
    drugs = list(drug_create_al(), drug_create_asaq()),
    temporal_cycling = 3, partner_drug_ratios = c(1, 0)
  )
  expect_equal(magenta:::drug_list_update(tc, year = 2, tf = NA)$drug_choice, 0)
  tc2 <- magenta:::drug_list_update(tc, year = 3, tf = NA)
  expect_equal(tc2$drug_choice, 1)
  expect_equal(tc2$next_temporal_cycle, 6)
  expect_equal(magenta:::drug_list_update(tc2, year = 6, tf = NA)$drug_choice, 0)

  # no resistance or mft: never switches
  expect_equal(magenta:::drug_list_update(drug_list_create(), 1, 1)$drug_choice, 0)
  mft <- drug_list_create(
    resistance_flag = TRUE, mft_flag = TRUE, number_of_drugs = 2,
    drugs = list(drug_create_al(), drug_create_asaq()), temporal_cycling = 1
  )
  expect_equal(magenta:::drug_list_update(mft, 1, 1)$drug_choice, 0)
})

test_that("pop_alf and lineages summarise raw barcode matrices", {
  # two individuals: barcodes (1,0) and (1,1), then (0,0), over two loci
  nums <- list(
    matrix(as.raw(c(1, 0, 1, 1)), ncol = 2, byrow = TRUE),
    matrix(as.raw(c(0, 0)), ncol = 2)
  )
  af <- magenta:::pop_alf(nums, nl = 2)
  expect_equal(af, c(2 / 3, 1 / 3))
  afw <- magenta:::pop_alf(nums, nl = 2, weighted = TRUE)
  expect_equal(afw, c(0.5, 0.25))
  lin <- magenta:::lineages(nums, nl = 2)
  expect_length(lin, 4)
  expect_equal(as.integer(lin), c(1, 1, 0, 1))
  expect_true(all(is.na(magenta:::pop_alf(list(), nl = 2))))
  expect_true(all(is.na(magenta:::lineages(list(), nl = 2))))
})

test_that("spl_grab pads year ranges outside 2000-2017", {
  early <- suppressMessages(spl_grab("Democratic Republic of the Congo", "Kinshasa", 1998:2001))
  expect_equal(nrow(early$incidence), 4)
  expect_true(all(early$incidence[1:2, ] == 0))
  late <- suppressMessages(spl_grab("Democratic Republic of the Congo", "Kinshasa", 2016:2019))
  expect_equal(nrow(late$incidence), 4)
  expect_equal(late$mosquitoFOI[3, ], late$mosquitoFOI[4, ])
  expect_equal(ncol(late$incidence), ncol(late$mosquitoFOI))
})

test_that("intervention_grab extends, truncates and overrides the historic series", {
  cc <- "Democratic Republic of the Congo"
  late <- suppressMessages(intervention_grab(cc, "Kinshasa", 2010:2020,
    final_itn_cov = 0.9, final_irs_cov = 0.1, final_ft = 0.7
  ))
  expect_length(late$itn_cov, 11)
  expect_equal(late$itn_cov[11], 0.9)
  expect_equal(late$irs_cov[11], 0.1)
  expect_equal(late$ft[11], 0.7)
  # linear interpolation between 2015 and the end year
  expect_true(all(diff(late$ft[6:11]) >= 0) || all(diff(late$ft[6:11]) <= 0))

  early <- suppressMessages(intervention_grab(cc, "Kinshasa", 1998:2002))
  expect_length(early$itn_cov, 5)
  expect_equal(early$itn_cov[1:2], c(0, 0))

  none <- suppressMessages(intervention_grab(cc, "Kinshasa", 1990:1995))
  expect_equal(none$itn_cov, rep(0, 6))
})

test_that("the schedule helpers read a period and recent treatment failure", {
  s <- build_schedule(years = 1, update_length = 100, ft = 0.3, plaf = c(0.1, 0.2))
  p <- magenta:::schedule_period(s, 2)
  expect_equal(p$i, 2)
  expect_equal(p$n, 4)
  expect_equal(p$day_start, 101)
  expect_equal(p$ft, 0.3)
  expect_equal(p$plaf, c(0.1, 0.2))
  expect_null(p$importation_human)

  res <- list(list(overall_treatment_failure = 0.2), list(overall_treatment_failure = 0.4))
  expect_equal(magenta:::recent_treatment_failure(res, 1), 0.2)
  expect_equal(magenta:::recent_treatment_failure(res, 2), 0.3)
  expect_true(is.na(magenta:::recent_treatment_failure(list(list(a = 1)), 1)))

  expect_equal(magenta:::spatial_type_code(NULL), 0)
  expect_equal(magenta:::spatial_type_code("island"), 1)
  expect_error(magenta:::spatial_type_code("metapop"), "not implemented")
})

test_that("intervention_grab falls back to zero IRS where there is no IRS series", {
  togo <- suppressMessages(intervention_grab("Togo", "Maritime", 2010:2012))
  expect_equal(togo$irs_cov, rep(0, 3))
  expect_length(togo$itn_cov, 3)
  expect_true(all(togo$itn_cov > 0))
})

test_that("spatial types map to the codes the C++ reads", {
  expect_equal(magenta:::spatial_type_code(NULL), 0)
  expect_equal(magenta:::spatial_type_code("island"), 1)
  expect_equal(magenta:::spatial_type_code(1), 1)
  expect_error(magenta:::spatial_type_code("metapop"), "not implemented")
})
