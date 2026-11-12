# The parameter list constructors must produce exactly the fields recorded in
# cpp_contract, in that order, and their validators must reject bad input with
# a message that names the list and the field.

contract <- magenta:::cpp_contract

test_that("each constructor produces exactly the contracted fields", {
  expect_identical(names(housekeeping_list_create()), contract$housekeeping$fields)
  expect_identical(names(drug_list_create()), contract$drug$fields)
  expect_identical(names(drug_create()), contract$drug_entry$fields)
  expect_identical(
    names(drug_create_default_no_resistance()), contract$drug_entry$fields
  )
  expect_identical(names(drug_create_al()), contract$drug_entry$fields)
  expect_identical(names(drug_create_asaq()), contract$drug_entry$fields)
  expect_identical(names(drug_create_dhappq()), contract$drug_entry$fields)
  expect_identical(
    names(vector_adaptation_list_create(1)), contract$vector_adaptation$fields
  )
  expect_identical(names(nmf_list_create()), contract$nmf$fields)
  expect_identical(names(barcode_list_create()), contract$barcode$fields)
  expect_identical(names(model_param_list_create()), contract$core$fields)
  expect_identical(names(magenta:::spl_create()), contract$spatial$fields)
})

test_that("the contracted C++ fields are a subset of the constructor fields", {
  for (nm in c(
    "housekeeping", "drug", "drug_entry", "vector_adaptation", "nmf",
    "barcode", "spatial", "eqss", "core"
  )) {
    expect_true(all(contract[[nm]]$cpp %in% contract[[nm]]$fields), info = nm)
  }
})

test_that("the C++ sources read exactly the contracted fields", {
  src <- system.file("..", "src", package = "magenta")
  skip_if(!dir.exists(src), "C++ sources not available")
  read_fields <- function(file, list_name) {
    lines <- readLines(file.path(src, file), warn = FALSE)
    hits <- regmatches(
      lines,
      gregexpr(sprintf('%s\\["([A-Za-z_0-9]+)"\\]', list_name), lines)
    )
    unique(sub('.*\\["([A-Za-z_0-9]+)"\\]', "\\1", unlist(hits)))
  }
  expect_setequal(read_fields("main_init.cpp", "param_list"), contract$init)
  expect_setequal(read_fields("main_init.cpp", "housekeeping_list"), contract$housekeeping$cpp)
  expect_setequal(read_fields("main_init.cpp", "barcode_list"), contract$barcode$cpp)
  expect_setequal(read_fields("main_init.cpp", "drug_list"), contract$drug$cpp)
  expect_setequal(read_fields("main_init.cpp", "drug"), contract$drug_entry$cpp)
  expect_setequal(
    read_fields("main_init.cpp", "vector_adaptation_list"), contract$vector_adaptation$cpp
  )
  expect_setequal(read_fields("main_init.cpp", "nmf_list"), contract$nmf$cpp)
  expect_setequal(read_fields("main_init.cpp", "core_parameter_list"), contract$core$cpp)
  expect_setequal(read_fields("main_init.cpp", "eqSS"), contract$eqss$cpp)
  expect_true(all(read_fields("main_init.cpp", "spatial_list") %in% contract$spatial$cpp))
  expect_true(all(read_fields("main_update.cpp", "spatial_list") %in% contract$spatial$cpp))
  expect_true(all(read_fields("main_update.cpp", "drug_list") %in% contract$drug$cpp))
  expect_true(all(read_fields("main_update.cpp", "barcode_list") %in% contract$barcode$cpp))
  expect_setequal(read_fields("main_saved_init.cpp", "savedState"), contract$saved_state)
})

test_that("constructors return classed lists that still behave as lists", {
  hk <- housekeeping_list_create()
  expect_s3_class(hk, "magenta_housekeeping_list")
  expect_s3_class(hk, "list")
  expect_true(is.list(hk))
  expect_identical(hk$quiet_print, TRUE)
  expect_output(print(hk), "<magenta_housekeeping_list> with 4 fields")

  dl <- drug_list_create()
  expect_s3_class(dl, "magenta_drug_list")
  expect_s3_class(dl$drugs[[1]], "magenta_drug")
})

test_that("validators name the list and field that is wrong", {
  hk <- housekeeping_list_create()
  hk$cluster <- "yes"
  expect_error(validate_housekeeping_list(hk), "housekeeping_list\\$cluster must be TRUE or FALSE")
  hk$cluster <- NULL
  expect_error(validate_housekeeping_list(hk), "missing the field\\(s\\): cluster")
  expect_error(validate_housekeeping_list(c(unclass(housekeeping_list_create()), x = 1)),
    "unexpected field\\(s\\): x"
  )
  expect_error(validate_housekeeping_list("a"), "housekeeping_list must be a list")

  dl <- drug_list_create()
  dl$drug_choice <- 5
  expect_error(validate_drug_list(dl), "drug_list\\$drug_choice must be a single whole number")
  dl <- drug_list_create()
  dl$partner_drug_ratios <- c(0.5, 0.5)
  expect_error(validate_drug_list(dl), "drug_list\\$partner_drug_ratios must have")
  dl <- drug_list_create()
  dl$drugs[[1]]$lpf <- 2
  expect_error(validate_drug_list(dl), "drug_list\\$drugs\\[\\[1\\]\\]\\$lpf must be")
  expect_error(drug_list_create(temporal_cycling = 1, sequential_cycling = 1),
    "Both sequential and temporal"
  )

  expect_error(barcode_list_create(num_loci = 4, mutation_rate = rep(0, 3)),
    "barcode_list\\$mutation_rate must have length 4"
  )
  expect_error(barcode_list_create(num_loci = 4, plaf = rep(2, 4)),
    "barcode_list\\$plaf must be a number vector between 0 and 1"
  )

  expect_error(nmf_list_create(prob_of_testing_nmf = 1.5), "nmf_list\\$prob_of_testing_nmf")
  expect_error(nmf_list_create(mean_nmf_frequency = 1:3, nmf_age_brackets = 1:3),
    "nmf_list\\$nmf_age_brackets must have length 4"
  )

  expect_error(vector_adaptation_list_create(1, local_oocyst_advantage = -1),
    "vector_adaptation_list\\$local_oocyst_advantage"
  )

  expect_error(model_param_list_create(itn_cov = 2), "model_param_list\\$itn_cov")
  mpl <- model_param_list_create()
  mpl$eta <- "a"
  expect_error(validate_model_param_list(mpl), "model_param_list\\$eta must be a single number")

  expect_error(magenta:::spl_create(spatial_type = 1, plaf = c(0.5, 2)), "spatial_list\\$plaf")
})

test_that("model_param_list_create derives the same quantities as before", {
  mpl <- model_param_list_create()
  expect_equal(mpl$fv0, 1 / (0.69 + 2.31))
  expect_equal(mpl$av0, 0.92 * mpl$fv0)
  expect_equal(mpl$Surv0, exp(-0.132 * 10))
  expect_equal(mpl$cov, c(1, 0, 0, 0))
  expect_equal(mpl$pop_split, rep(0.25, 4))
  expect_equal(mpl$irs_loss, log(2) / (0.5 * 365))
  expect_equal(mpl$lambda, -0.5 * mpl$b_lambda + sqrt(0.25 * mpl$b_lambda^2 +
    13.25 * 21.2 * 0.0348 * 6.64 / (2 * 0.0338 * 0.132 * 3.72 * (1 + 0.643 * 0.249))))
  mpl2 <- model_param_list_create(itn_cov = 0.5, irs_cov = 0.2)
  expect_equal(mpl2$cov, c(0.5 * 0.8, 0.5 * 0.8, 0.5 * 0.2, 0.5 * 0.2))
})

test_that("the equilibrium state validator catches the wrong shapes", {
  mpl <- model_param_list_create()
  eq <- suppressMessages(equilibrium_init_create(
    age_vector = magenta:::age_brackets(100, 10, TRUE), het_brackets = 3,
    ft = 0.4, EIR = 10, model_param_list = mpl, quiet = TRUE
  ))
  ss <- magenta:::equilibrium_ss_create(eq)
  expect_s3_class(ss, "magenta_eqss")
  expect_identical(names(ss), contract$eqss$fields)
  ss$Smat <- ss$Smat[1:2, ]
  expect_error(magenta:::validate_eqss(ss), "eqSS\\$Smat must be a 10 by 3 numeric matrix")
})

test_that("the init and update parameter lists are validated and ordered", {
  mpl <- model_param_list_create()
  eq <- suppressMessages(equilibrium_init_create(
    age_vector = magenta:::age_brackets(100, 10, TRUE), het_brackets = 3,
    ft = 0.4, EIR = 10, model_param_list = mpl, quiet = TRUE
  ))
  ss <- magenta:::equilibrium_ss_create(eq)
  pl <- magenta:::param_list_simulation_init_create(
    N = 10, eqSS = ss, barcode_list = barcode_list_create(),
    spatial_list = magenta:::spl_create(), housekeeping_list = housekeeping_list_create(),
    drug_list = drug_list_create(), nmf_list = nmf_list_create(),
    vector_adaptation_list = vector_adaptation_list_create(1), mpl = mpl
  )
  expect_identical(names(pl), contract$init)
  expect_error(
    magenta:::param_list_simulation_init_create(
      N = 10, eqSS = unclass(ss)[-1], barcode_list = barcode_list_create(),
      spatial_list = magenta:::spl_create(), housekeeping_list = housekeeping_list_create(),
      drug_list = drug_list_create(), nmf_list = nmf_list_create(),
      vector_adaptation_list = vector_adaptation_list_create(1), mpl = mpl
    ),
    "eqSS is missing the field\\(s\\): age_brackets"
  )
  expect_error(
    magenta:::param_list_simulation_update_create(
      years = 1, statePtr = "not a pointer", spatial_list = magenta:::spl_create(),
      drug_list = drug_list_create(), barcode_list = barcode_list_create()
    ),
    "statePtr must be an external pointer"
  )
})
