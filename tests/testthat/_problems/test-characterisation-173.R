# Extracted from test-characterisation.R:173

# test -------------------------------------------------------------------------
skip_characterisation_if_unsupported()
mpl <- model_param_list_create()
age_vector <- magenta:::age_brackets(100, 40, TRUE)
eq <- suppressMessages(equilibrium_init_create(
  age_vector = age_vector, het_brackets = 5, ft = 0.4, EIR = 10,
  model_param_list = mpl, quiet = TRUE
))
ss <- magenta:::equilibrium_ss_create(eq)
expect_snapshot_value(
  list(
    mv0 = eq$mv0, FOIv_eq = eq$FOIv_eq, omega = eq$omega,
    state_totals = vapply(ss[grep("mat$", names(ss))], sum, numeric(1)),
    Sv = ss$Sv, Ev = ss$Ev, Iv = ss$Iv
  ),
  style = "json2", tolerance = 1e-8
)
