## Submit the comparison: 2 units x 3 seeds x 2 models, N 50k, 20 y burn-in,
## 2000-2024 (each task about an hour). Check first with a tiny task per model.
library(hipercow)
source("funs.R")

test <- task_create_bulk_expr(
  if (model == "magenta") run_magenta(unit, 1L, 1000L) else run_malsim(unit, 1L, 1000L),
  data.frame(unit = "Busia urban", model = c("magenta", "malariasimulation")),
  environment = "site", resources = hipercow_resources(queue = "Testing"))
stopifnot(hipercow_bundle_wait(test))
print(do.call(rbind, hipercow_bundle_result(test))[c(1:2, 26:27), ])

grid <- expand.grid(unit = vapply(units, `[[`, "", "name"), seed = 1:3,
                    model = c("magenta", "malariasimulation"), stringsAsFactors = FALSE)
grid$N <- 50000L
b <- task_create_bulk_expr(
  if (model == "magenta") run_magenta(unit, seed, N) else run_malsim(unit, seed, N),
  grid, environment = "site", bundle_name = "magenta_site_vignette",
  resources = hipercow_resources(cores = 1, max_runtime = "1d", memory_per_node = "8G"))
print(table(hipercow_bundle_status(b)))
