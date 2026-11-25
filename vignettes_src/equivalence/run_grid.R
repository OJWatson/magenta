# run from this folder (hipercow root): Rscript run_grid.R
library(hipercow)
N <- 1e5; YEARS <- 10; REPS <- 10
EIRS <- c(0.5, 1, 2, 4, 8, 16, 32, 64)
jobs <- expand.grid(rep = seq_len(REPS), eir = rev(EIRS),
                    model = c("magenta", "malariasimulation", "malariasimulation_p25"), stringsAsFactors = FALSE)
jobs$seed <- 50000L + 100L * match(jobs$eir, EIRS) + jobs$rep + c(magenta = 0L, malariasimulation = 5000L, malariasimulation_p25 = 10000L)[jobs$model]
if (file.exists("grid_ids.rds")) stop("grid already submitted; see grid_ids.rds")
r <- hipercow_resources(cores = 1L, max_runtime = "12h")
b <- task_create_bulk_expr(run_forensic(model, eir, rep, seed, N, YEARS, "results"),
                           jobs, resources = r, environment = "h2h")
saveRDS(list(bundle = b, jobs = jobs, submitted = Sys.time()), "grid_ids.rds")
print(b); print(length(b$ids))
