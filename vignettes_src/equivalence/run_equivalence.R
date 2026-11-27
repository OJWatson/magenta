## Runs behind the "Equivalence of magenta and malariasimulation" vignette:
## magenta and malariasimulation (with and without a 25-day prophylaxis
## period on the clinical treatment) at 8 EIRs from 0.5 to 64, 10 replicates
## each, 100,000 people for 10 years, and the deterministic Griffin
## equilibrium as the reference. The runs take 5 to 60 minutes each, so they are
## run once by hand on the DIDE cluster with hipercow, not when the package or
## the vignette is built. The vignette only reads the output,
## equivalence_summary.rds, which this script writes next to it.
##
## To rerun: copy this folder (funs.R, analysis.R and this script) to a folder
## on a network share the cluster can see, run this script from there, and
## copy equivalence_summary.rds back to vignettes_src/equivalence/. Then
## rebuild the vignette with
## `Rscript vignettes_src/build_vignettes.R Equivalence`.
##
## The runs install magenta and malariasimulation from their current versions
## (r-universe and GitHub), so the outputs come from the package versions at the
## time they were made (recorded in the meta element of equivalence_summary.rds).
##
## Set `steps` before sourcing this script to run only some of it: "run"
## (cluster setup, runs and deterministic reference) and "analyse" (summarise
## the saved results). Analysing needs only the files in results/, and
## deterministic.rds and versions.rds written by "run".

if (!exists("steps")) steps <- c("run", "analyse")
if (!exists("results_files")) results_files <- list.files("results", pattern = "\\.rds$", full.names = TRUE)
if (!exists("deterministic_file")) deterministic_file <- "deterministic.rds"
if (!exists("versions_file")) versions_file <- "versions.rds"
if (!exists("summary_file")) summary_file <- "equivalence_summary.rds"

if ("run" %in% steps) {
  library(hipercow)

  ## 1. Cluster environment ----
  writeLines(c("repo::https://mrc-ide.r-universe.dev",
               "OJWatson/magenta",
               "mrc-ide/malariasimulation"),
             "pkgdepends.txt")
  hipercow_init(driver = "dide-windows", r_version = "4.5.3")
  stopifnot(dide_check())
  hipercow_environment_create(name = "equivalence", sources = "funs.R", overwrite = TRUE)
  hipercow_provision(environment = "equivalence")

  ## 2. Runs ----
  ## 8 EIRs x 10 replicates x 3 models = 240 tasks. Seeds are paired across
  ## EIRs and replicates, and differ between models.
  N <- 1e5; YEARS <- 10; REPS <- 10
  EIRS <- c(0.5, 1, 2, 4, 8, 16, 32, 64)
  jobs <- expand.grid(rep = seq_len(REPS), eir = rev(EIRS),
                      model = c("magenta", "malariasimulation", "malariasimulation_p25"),
                      stringsAsFactors = FALSE)
  jobs$seed <- 50000L + 100L * match(jobs$eir, EIRS) + jobs$rep +
    c(magenta = 0L, malariasimulation = 5000L, malariasimulation_p25 = 10000L)[jobs$model]
  b <- task_create_bulk_expr(
    run_forensic(model, eir, rep, seed, N, YEARS, "results"),
    jobs, environment = "equivalence", bundle_name = "magenta_equivalence",
    resources = hipercow_resources(cores = 1L, max_runtime = "12h"))

  ## the deterministic reference and the package versions the runs used
  det_task <- task_create_expr(run_deterministic(), environment = "equivalence")
  ver_task <- task_create_expr(
    list(magenta = as.character(packageVersion("magenta")),
         malariasimulation = as.character(packageVersion("malariasimulation"))),
    environment = "equivalence")

  stopifnot(hipercow_bundle_wait(b, timeout = 24 * 3600, progress = FALSE))
  task_wait(det_task); task_wait(ver_task)
  saveRDS(task_result(det_task), deterministic_file)
  saveRDS(task_result(ver_task), versions_file)
}

## 3. Analyse the saved results ----
if ("analyse" %in% steps) {
  source("analysis.R")
  det <- readRDS(deterministic_file)
  res <- analyse_equivalence(results_files, det$global, det$bands,
                             versions = readRDS(versions_file), Y0 = 4, Y1 = 10)
  saveRDS(res$summary, summary_file, compress = "xz")

  h1 <- res$headline
  op <- options(width = 250)
  print(h1[h1$outcome %in% c("realised_EIR", "kappa", "PfPR2-10", "clinical_all", "clinical_u5", "severe_all", "severe_u5"),
           c("EIR", "outcome", "mean_magenta", "mean_ms_p25", "mean_ms", "det",
             "pct_magenta_vs_ms_p25", "signif_magenta_vs_ms_p25", "pct_magenta_vs_ms", "signif_magenta_vs_ms")],
        digits = 3, row.names = FALSE)
  cat("\nBurn-in drift, |z| > 3\n")
  print(res$drift[abs(res$drift$z) > 3, ], digits = 3, row.names = FALSE)
  options(op)
}
