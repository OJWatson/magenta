## Runs behind the "Running from a site file" vignette: magenta and
## malariasimulation for two Kenya site file units, 3 seeds each, 50,000
## people, a 20 year burn-in and then 2000-2024. The runs take up to an hour
## each, so they are run once by hand on the DIDE cluster with hipercow, not
## when the package or the vignette is built. The vignette only reads the
## output, site_summary.rds, which this script writes next to it.
##
## To rerun: copy this folder (funs.R, kenya_units.rds and this script) to a
## folder on a network share the cluster can see, run this script from there,
## and copy site_summary.rds back to vignettes_src/site/. Then rebuild the
## vignette with `Rscript vignettes_src/build_vignettes.R SiteFiles`.
##
## The runs install magenta from its GitHub default branch, so the outputs come
## from the magenta version in DESCRIPTION at the time they were made (recorded
## in site_summary.rds).

library(hipercow)
source("funs.R")

## 1. Site file ----
## kenya_units.rds holds the two units used, cut from the Kenya site file
## (version malariaverse_06_2026) in the foresite calibration packet
## 20260612-145526-c809460a (calibrated_scaled_site.rds, 12.7 MB), with the
## boundary shapes dropped. To remake it, set `full_site` to that file.
full_site <- NULL
if (!is.null(full_site)) {
  s <- readRDS(full_site)
  cut <- lapply(units, function(u) {
    v <- strsplit(u, " ")[[1]]
    k <- s$sites$name_1 == v[1] & s$sites$urban_rural == v[2]
    r <- site::subset_site(s, s$sites[k, ])
    r$shape <- NULL
    r
  })
  saveRDS(stats::setNames(cut, units), unit_sites, compress = "xz")
}

## 2. Cluster environment ----
writeLines(c("repo::https://mrc-ide.r-universe.dev",
             "OJWatson/magenta",
             "mrc-ide/malariasimulation",
             "mrc-ide/site",
             "mrc-ide/netz"),
           "pkgdepends.txt")
hipercow_init(driver = "dide-windows", r_version = "4.5.3")
hipercow_environment_create(name = "site", sources = "funs.R", overwrite = TRUE)
hipercow_provision(environment = "site")

## 3. Runs ----
## a tiny task per model first
test <- task_create_bulk_expr(
  run_model(model, "Busia urban", 1L, 1000L),
  data.frame(model = models),
  environment = "site", resources = hipercow_resources(queue = "Testing"))
stopifnot(hipercow_bundle_wait(test))

grid <- expand.grid(unit = units, seed = 1:3, model = models, stringsAsFactors = FALSE)
grid$N <- 50000L
b <- task_create_bulk_expr(
  run_model(model, unit, seed, N),
  grid, environment = "site", bundle_name = "magenta_site_comparison",
  resources = hipercow_resources(cores = 1, max_runtime = "1d", memory_per_node = "8G"))
stopifnot(hipercow_bundle_wait(b, timeout = 8 * 3600, progress = FALSE))

## 4. Collect the yearly summaries, with the site's own inputs and PfPR2-10 ----
runs <- do.call(rbind, hipercow_bundle_result(b))
site <- do.call(rbind, lapply(units, function(u) {
  s <- read_unit(u)
  pick <- function(d, col) d[[col]][match(years, d$year)]
  data.frame(unit = u, year = years,
             site_pfpr = pick(s$prevalence, "pfpr"),
             tx_cov = pick(s$interventions$treatment$implementation, "tx_cov"),
             itn_use = pick(s$interventions$itn$use, "itn_use"),
             irs_cov = pick(s$interventions$irs$implementation, "irs_cov"),
             eir = s$eir$eir[s$eir$sp == "pf"])
}))
saveRDS(list(runs = runs, site = site), "site_summary.rds")
