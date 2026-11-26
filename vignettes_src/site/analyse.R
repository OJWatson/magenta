## Collect the runs, with the site file's own inputs and PfPR2-10, into
## site_summary.rds for the SiteFiles vignette (copy it to vignettes_src/).
library(hipercow)
source("funs.R")
res <- do.call(rbind, hipercow_bundle_result("magenta_site_vignette"))
site <- do.call(rbind, lapply(units, function(u) {
  s <- readRDS(unit_file(u$name))
  pick <- function(d, col) d[[col]][match(years, d$year)]
  data.frame(unit = u$name, year = years,
             site_pfpr = pick(s$prevalence, "pfpr"),
             tx_cov = pick(s$interventions$treatment$implementation, "tx_cov"),
             itn_use = pick(s$interventions$itn$use, "itn_use"),
             irs_cov = pick(s$interventions$irs$implementation, "irs_cov"),
             eir = s$eir$eir[s$eir$sp == "pf"])
}))
saveRDS(list(runs = res, site = site,
             versions = list(magenta = readLines("pkgdepends.txt")[2])),
        "site_summary.rds")
