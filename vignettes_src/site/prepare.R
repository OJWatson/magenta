## Cut the two units out of the Kenya site file (foresite calibration packet,
## calibrated_scaled_site.rds) so the cluster tasks only read small files.
## Run locally, from this folder, with the Malaria share mounted.

source("funs.R")
path <- file.path("/mnt/dide/Malaria/Pete/foresite-orderly/archive/calibration",
                  "20260612-145526-c809460a/calibrated_scaled_site.rds")
s <- readRDS(path)
dir.create("data", showWarnings = FALSE)
for (u in units) {
  k <- s$sites$name_1 == u$name_1 & s$sites$urban_rural == u$urban_rural
  saveRDS(site::subset_site(s, s$sites[k, ]), unit_file(u$name))
}
