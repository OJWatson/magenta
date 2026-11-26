## Provision the cluster library (pkgdepends.txt) and the "site" environment.
## Run from a copy of this folder on the Malaria share.
library(hipercow)
hipercow_init(driver = "dide-windows", r_version = "4.5.3")
hipercow_environment_create(name = "site", sources = "funs.R", overwrite = TRUE)
hipercow_provision(environment = "site")
