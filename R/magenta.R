#' magenta: Individual-based simulation model of malaria epidemiology and genomics.
#'
#' @description
#' magenta is an individual-based simulation model of malaria epidemiology and genomics.
#' magenta extends the imperial malaria model by tracking the infection history of individuals.
#' With this additional genetic characteristics of the parasite can be assessed.
#'
#' @docType package
#' @name magenta
#'
#' @useDynLib magenta
#' @importFrom Rcpp sourceCpp
#' @importFrom stats runif qt sd approx qexp qlnorm rbinom rgeom rnbinom
#' @importFrom stats weighted.mean na.omit optim pgamma dexp time
#' @importFrom utils read.csv tail packageVersion
#' @importFrom rlang .data
#' @importFrom stats heatmap
#' @importFrom graphics par plot
#' @importFrom grDevices dev.off tiff hcl
#' @importFrom odin odin
#' @importFrom dde dopri
#' @keywords internal
"_PACKAGE"

globalVariables(c("admin_units_seasonal", "irs_2000_2015", "itn_2000_2015"))
