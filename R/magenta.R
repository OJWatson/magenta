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
#' @importFrom stats rlnorm heatmap runif qt sd approx qexp qlnorm rbinom rgeom
#' @importFrom stats rnbinom time weighted.mean na.omit optim pgamma dexp
#' @importFrom utils adist read.csv tail packageVersion
#' @importFrom ggplot2 ggplot
#' @importFrom graphics par plot
#' @importFrom rlang .data
#' @importFrom grDevices dev.off tiff dev.new hcl
#' @importFrom ggplot2 geom_line scale_color_manual scale_shape_manual theme
#' @importFrom ggplot2 guides guide_legend scale_color_discrete ylab
#' @importFrom ggplot2 element_blank
#' @importFrom odin odin
#' @keywords internal
"_PACKAGE"

globalVariables(c("admin_units_seasonal", "irs_2000_2015", "itn_2000_2015", ""))
