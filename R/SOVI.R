#' SOVI
#'
#' Scatterplot of observed versus imputed at nonzeros.
#' 
#' @param observed A vector of observed single cell.
#' @param imputed A vector of imputed single cell.
#'
#' @return The scatterplot of observed versus imputed.
#' @export
#'
#' @examples 
#' data("GSE117874_imp")
#' SOVI(observed=GSE117874_chr1_wo_diag[,1], imputed = GSE117874_imp[,1])

SOVI <- function(observed, imputed) {
  plot(observed[observed>0], imputed[observed>0], pch=20, cex=0.4, col="blue", xlab = "Observed", ylab = "Imputed",main = "SOVI")
}