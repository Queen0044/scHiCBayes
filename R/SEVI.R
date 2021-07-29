#' SEVI
#'
#' @param observed A vector of observed single cell.
#' @param expected A vector of expected single cell.
#' @param imputed A vector of imputed single cell.
#'
#' @return The scatterplot of expected versus imputed, with read dots being the observed zero pairs.
#' @export
#'
#' @examples
#' SEVI(observed=K562_T1_7k[,1], expected=K562_1_true[,1], imputed=T1_7k_imp[,1] )

SEVI <- function(observed, expected, imputed) {
   plot(expected, imputed, pch=20, cex=0.4, col="blue", xlab = "Expected", ylab = "Imputed",main = "SEVI")
   points(expected[observed==0], imputed[observed==0], pch=20, cex=0.4, col="red")
}

