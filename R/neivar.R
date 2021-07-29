#' Calculate neighborhood standard deviation.
#' 
#' This function calculates the standard deviation of nonzero observations among single cells.
#'
#' @param single matrix with each column being the upper triangular of each single cell
#' @param nei neighborhood size
#' @param n matrix size
#'
#' @return standard deviation matrix with each position being the standard deviation in the neighborhood among all single cells.
#' @export
#'
#' @examples
#' neivar(m,nei=2) 
#' 
neivar <- function(single, nei, n){
  single_mat <- list()
  for (k in 1:dim(single)[2]) {
    m <- matrix(NA,n,n)
    m[upper.tri(m)] <- single[,k]
    single_mat[[k]] <- m
  }
  correct = matrix(0 ,nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      neighbor <- NULL
      for (l in 1:dim(single)[2]) {
        neighbor = c(neighbor,as.vector(single_mat[[l]][max(1,i-nei):min(n,i+nei), max(1,j-nei):min(n,j+nei)]))
      }
      c <- neighbor[!is.na(neighbor)]
      cc <- c[c>0]
      correct[i,j] = sd(cc)
    }
  }
  return(correct)
  }