#'Calculate correction factor matrix from bulk data.
#'
#'This function calculates the correction factor that borrows information from bulk data. If bulk data is available, it can be used to modify the prior setting in Bayesian hierarchical model.
#'
#' @param matrix bulk matrix
#' @param nei neighborhood size
#'
#' @return a matrix of bulk data correction factor
#' @export
#'
#' @examples
#' correctfac(m, nei=2) #calculate bulk data correction factor, m is the bulk matrix
correctfac <- function(type_bulk,nei){
  n = dim(type_bulk)[1]
  mean_bulk = sum(type_bulk[upper.tri(type_bulk)])/sum(upper.tri(type_bulk))
  correct = matrix(rep(0,n*n),nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m <- type_bulk[max(1,i-nei):min(n,i+nei), max(1,j-nei):min(n,j+nei)]
      c <- mean(m[!is.na(m)])
      correct[i,j] = c
    }
  }
  return(correct)
}