
#' Converts the upper triangular matrix to a vector.
#'
#'This function converts the upper triangular part of a matrix to a vector
#'
#' @param  mat A matrix.
#'
#' @return A vector of the upper triangular of the matrix.
#' @export
#'
#' @examples
#' mattovec(m)
mattovec <- function(mat){
  vec <- mat[upper.tri(mat, diag = FALSE)]
  return(vec)
}