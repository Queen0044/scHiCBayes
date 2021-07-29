#' Calculate the position in the vector
#' 
#' This function calculates the position of a pair in the vector when we transform the upper triangular matrix to a vector.
#'
#' @param x bin1 
#' @param y bin2
#'
#' @return position in the vector
#' @export
#'
#' @examples
#' matrow(2,3) #calculate the position of bin2 and bin3 in the vector.
matrow <- function(x,y) {   
  r <- x + (y-1)*(y-2)/2
  return(r)
}