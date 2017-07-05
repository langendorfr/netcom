#' @title Gini coefficient
#'
#' @description Takes a matrix and returns the Gini coefficient of each row.
#'
#' @param input A row-stochastic matrix, meaning each row is a multinomial distribution.
#' 
#' @param byrow Defaults to FALSE. Set to TRUE to calculate the Gini coefficient of each row.
#'
#' @return A vector of the Gini coefficient of each row.
#'
#' @examples
#' gini(matrix(runif(25, 0, 1), nrow = 5, ncol = 5))
#' 
#' @export

gini <- function(input, byrow = FALSE)
{
  # Convert row-oriented input to columns 
  if (byrow == TRUE) {
    input <- t(input)
  }
  
  # Number of elements in each Gini coefficient
  size <- nrow(input)
  
  # Convert the input to a column-stochastic matrix (won't affect inputs already in that format)
  matrix <- sweep(input, 2, Matrix::colSums(input), FUN = "/")

  # Sort each column into ascending order (a CDF)
  sorted <- apply(matrix, 2, sort)

  cdf <- matrix(NA, nrow = size, ncol = ncol(input))
  cdf[1, ] <- sorted[1, ]
  for (i in 2:size) {
    cdf[i, ] = cdf[i - 1, ] + sorted[i, ]
  }

  # The Gini coefficient is relative to a perfectly equitable distribution
  equality <- seq(from = 1/size, to = 1, by = 1/size)
  output <- colSums(abs(sweep(cdf, 1, equality))) / size

  # The output is a vector of Gini coefficients for each column (row if byrow = TRUE) of the input matrix
  return(output)
}
