#' @title Gini coefficient
#'
#' @description Takes a matrix and returns the Gini coefficient of each column.
#'
#' @param input A matrix where the Gini coefficient will be calculated on each column.
#' 
#' @param byrow Defaults to FALSE. Set to TRUE to calculate the Gini coefficient of each row.
#'
#' @return A vector of the Gini coefficients of each column.
#' 
#' @references Gini, C. (1912). Variabilita e mutabilita. Reprinted in Memorie di metodologica statistica (Ed. Pizetti E, Salvemini, T). Rome: Libreria Eredi Virgilio Veschi.
#'
#' @examples
#' gini(matrix(runif(20, 0, 1), nrow = 10, ncol = 2))
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
  output <- colSums(abs(sweep(cdf, 1, equality))) / (size - 1) # Divide by size - 1 because the last value in both equality and output will always be the same (1)

  # The output is a vector of Gini coefficients for each column (row if byrow = TRUE) of the input matrix
  return(output)
}
